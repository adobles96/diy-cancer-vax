"""Reproduce the six-patient retrospective blog benchmark from pinned source data."""

import csv
import datetime
import html
import json
import statistics
from dataclasses import replace
from pathlib import Path

from neovax.assay_auc import assay_labels, roc_metrics
from neovax.blog_figures import render_auc, render_recall
from neovax.cohort import DEFAULT_RECALL_K, METHODS, audit_components, rank_targets, recall_at_k
from neovax.config import Config
from neovax.data import actual_vaccine, load_reference, patient_data, verify_sources
from neovax.deepimmuno import predict_immunogenicity
from neovax.iedb import METHOD, BindingCache, Cancelled, InferenceError, predict
from neovax.io import digest, write_csv, write_json
from neovax.peptides import enumerate_epitopes, reconstruct
from neovax.reference_filter import filter_reference_peptides, reference_filter_source

LABELS = {
    "neovax": "NeoVax",
    "netmhcpan": "NetMHCpan",
    "deepimmuno": "DeepImmuno",
    "mean_ensemble": "Ensemble",
}

RECALL_METHODOLOGY = f"""Top {DEFAULT_RECALL_K} distinct genomic mutation targets from independently reconstructed missense candidates,
representing a budget of one vaccine component per mutation, not 20 overlapping short epitopes.
using all mutation-containing 9/10-mers and the patient's HLA-A/B alleles. No binding threshold.
Short candidates matching any protein in the pinned hg19 reference proteome are excluded before scoring
and ranking. Matching transcript IDs and positions are recorded in reference-filtered-epitopes.csv.
NetMHCpan BA 4.1 target score = minimum IC50 over peptide–HLA pairs. DeepImmuno target score = maximum
immunogenicity score over those pairs. Each model's target scores are converted to fractional ranks
within the patient, averaging ties (0 best, 1 worst); ensemble cost is their arithmetic mean.
The two models may select different driving peptide–HLA pairs for a target.
Final ranking ties use genomic mutation ID, without consulting vaccine membership or response labels.
Vaccine recall counts published IMP components whose mutation is selected, divided by all published
IMP components for that patient. This measures target recovery, not exact long-peptide sequence recovery.
Immunogenic recall restricts that denominator to components with any positive week-16 peptide assay:
CD4 OR CD8, ex vivo OR expanded cells. No requirement for mutant selectivity, minigene response,
tumor recognition or killing. Later observations are excluded to avoid mixing in subsequent therapy.
Unreconstructable components remain in the full denominator; conditional recall and coverage ceilings
are also reported. Untested candidates are not assigned negative labels.
These class-I models do not directly predict CD4 responses; this broad endpoint is descriptive.
The source tables contain 98 distinct patient/IMP components, whereas the paper states 97 unique
neoantigens. We retain all table components and disclose this unresolved denominator discrepancy.
This is an exploratory retrospective benchmark with incompletely audited model training overlap."""


AUC_METHODOLOGY = """Evaluation is restricted to actual vaccinated mutations with independently reconstructed
and scored candidates and usable week-16 peptide-response labels. Each genomic mutation counts once.
Positive = any positive CD4 or CD8 peptide response, ex vivo or after expansion, from any associated IMP.
Specifically: IFN-γ ELISpot against peptide-pulsed autologous antigen-presenting cells, using CD8
cells after pre-stimulation, CD4 responses directly ex vivo, or CD4 cells after pre-stimulation.
We use the authors' positive/negative calls; we do not impose a new threshold on spot counts.
Negative = both CD4 and CD8 have measurements for each associated IMP, and all eligible measurements
are explicitly negative. Missing, not-determined, and ambiguous observations do not become negatives.
Minigene, tumor, and later-timepoint observations do not enter the label. This is observed response
under the reported assays, not proof of biological non-immunogenicity.
Scores are unchanged from the previous cohort analysis: minimum NetMHCpan IC50 or maximum DeepImmuno
score over all reconstructed mutation-containing 9/10-mer–HLA-A/B pairs for each mutation.
The candidate set inherits the reference-proteome exact-match filter recorded in the input provenance.
The ensemble averages the two within-patient fractional ranks over the full discovery candidate set;
it is not refitted or renormalized on the vaccine evaluation set. No top-K or binding cutoff is used.
These are mutation-level scores; the best short peptide need not be the experimentally tested peptide.
AUC is the fraction of positive–negative target pairs correctly ordered, giving score ties half credit.
ROC curves group equal scores; genomic-ID tie breaks from shortlist generation are not used for AUC.
Across patients, report the equally weighted mean AUC and a pair-weighted within-patient AUC.
The latter pools concordance counts but never compares a target from one patient with another patient's
target. It is not a conventional AUC computed by concatenating all patient scores.
The estimand is discrimination conditional on vaccination and successful reconstruction. Historical
binding-based selection still limits generalization to all mutations. These class-I models do not
directly predict CD4 responses. HLA rescue and training-overlap limitations from the prior run persist.
With only 9–15 evaluated mutations per patient, these descriptive estimates have substantial uncertainty."""


def artifact_hashes(output):
    write_json(
        output / "artifact-hashes.json",
        {p.name: digest(p) for p in output.iterdir() if p.is_file() and p.name != "artifact-hashes.json"},
    )


def predict_retained_binding(cache, epitopes, alleles, config, progress, cancel):
    """Score only retained candidates against actual HLA, one serial request stream."""
    binding, statistics = {}, []
    for length in (9, 10):
        if cancel():
            raise Cancelled("Blog reproduction cancelled before binding inference")
        sequences = sorted({r["peptide"] for r in epitopes if len(r["peptide"]) == length})
        if not sequences:
            continue
        scores, stats, status, errors = predict(
            cache,
            sequences,
            alleles,
            (length,),
            batch_size=config.batch_size,
            cached_only=config.cached_only,
            progress=progress,
            cancel=cancel,
        )
        statistics.append({"length": length, **stats})
        if status == "cancelled":
            raise Cancelled("Blog reproduction cancelled during binding inference")
        if status != "complete":
            detail = "; ".join(errors) or f"{stats['remaining_requests']} binding requests are not cached"
            raise InferenceError(f"Cannot rank an incomplete blog cohort: {detail}")
        binding.update(
            {(sequence, allele): scores[(sequence, allele)] for sequence in sequences for allele in alleles}
        )
    return binding, statistics


def write_run_report(output, title, figure, methodology):
    downloads = " · ".join(
        f'<a href="{html.escape(p.name)}">{html.escape(p.name)}</a>'
        for p in sorted(output.iterdir())
        if p.suffix in {".csv", ".json", ".png", ".svg", ".pdf"}
    )
    (output / "report.html").write_text(
        f'<!doctype html><html lang="en"><meta charset="utf-8"><title>{title}</title>'
        f"<style>body{{font:16px/1.6 system-ui;max-width:1150px;margin:40px auto;padding:0 24px}}"
        f"img{{width:100%}}pre{{white-space:pre-wrap}}</style><h1>{title}</h1>"
        f'<img src="{figure}.svg" alt="{title}"><p>{downloads}</p>'
        f"<h2>Definitions and limits</h2><pre>{html.escape(methodology)}</pre></html>"
    )


def run_cohort_recall(
    base_config, output_dir, progress=lambda message: None, cancel=lambda: False, *, allow_hla_rescue=False
):
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    metrics, curves, all_targets, all_components, metadata = [], [], [], [], []
    all_filtered, filter_summaries = [], []
    for patient in range(1, 7):
        if cancel():
            raise Cancelled("Blog reproduction cancelled between patients")
        progress(f"Patient {patient}/6: reconstructing independent candidate targets")
        config = replace(
            base_config,
            patient=patient,
            lengths=(9, 10),
            threshold=None,
            budget=20,
            selection_unit="mutation",
            alpha=0.5,
            rescue_unknown_hla=allow_hla_rescue,
            output_dir=str(output),
        )
        variants, alleles, published, assays = patient_data(config.data_dir, patient)
        reference = load_reference(config.data_dir)
        proteins, skipped = reconstruct(variants, reference)
        epitopes = enumerate_epitopes(proteins, config.lengths)
        epitopes, filtered, filter_summary = filter_reference_peptides(epitopes, reference)
        filter_summary.update(patient=patient, source=reference_filter_source())
        all_filtered.extend(filtered)
        filter_summaries.append(filter_summary)
        progress(
            f"Patient {patient}/6: {filter_summary['unique_peptides_removed']} reference peptides removed"
        )
        cache = BindingCache(config.cache_dir)
        try:
            binding, binding_statistics = predict_retained_binding(
                cache, epitopes, alleles, config, progress, cancel
            )
        finally:
            cache.close()
        immuno, provenance = predict_immunogenicity(
            epitopes,
            alleles,
            config,
            progress=progress,
            cancel=cancel,
            allow_hla_rescue=allow_hla_rescue,
        )
        targets = rank_targets(epitopes, alleles, binding, immuno, provenance["allele_representations"])
        # Labels enter only after all target scores and rankings are fixed.
        components = audit_components(actual_vaccine(published, variants, proteins), proteins, assays)
        by_target = {r["mutation_id"]: r for r in targets}
        for component in components:
            row = by_target.get(component["mutation_id"], {})
            component["rankable"] = bool(row)
            for method in METHODS:
                component[f"{method}_rank"] = row.get(f"{method}_rank")
        patient_metrics = [{"patient": patient, **r} for r in recall_at_k(targets, components)]
        metrics.extend(patient_metrics)
        for k in [10, 20, 30, 40, 50, 75, 100]:
            curves.extend({"patient": patient, **r} for r in recall_at_k(targets, components, k))
        all_targets.extend({"patient": patient, **r} for r in targets)
        all_components.extend(components)
        metadata.append(
            {
                "patient": patient,
                "config": config.to_dict(),
                "alleles": alleles,
                "deepimmuno": provenance,
                "binding": binding_statistics,
                "skipped_variants": skipped,
                "reference_filter": filter_summary,
            }
        )
        write_json(
            output / f"patient-{patient}.json", {"metrics": patient_metrics, "provenance": metadata[-1]}
        )
        progress(f"Patient {patient}/6: complete; {len(targets)} ranked targets")
    pooled = []
    for endpoint in ["vaccine", "immunogenic"]:
        for method in METHODS:
            group = [r for r in metrics if r["endpoint"] == endpoint and r["method"] == method]
            hits, denominator = sum(r["hits"] for r in group), sum(r["denominator"] for r in group)
            available = sum(r["reconstructable_denominator"] for r in group)
            pooled.append(
                {
                    "method": method,
                    "endpoint": endpoint,
                    "hits": hits,
                    "denominator": denominator,
                    "micro_recall": hits / denominator,
                    "macro_recall": statistics.mean(r["recall"] for r in group),
                    "reconstructable_denominator": available,
                    "reconstructable_recall": hits / available,
                }
            )
    result = {
        "created_at": datetime.datetime.now(datetime.UTC).isoformat(),
        "methodology": RECALL_METHODOLOGY,
        "recall_k": DEFAULT_RECALL_K,
        "binding_model": METHOD,
        "allow_hla_rescue": allow_hla_rescue,
        "sources": verify_sources(Path(config.data_dir)),
        "patients": metadata,
        "metrics": metrics,
        "pooled": pooled,
        "reference_filter": filter_summaries,
    }
    for name, rows in [
        ("metrics", metrics),
        ("pooled", pooled),
        ("target-rankings", all_targets),
        ("vaccine-component-audit", all_components),
        ("recall-curves", curves),
        ("reference-filtered-epitopes", all_filtered),
    ]:
        write_csv(output / f"{name}.csv", rows)
    write_json(output / "analysis.json", result)
    write_json(output / "reference-filter-summary.json", filter_summaries)
    render_recall(output, metrics, pooled)
    write_run_report(output, "Cohort recall at 20", "vaccine-recall-at-20", RECALL_METHODOLOGY)
    artifact_hashes(output)
    return result


def run_cohort_auc(input_dir, output_dir, data_dir=None, progress=lambda message: None, cancel=lambda: False):
    source, output = Path(input_dir).resolve(), Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    hashes = json.loads((source / "artifact-hashes.json").read_text())
    inputs = ["analysis.json", "target-rankings.csv", "vaccine-component-audit.csv"]
    for name in inputs:
        if digest(source / name) != hashes[name]:
            raise ValueError(f"Input file changed since the cohort run: {name}")
    provenance = json.loads((source / "analysis.json").read_text())
    with (source / "target-rankings.csv").open() as handle:
        targets = list(csv.DictReader(handle))
    with (source / "vaccine-component-audit.csv").open() as handle:
        components = list(csv.DictReader(handle))
    for target in targets:
        for method in METHODS:
            target[f"{method}_cost"] = float(target[f"{method}_cost"])
    metrics, all_labeled, all_excluded = [], [], []
    for patient in range(1, 7):
        config = next(p["config"] for p in provenance["patients"] if p["patient"] == patient)
        if cancel():
            raise Cancelled("Blog reproduction cancelled during assay evaluation")
        progress(f"Patient {patient}/6: evaluating week-16 peptide recognition")
        _, _, _, observations = patient_data(data_dir or config["data_dir"], patient)
        labeled, excluded = assay_labels(
            [r for r in components if int(r["patient"]) == patient],
            [r for r in targets if int(r["patient"]) == patient],
            observations,
        )
        all_labeled.extend(labeled)
        all_excluded.extend(excluded)
        for method in METHODS:
            metrics.append(
                {
                    "patient": patient,
                    "method": method,
                    **roc_metrics([r["label"] for r in labeled], [r[f"{method}_cost"] for r in labeled]),
                }
            )
    aggregate = []
    for method in METHODS:
        group = [r for r in metrics if r["method"] == method and r["roc_auc"] is not None]
        pairs = sum(r["positive_negative_pairs"] for r in group)
        credits = sum(r["correctly_ordered_pairs"] + 0.5 * r["tied_pairs"] for r in group)
        aggregate.append(
            {
                "method": method,
                "patients_with_defined_auc": len(group),
                "macro_auc": statistics.mean(r["roc_auc"] for r in group) if group else None,
                "within_patient_pair_weighted_auc": credits / pairs if pairs else None,
                "correct_pair_credits": credits,
                "within_patient_pairs": pairs,
            }
        )
    result = {
        "created_at": datetime.datetime.now(datetime.UTC).isoformat(),
        "methodology": AUC_METHODOLOGY,
        "inputs": {str(source / name): hashes[name] for name in inputs},
        "score_provenance": provenance,
        "labeled_targets": len(all_labeled),
        "metrics": metrics,
        "aggregate": aggregate,
        "excluded": all_excluded,
    }
    write_json(output / "analysis.json", result)
    write_csv(output / "metrics.csv", [{k: v for k, v in r.items() if k != "roc_curve"} for r in metrics])
    write_csv(output / "aggregate.csv", aggregate)
    write_csv(output / "labeled-targets.csv", all_labeled)
    write_csv(output / "excluded.csv", all_excluded)
    render_auc(output, metrics, aggregate)
    write_run_report(output, "Peptide recognition ROC AUC", "peptide-recognition-roc-auc", AUC_METHODOLOGY)
    artifact_hashes(output)
    return result


TUMOR_METHODOLOGY = """Any positive author-called week-16 CD4 or CD8 autologous tumor recognition assay,
after expansion. Tested requires at least one explicit positive or negative result. Missing,
not-determined and ambiguous observations are not negatives. Components are counted once;
top-20 model overlap is matched by genomic mutation. The denominator is tested components
within each method's recovered published vaccine set, not all 20 selected candidates.
Tumor recognition does not demonstrate tumor-cell killing or clinical benefit."""


def recognition_tables(components, observations):
    """Summarize published components and model overlap without labeling untested candidates."""
    evidence = []
    for component in components:
        tumor = [
            o
            for o in observations
            if o["patient"] == component["patient"]
            and o["imp_id"] == component["imp_id"]
            and o["timepoint"] == "week 16"
            and o["assay_type"] == "tumor"
            and o["cell_type"] in {"CD4", "CD8"}
        ]
        selected = ["neovax"] + [
            method
            for method in METHODS
            if component.get(f"{method}_rank") and int(component[f"{method}_rank"]) <= DEFAULT_RECALL_K
        ]
        evidence.append(
            {
                **component,
                "selected_by": selected,
                "tumor_tested": any(o["result"] in {"positive", "negative"} for o in tumor),
                "tumor_positive": any(o["result"] == "positive" for o in tumor),
                "tumor_observations": tumor,
            }
        )
    peptide, tumor = [], []
    for patient in [*range(1, 7), "Total"]:
        patients = [r for r in evidence if patient == "Total" or r["patient"] == patient]
        for method in ("neovax", *METHODS):
            included = [r for r in patients if method in r["selected_by"]]
            tested = [r for r in included if r["tumor_tested"]]
            peptide.append(
                {
                    "patient": patient,
                    "method": method,
                    "positive": sum(r["immunogenic"] for r in included),
                    "included": len(included),
                    "fraction": sum(r["immunogenic"] for r in included) / len(included) if included else None,
                }
            )
            tumor.append(
                {
                    "patient": patient,
                    "method": method,
                    "positive": sum(r["tumor_positive"] for r in tested),
                    "tested": len(tested),
                    "included": len(included),
                    "untested": len(included) - len(tested),
                    "fraction": sum(r["tumor_positive"] for r in tested) / len(tested) if tested else None,
                }
            )
    return peptide, tumor, evidence


def count_table(rows, denominator):
    table = []
    for patient in [*range(1, 7), "Total"]:
        row = {"Patient": patient}
        for method in ("neovax", *METHODS):
            item = next(r for r in rows if r["patient"] == patient and r["method"] == method)
            row[LABELS[method]] = f"{item['positive']}/{item[denominator]}"
        table.append(row)
    return table


def html_table(rows):
    fields = list(rows[0])
    head = "".join(f"<th>{html.escape(str(field))}</th>" for field in fields)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row[field]))}</td>" for field in fields) + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def reproduce_blog(
    config: Config, output_dir: Path, progress=lambda message: None, cancel=lambda: False
) -> dict:
    """Produce the blog's charts and audited tables without relying on prior run outputs."""
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    progress("Blog reproduction: reconstructing, filtering and scoring all six patients")
    recall = run_cohort_recall(config, output / "cohort-recall", progress, cancel, allow_hla_rescue=True)
    auc = run_cohort_auc(output / "cohort-recall", output / "cohort-auc", config.data_dir, progress, cancel)
    if cancel():
        raise Cancelled("Blog reproduction cancelled before rendering")
    progress("Building peptide recognition and tested-only tumor recognition tables")
    with (output / "cohort-recall" / "vaccine-component-audit.csv").open() as handle:
        components = list(csv.DictReader(handle))
    for row in components:
        row["patient"] = int(row["patient"])
        row["immunogenic"] = row["immunogenic"] == "True"
    observations = [o for patient in range(1, 7) for o in patient_data(config.data_dir, patient)[3]]
    peptide, tumor, evidence = recognition_tables(components, observations)
    peptide_table, tumor_table = count_table(peptide, "included"), count_table(tumor, "tested")
    with (output / "cohort-auc" / "labeled-targets.csv").open() as handle:
        labeled = list(csv.DictReader(handle))
    auc_table = []
    for patient in range(1, 7):
        group = [r for r in auc["metrics"] if r["patient"] == patient]
        auc_table.append(
            {
                "Patient": patient,
                "Tested targets": group[0]["n"],
                "Positive": group[0]["positives"],
                "Negative": group[0]["negatives"],
                **{LABELS[r["method"]]: r["roc_auc"] for r in group},
            }
        )
    positives = sum(r["Positive"] for r in auc_table)
    negatives = sum(r["Negative"] for r in auc_table)
    for key, label in [
        ("macro_auc", "Patient mean"),
        ("within_patient_pair_weighted_auc", "Within-patient pair-weighted"),
    ]:
        auc_table.append(
            {
                "Patient": label,
                "Tested targets": len(labeled),
                "Positive": positives,
                "Negative": negatives,
                **{LABELS[r["method"]]: r[key] for r in auc["aggregate"]},
            }
        )
    for filename, rows in [
        ("peptide-recognition-table", peptide_table),
        ("peptide-recognition-counts", peptide),
        ("tumor-recognition-tested-table", tumor_table),
        ("tumor-recognition-counts", tumor),
        ("roc-auc-table", auc_table),
    ]:
        write_csv(output / f"{filename}.csv", rows)
    write_json(
        output / "tumor-recognition-audit.json", {"definition": TUMOR_METHODOLOGY, "components": evidence}
    )
    progress("Rendering charts and the results report")
    render_recall(output, recall["metrics"], recall["pooled"])
    render_auc(output, auc["metrics"], auc["aggregate"])
    rescue = {f"P{p['patient']}": p["deepimmuno"]["rescued_alleles"] for p in recall["patients"]}
    summary = {
        "status": "complete",
        "report_path": str(output / "report.html"),
        "output_dir": str(output),
        "recall_k": DEFAULT_RECALL_K,
        "vaccine_recall": [r for r in recall["pooled"] if r["endpoint"] == "vaccine"],
        "peptide_recognition": [r for r in peptide if r["patient"] == "Total"],
        "tumor_recognition": [r for r in tumor if r["patient"] == "Total"],
        "auc": auc["aggregate"],
        "labeled_targets": len(labeled),
        "positive_labels": positives,
        "negative_labels": negatives,
        "allow_hla_rescue": True,
        "rescued_alleles": rescue,
        "binding_model": METHOD,
        "deepimmuno_model": recall["patients"][0]["deepimmuno"]["model"],
    }
    write_json(output / "summary.json", summary)
    downloads = " · ".join(
        f'<a href="{p.name}">{p.name}</a>'
        for p in sorted(output.iterdir())
        if p.suffix in {".csv", ".json", ".png", ".svg", ".pdf"}
    )
    formatted_auc = [
        {key: f"{value:.6f}" if isinstance(value, float) else value for key, value in row.items()}
        for row in auc_table
    ]
    page = f"""<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>DIY Cancer Vaccine — blog results</title>
<style>body{{font:16px/1.6 system-ui;max-width:1280px;margin:40px auto;padding:0 24px;background:#FAF8F3;color:#23332F}}
img{{width:100%;height:auto}}table{{border-collapse:collapse;width:100%;margin:20px 0 36px;background:white}}
td,th{{text-align:left;padding:10px;border-bottom:1px solid #D5DED6}}pre{{white-space:pre-wrap;padding:18px;background:#E8EEE9}}
a{{color:#176E60}}.table{{overflow:auto}}</style>
<h1>DIY Cancer Vaccine: blog results</h1>
<p>Six patients · top 20 distinct mutation targets per patient · actual patient HLA-A/B · no binding cutoff.</p>
<p>These results recover published vaccine targets by mutation. Peptide-response and tumor-response estimates
are restricted to experimentally observed vaccine components; untested new predictions have no response label.</p>
<h2>Vaccine target recall@20</h2><img src="vaccine-recall-at-20.svg" alt="Recall@20 for six patients and overall">
<p>Overall recovery: {", ".join(f"{LABELS[r['method']]} {r['hits']}/{r['denominator']}" for r in summary["vaccine_recall"])}.
All 98 patient/IMP source-table components remain in the denominator, including unreconstructable components.
The paper reports 97 unique neoantigens; this source-table denominator discrepancy remains unresolved.</p>
<h2>Peptide recognition</h2><p>Any author-called positive week-16 CD4 or CD8 peptide response, directly ex vivo
or after expansion. Numerator = responding components; denominator = published vaccine components recovered
by each method, or all published components for NeoVax. A zero denominator has an undefined rate.</p>
<div class="table">{html_table(peptide_table)}</div>
<h2>Peptide recognition ROC AUC</h2><img src="peptide-recognition-roc-auc.svg" alt="Per-patient AUC, patient mean, and within-patient pair-weighted aggregate">
<p>{len(labeled)} reconstructed vaccine mutations: {positives} positive, {negatives} negative. Negatives require
explicitly negative CD4 and CD8 measurements. AUC uses no top-20 cutoff. Patient mean and within-patient
pair-weighted AUC never compare different patients.</p>
<div class="table">{html_table(formatted_auc)}</div>
<h2>Tumor recognition: tested components only</h2><p>Numerator = components with any positive week-16 CD4/CD8
autologous tumor assay. Denominator = components with an explicit positive or negative tumor assay in that
method's recovered vaccine set. Untested components are excluded from this denominator.</p>
<div class="table">{html_table(tumor_table)}</div>
<h2>HLA substitutions and model provenance</h2><p>Binding uses every actual patient HLA allele. DeepImmuno uses
the pinned upstream numerical HLA rescue where an exact allele is unavailable. Rescued scores are approximations.</p>
<pre>{html.escape(json.dumps(rescue, indent=2))}</pre>
<p>Model versions, input hashes, reference filtering, target ranks, assay evidence and training-overlap audits
are recorded in <a href="cohort-recall/analysis.json">recall provenance</a> and
<a href="cohort-auc/analysis.json">AUC provenance</a>.</p>
<h2>Downloads</h2><p>{downloads}</p>
<p><a href="cohort-recall/report.html">Full recall data and methodology</a> ·
<a href="cohort-auc/report.html">Full AUC data and methodology</a></p>
<h2>Definitions and limits</h2><pre>{html.escape(RECALL_METHODOLOGY)}</pre>
<pre>{html.escape(AUC_METHODOLOGY)}</pre><pre>{html.escape(TUMOR_METHODOLOGY)}</pre></html>"""
    (output / "report.html").write_text(page)
    artifact_hashes(output)
    progress(f"Blog reproduction complete: {output / 'report.html'}")
    return summary
