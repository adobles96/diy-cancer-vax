"""Run: .venv/bin/python scripts/cohort_recall.py --rescue-unknown-hla"""

import argparse
import datetime
import html
import json
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from neovax.cohort import DEFAULT_RECALL_K, METHODS, audit_components, rank_targets, recall_at_k
from neovax.config import Config
from neovax.data import actual_vaccine, load_reference, patient_data, verify_sources
from neovax.deepimmuno import predict_immunogenicity
from neovax.iedb import METHOD, BindingCache
from neovax.io import digest, write_csv, write_json
from neovax.peptides import enumerate_epitopes, reconstruct
from neovax.reference_filter import filter_reference_peptides, reference_filter_source

LABELS = {"netmhcpan": "NetMHCpan", "deepimmuno": "DeepImmuno", "mean_ensemble": "Mean ensemble"}
METHODOLOGY = f"""Top {DEFAULT_RECALL_K} distinct genomic mutation targets from independently reconstructed missense candidates,
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


def render(output, metrics, pooled, metadata):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True, layout="constrained")
    colors = ["#3266a8", "#bb6233", "#38816a"]
    for ax, endpoint in zip(axes, ["vaccine", "immunogenic"]):
        for index, method in enumerate(METHODS):
            values = [
                next(
                    r["recall"]
                    for r in metrics
                    if r["patient"] == p and r["method"] == method and r["endpoint"] == endpoint
                )
                for p in range(1, 7)
            ]
            ax.bar(
                [p + (index - 1) * 0.24 for p in range(1, 7)],
                [100 * v for v in values],
                width=0.23,
                color=colors[index],
                label=LABELS[method],
            )
        ax.set(
            title=f"{'Vaccine components' if endpoint == 'vaccine' else 'Any CD4/CD8 peptide response'}",
            xlabel="Patient",
            xticks=range(1, 7),
            ylim=(0, 105),
        )
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.15)
    axes[0].set_ylabel(f"Recall at {DEFAULT_RECALL_K} targets (%) · full denominator")
    axes[1].legend(frameon=False, loc="upper right")
    fig.suptitle(f"NeoVax: top-{DEFAULT_RECALL_K} target recovery", fontsize=16)
    fig.savefig(output / f"recall-at-{DEFAULT_RECALL_K}.png", dpi=180)
    fig.savefig(output / f"recall-at-{DEFAULT_RECALL_K}.svg")
    plt.close(fig)
    tables = []
    for endpoint in ["vaccine", "immunogenic"]:
        body = []
        for patient in range(1, 7):
            group = [r for r in metrics if r["patient"] == patient and r["endpoint"] == endpoint]
            cells = "".join(
                f"<td>{r['hits']}/{r['denominator']} ({r['recall']:.1%})<small>Within reconstructable: {r['hits']}/{r['reconstructable_denominator']} ({r['reconstructable_recall']:.1%})</small></td>"
                for r in group
            )
            body.append(
                f"<tr><th>P{patient}</th>{cells}<td>{group[0]['reconstructable_denominator']}/{group[0]['denominator']}</td></tr>"
            )
        group = [r for r in pooled if r["endpoint"] == endpoint]
        cells = "".join(
            f"<td>{r['hits']}/{r['denominator']} ({r['micro_recall']:.1%})<small>Patient mean: {r['macro_recall']:.1%}</small></td>"
            for r in group
        )
        body.append(
            f"<tr><th>Pooled</th>{cells}<td>{group[0]['reconstructable_denominator']}/{group[0]['denominator']}</td></tr>"
        )
        tables.append(
            f"<h2>{'Vaccine component recall' if endpoint == 'vaccine' else 'Immunogenic component recall'}</h2><table><thead><tr><th>Patient</th><th>NetMHCpan</th><th>DeepImmuno</th><th>Mean ensemble</th><th>Reconstructable</th></tr></thead><tbody>{''.join(body)}</tbody></table>"
        )
    rescue = {f"P{p['patient']}": p["deepimmuno"]["rescued_alleles"] for p in metadata}
    filter_rows = "".join(
        f"<tr><th>P{p['patient']}</th><td>{p['reference_filter']['unique_peptides_before']}</td><td>{p['reference_filter']['unique_peptides_removed']}</td><td>{p['reference_filter']['unique_peptides_after']}</td><td>{len(p['reference_filter']['mutations_fully_removed'])}</td></tr>"
        for p in metadata
    )
    page = f"""<!doctype html><html lang="en"><meta charset="utf-8"><title>NeoVax cohort recall</title>
<style>body{{font:16px/1.5 system-ui;max-width:1200px;margin:40px auto;padding:0 24px;color:#223042}}table{{border-collapse:collapse;width:100%;margin:20px 0 36px}}td,th{{text-align:left;padding:12px;border-bottom:1px solid #ddd}}small{{display:block;color:#637080;font-size:12px}}pre{{white-space:pre-wrap;background:#f2f5f8;padding:20px}}img{{width:100%}}h1,h2{{color:#163a5f}}</style>
<h1>NeoVax cohort recall at {DEFAULT_RECALL_K}</h1><p>Six patients · {DEFAULT_RECALL_K} distinct mutation targets per method · NetMHCpan BA 4.1 · DeepImmuno with disclosed HLA rescue · equal-weight rank ensemble</p>
<img src="recall-at-{DEFAULT_RECALL_K}.svg" alt="Recall by patient and model">{"".join(tables)}
<h2>Reference-proteome filter</h2><p>Exact matches to any sequence in the pinned hg19 reference protein set are excluded before scoring and ranking. These counts refer to distinct short peptide sequences per patient; multiple windows or transcripts do not multiply the count.</p>
<table><thead><tr><th>Patient</th><th>Before</th><th>Removed</th><th>Retained</th><th>Mutations losing all candidates</th></tr></thead><tbody>{filter_rows}</tbody></table>
<p><a href="reference-filtered-epitopes.csv">Excluded peptides with matching reference transcripts and positions</a> · <a href="reference-filter-summary.json">Filter counts and reference provenance</a>. This is a reference-proteome approximation, not a patient-specific germline screen.</p>
<h2>Definitions and limits</h2><pre>{html.escape(METHODOLOGY)}</pre>
<h2>Actual HLA → substituted DeepImmuno HLA</h2><pre>{html.escape(json.dumps(rescue, indent=2))}</pre>
<p>Binding predictions always use the actual HLA. A rescued DeepImmuno score is an approximation.</p>
<p>Downloads: <a href="metrics.csv">Metrics</a> · <a href="vaccine-component-audit.csv">Component evidence and ranks</a> · <a href="target-rankings.csv">All target rankings</a> · <a href="recall-curves.csv">Recall at 10–100</a> · <a href="analysis.json">Full provenance</a></p></html>"""
    (output / "report.html").write_text(page)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rescue-unknown-hla", action="store_true")
    parser.add_argument("--cached-only", action="store_true")
    parser.add_argument("--output-dir", default="data/outputs/cohort-recall")
    args = parser.parse_args()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    metrics, curves, all_targets, all_components, metadata = [], [], [], [], []
    all_filtered, filter_summaries = [], []
    for patient in range(1, 7):
        print(f"Patient {patient}: reconstructing independent candidate targets", flush=True)
        config = Config(
            patient=patient,
            threshold=None,
            cached_only=args.cached_only,
            deepimmuno_manifest="data/outputs/models/deepimmuno/manifest.json",
        )
        variants, alleles, published, assays = patient_data(config.data_dir, patient)
        reference = load_reference(config.data_dir)
        proteins, skipped = reconstruct(variants, reference)
        epitopes = enumerate_epitopes(proteins, config.lengths)
        epitopes, filtered, filter_summary = filter_reference_peptides(epitopes, reference)
        filter_summary.update(patient=patient, source=reference_filter_source())
        all_filtered.extend(filtered)
        filter_summaries.append(filter_summary)
        print(
            f"Reference filter: {filter_summary['unique_peptides_removed']} distinct peptides removed",
            flush=True,
        )
        cache = BindingCache(config.cache_dir)
        try:
            binding = cache.lookup(alleles)
        finally:
            cache.close()
        immuno, provenance = predict_immunogenicity(
            epitopes,
            alleles,
            config,
            progress=lambda message: print(message, flush=True),
            allow_hla_rescue=args.rescue_unknown_hla,
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
                "skipped_variants": skipped,
                "reference_filter": filter_summary,
            }
        )
        write_json(
            output / f"patient-{patient}.json", {"metrics": patient_metrics, "provenance": metadata[-1]}
        )
        print(f"Patient {patient}: complete; {len(targets)} ranked targets", flush=True)
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
        "methodology": METHODOLOGY,
        "recall_k": DEFAULT_RECALL_K,
        "binding_model": METHOD,
        "allow_hla_rescue": args.rescue_unknown_hla,
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
    render(output, metrics, pooled, metadata)
    write_json(
        output / "artifact-hashes.json",
        {p.name: digest(p) for p in output.iterdir() if p.is_file() and p.name != "artifact-hashes.json"},
    )
    print(json.dumps(pooled, indent=2), flush=True)
    print(f"Report: {output / 'report.html'}", flush=True)


if __name__ == "__main__":
    main()
