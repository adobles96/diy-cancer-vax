"""Compare BigMHC IM with the saved six-patient cohort without changing defaults."""

import argparse
import csv
import datetime
import html
import json
import statistics
from pathlib import Path

from neovax.assay_auc import assay_labels, roc_metrics
from neovax.bigmhc import exact_alleles, predict_immunogenicity
from neovax.cohort import METHODS, audit_components, rank_targets
from neovax.config import Config
from neovax.data import actual_vaccine, load_reference, patient_data, verify_sources
from neovax.deepimmuno import predict_immunogenicity as predict_deepimmuno
from neovax.iedb import BindingCache
from neovax.io import digest, write_csv, write_json
from neovax.peptides import enumerate_epitopes, reconstruct
from neovax.reference_filter import filter_reference_peptides

PILOT_METHODS = (*METHODS, "bigmhc", "bigmhc_ensemble")
LABELS = {
    "netmhcpan": "NetMHCpan",
    "deepimmuno": "DeepImmuno",
    "mean_ensemble": "NetMHCpan + DeepImmuno",
    "bigmhc": "BigMHC IM",
    "bigmhc_ensemble": "NetMHCpan + BigMHC IM",
}
METHODOLOGY = """All six patients; same independently reconstructed, reference-proteome-filtered
mutation-containing 9/10-mers and actual HLA-A/B as the existing cohort baseline. No binding cutoff.
Every regenerated baseline target, pair count, raw score, cost and rank must match the saved cohort.
BigMHC IM is the unmodified official seven-model mean of sigmoid outputs, not BigMHC EL.
Exact patient HLA pseudosequences are required; no fuzzy substitution is permitted.
Target score = maximum IM score across peptide-HLA pairs; NetMHCpan = minimum IC50.
The ensemble averages target fractional ranks within each patient's entire candidate population,
averaging tied scores; 0 best, 1 worst. Top-20 ties use genomic ID. Labels join after ranking.
Vaccine recall denominator = all 98 published patient/IMP table components (13,17,14,14,20,20),
including unreconstructable components. Paper reports 97 unique neoantigens; discrepancy persists.
Peptide recognition among selected vaccine overlap is the week-16 positive component count divided
by selected vaccine components with usable labels. A positive is any CD4 or CD8 peptide IFNg ELISpot,
ex vivo or expanded, using author calls; a negative requires explicit negatives in both cell types.
This does not measure minigene recognition, tumor recognition, mutant selectivity, or cell killing.
ROC AUC includes all reconstructable tested vaccinated mutations, with no top-K restriction.
ROC ties get half credit, without genomic-ID tie breaking. Macro AUC averages patient AUCs;
pair-weighted AUC pools only within-patient positive-negative concordance pairs.
These class-I models do not directly predict CD4 responses. The benchmark is retrospective and
conditional on historical vaccine selection and reconstruction. Training overlap must be audited;
results are descriptive and do not establish independent predictive performance or clinical utility.
"""


def read_csv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle))


def bigmhc_targets(epitopes, alleles, binding, predictions):
    """Reuse precisely the established best-pair and within-patient rank implementation."""
    rows = rank_targets(epitopes, alleles, binding, predictions, {a: a for a in alleles})
    return [
        {
            key.replace("deepimmuno", "bigmhc").replace("mean_ensemble", "bigmhc_ensemble"): value
            for key, value in row.items()
        }
        for row in rows
    ]


def verify_baseline(targets, saved):
    by_mutation = {row["mutation_id"]: row for row in saved}
    if {row["mutation_id"] for row in targets} != set(by_mutation):
        raise ValueError("Regenerated baseline mutation population changed")
    for row in targets:
        previous = by_mutation[row["mutation_id"]]
        for key, value in row.items():
            expected = previous[key]
            if isinstance(value, bool):
                expected = expected == "True"
            elif isinstance(value, (int, float)):
                expected = float(expected)
            if value != expected:
                raise ValueError(f"Regenerated baseline mismatch: {row['mutation_id']} / {key}")


def score_summary(values):
    values = sorted(values)
    quantiles = statistics.quantiles(values, n=100, method="inclusive")
    return {
        "n": len(values),
        "min": values[0],
        "max": values[-1],
        "mean": statistics.mean(values),
        "p01": quantiles[0],
        "p25": quantiles[24],
        "median": statistics.median(values),
        "p75": quantiles[74],
        "p99": quantiles[98],
        "unique": len(set(values)),
    }


def evaluate_patient(patient, targets, components, observations):
    labeled, excluded = assay_labels(components, targets, observations)
    labels = {row["mutation_id"]: row["label"] for row in labeled}
    recall, auc = [], []
    available = {row["mutation_id"] for row in targets}
    for method in PILOT_METHODS:
        selected = {row["mutation_id"] for row in targets if row[f"{method}_rank"] <= 20}
        overlap = [row for row in components if row["mutation_id"] in selected]
        known = [row for row in overlap if row["mutation_id"] in labels]
        recognized = sum(labels[row["mutation_id"]] for row in known)
        recall.append(
            {
                "patient": patient,
                "method": method,
                "k": 20,
                "candidate_targets": len(targets),
                "selected_targets": len(selected),
                "vaccine_hits": len(overlap),
                "vaccine_denominator": len(components),
                "vaccine_recall": len(overlap) / len(components),
                "reconstructable_vaccine_components": sum(r["mutation_id"] in available for r in components),
                "immunogenic_hits": sum(row["immunogenic"] for row in overlap),
                "immunogenic_denominator": sum(row["immunogenic"] for row in components),
                "recognized_selected_overlap": recognized,
                "tested_selected_overlap": len(known),
                "unknown_selected_overlap": len(overlap) - len(known),
                "recognition_fraction": recognized / len(known) if known else None,
                "hit_imp_ids": [row["imp_id"] for row in overlap],
            }
        )
        auc.append(
            {
                "patient": patient,
                "method": method,
                **roc_metrics([r["label"] for r in labeled], [r[f"{method}_cost"] for r in labeled]),
            }
        )
    shifts = []
    by_target = {row["mutation_id"]: row for row in targets}
    for old, new in [("deepimmuno", "bigmhc"), ("mean_ensemble", "bigmhc_ensemble")]:
        old_ids = {r["mutation_id"] for r in targets if r[f"{old}_rank"] <= 20}
        new_ids = {r["mutation_id"] for r in targets if r[f"{new}_rank"] <= 20}
        for mutation in sorted(old_ids | new_ids):
            row = by_target[mutation]
            shifts.append(
                {
                    "patient": patient,
                    "comparison": f"{old} -> {new}",
                    "mutation_id": mutation,
                    "gene": row["gene"],
                    "protein_change": row["protein_change"],
                    "membership": "retained"
                    if mutation in old_ids & new_ids
                    else "entered"
                    if mutation in new_ids
                    else "left",
                    "old_rank": row[f"{old}_rank"],
                    "new_rank": row[f"{new}_rank"],
                    "vaccine_imp_ids": [c["imp_id"] for c in components if c["mutation_id"] == mutation],
                }
            )
    return recall, auc, labeled, excluded, shifts


def render(output, aggregate, recall, auc, coverage, overlap):
    def table(rows, keys):
        head = "".join(f"<th>{html.escape(key.replace('_', ' '))}</th>" for key in keys)
        body = "".join(
            "<tr>" + "".join(f"<td>{html.escape(str(row[key]))}</td>" for key in keys) + "</tr>"
            for row in rows
        )
        return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"

    display = [
        {
            **row,
            "method": LABELS[row["method"]],
            "macro_auc": f"{row['macro_auc']:.4f}",
            "within_patient_pair_weighted_auc": f"{row['within_patient_pair_weighted_auc']:.4f}",
        }
        for row in aggregate
    ]
    auc_display = [
        {
            "patient": patient,
            **{
                LABELS[
                    method
                ]: f"{next(r['roc_auc'] for r in auc if r['patient'] == patient and r['method'] == method):.4f}"
                for method in PILOT_METHODS
            },
        }
        for patient in range(1, 7)
    ]
    overlap_summary = overlap.get("evaluation_summary", {})
    contamination = (
        (
            f"Training overlap: {overlap_summary['with_training_pair_overlap']}/{overlap_summary['targets']} evaluated mutations "
            f"have an exact candidate peptide–HLA pair in BigMHC IM training data; "
            f"{overlap_summary['with_training_or_validation_pair_overlap']}/{overlap_summary['targets']} overlap training or validation. "
            f"Removing these leaves {overlap_summary['remaining_positive']} positive and {overlap_summary['remaining_negative']} negative, "
            "so an uncontaminated AUC cannot be computed. This is a contaminated retrospective scoring comparison, "
            "not evidence that BigMHC improves independent prediction."
        )
        if overlap_summary
        else "Training overlap audit is unavailable; independent predictive performance is unestablished."
    )
    page = f"""<!doctype html><html lang="en"><meta charset="utf-8"><title>BigMHC IM six-patient pilot</title>
<style>body{{font:16px/1.5 system-ui;max-width:1250px;margin:40px auto;padding:0 24px;color:#223042}}
table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{text-align:left;padding:10px;border-bottom:1px solid #ddd}}
pre{{white-space:pre-wrap;background:#f2f5f8;padding:20px}}</style><h1>BigMHC IM six-patient pilot</h1>
<p>Optional retrospective comparison. Existing DeepImmuno pipeline unchanged. Research scores, not clinical recommendations.</p>
<p style="padding:18px;background:#fff1d6;border-left:5px solid #ad661c"><strong>{html.escape(contamination)}</strong>
<a href="audit/overlap.json">Exact overlap audit</a> · <a href="https://data.mendeley.com/datasets/dvmz6pkzvb/4">Released BigMHC datasets</a>.</p>
{table(display, ["method", "vaccine_hits", "vaccine_denominator", "recognized_selected_overlap", "tested_selected_overlap", "macro_auc", "within_patient_pair_weighted_auc"])}
<p>Recognition numerator/denominator refers only to selected vaccine overlaps with usable week-16 peptide labels.
It is not the fraction of all 120 proposed targets that are immunogenic: most unvaccinated targets were never tested.</p>
<h2>ROC AUC on all 72 reconstructable tested mutations</h2>
{table(auc_display, ["patient", *LABELS.values()])}
<h2>Per-patient top-20 recovery and recognition</h2>
{table(recall, ["patient", "method", "vaccine_hits", "vaccine_denominator", "recognized_selected_overlap", "tested_selected_overlap"])}
<h2>Coverage</h2>{table(coverage, ["patient", "targets", "unique_peptides", "pairs", "exact_hla_alleles"])}
<h2>Methodology</h2><pre>{html.escape(METHODOLOGY)}</pre>
<p><a href="analysis.json">Full results/provenance</a> · <a href="target-rankings.csv">Target rankings</a> ·
<a href="top20-shifts.csv">Top-20 membership shifts</a> · <a href="labeled-targets.csv">Assay labels</a> ·
<a href="score-summaries.csv">Raw score summaries</a> · <a href="inference/inference.json">Inference provenance</a></p></html>"""
    (output / "report.html").write_text(page)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", default="data/outputs/cohort-recall")
    parser.add_argument("--output-dir", default="data/outputs/bigmhc-pilot")
    parser.add_argument("--model-root", default="data/outputs/bigmhc-pilot/upstream")
    parser.add_argument("--python", default="data/outputs/bigmhc-pilot/runtime/bin/python")
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--device", choices=["cpu", "mps"], default="cpu")
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    source, output = Path(args.baseline).resolve(), Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    hashes = json.loads((source / "artifact-hashes.json").read_text())
    input_names = ["analysis.json", "target-rankings.csv", "vaccine-component-audit.csv"]
    for name in input_names:
        if digest(source / name) != hashes[name]:
            raise ValueError(f"Baseline input changed: {name}")
    baseline = json.loads((source / "analysis.json").read_text())
    saved_targets = read_csv(source / "target-rankings.csv")
    reference = load_reference(baseline["patients"][0]["config"]["data_dir"])
    patients, all_pairs = [], set()
    for record in baseline["patients"]:
        patient = record["patient"]
        print(f"P{patient}: preparing and verifying unchanged baseline", flush=True)
        config = Config(**{**record["config"], "cached_only": True})
        variants, alleles, published, assays = patient_data(config.data_dir, patient)
        if alleles != record["alleles"]:
            raise ValueError(f"Patient HLA changed: P{patient}")
        proteins, _ = reconstruct(variants, reference)
        epitopes, _, filter_summary = filter_reference_peptides(
            enumerate_epitopes(proteins, config.lengths), reference
        )
        for key, value in filter_summary.items():
            if record["reference_filter"][key] != value:
                raise ValueError(f"Reference-filter mismatch: P{patient} / {key}")
        cache = BindingCache(config.cache_dir)
        try:
            binding = cache.lookup(alleles)
        finally:
            cache.close()
        immuno, provenance = predict_deepimmuno(epitopes, alleles, config, allow_hla_rescue=True)
        targets = rank_targets(epitopes, alleles, binding, immuno, provenance["allele_representations"])
        verify_baseline(targets, [row for row in saved_targets if int(row["patient"]) == patient])
        components = audit_components(actual_vaccine(published, variants, proteins), proteins, assays)
        pairs = {(row["peptide"], allele) for row in epitopes for allele in alleles}
        all_pairs.update(pairs)
        patients.append(
            {
                "patient": patient,
                "alleles": alleles,
                "epitopes": epitopes,
                "binding": binding,
                "deepimmuno": immuno,
                "targets": targets,
                "components": components,
                "assays": assays,
                "pairs": pairs,
            }
        )
    with (output / "candidate-pairs.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["peptide", "allele"])
        writer.writerows(sorted(all_pairs))
    print(f"Verified all six baseline populations: {len(all_pairs)} distinct peptide-HLA pairs", flush=True)
    exact_alleles(args.model_root, [allele for _, allele in all_pairs])
    if args.prepare_only:
        return
    predictions, inference = predict_immunogenicity(
        all_pairs,
        args.model_root,
        args.python,
        output / "inference",
        args.batch_size,
        args.threads,
        args.device,
    )
    all_targets, all_components, all_labeled, all_excluded, recall, auc, shifts, coverage, summaries = (
        [] for _ in range(9)
    )
    for data in patients:
        patient = data["patient"]
        big = {
            r["mutation_id"]: r
            for r in bigmhc_targets(data["epitopes"], data["alleles"], data["binding"], predictions)
        }
        targets = [{"patient": patient, **row, **big[row["mutation_id"]]} for row in data["targets"]]
        results = evaluate_patient(patient, targets, data["components"], data["assays"])
        for destination, rows in zip([recall, auc, all_labeled, all_excluded, shifts], results):
            destination.extend(rows)
        all_targets.extend(targets)
        all_components.extend(data["components"])
        coverage.append(
            {
                "patient": patient,
                "targets": len(targets),
                "unique_peptides": len({p for p, _ in data["pairs"]}),
                "pairs": len(data["pairs"]),
                "exact_hla_alleles": data["alleles"],
                "scored_pairs": len(data["pairs"]),
                "rescued_hla_alleles": [],
                "unsupported_pairs": 0,
            }
        )
        for model, pair_values, target_key in [
            ("bigmhc", [predictions[pair] for pair in data["pairs"]], "bigmhc_score"),
            ("deepimmuno", [data["deepimmuno"][pair] for pair in data["pairs"]], "deepimmuno_score"),
            ("netmhcpan", [data["binding"][pair]["ic50"] for pair in data["pairs"]], "binding_ic50_nm"),
        ]:
            for scope, values in [
                ("peptide-HLA", pair_values),
                ("mutation", [r[target_key] for r in targets]),
            ]:
                summaries.append(
                    {"patient": patient, "model": model, "scope": scope, **score_summary(values)}
                )
    aggregate = []
    for method in PILOT_METHODS:
        ar = [r for r in auc if r["method"] == method and r["roc_auc"] is not None]
        rr = [r for r in recall if r["method"] == method]
        pairs = sum(r["positive_negative_pairs"] for r in ar)
        credits = sum(r["correctly_ordered_pairs"] + 0.5 * r["tied_pairs"] for r in ar)
        sums = {
            key: sum(r[key] for r in rr)
            for key in [
                "vaccine_hits",
                "vaccine_denominator",
                "immunogenic_hits",
                "immunogenic_denominator",
                "recognized_selected_overlap",
                "tested_selected_overlap",
            ]
        }
        aggregate.append(
            {
                "method": method,
                **sums,
                "vaccine_recall": sums["vaccine_hits"] / sums["vaccine_denominator"],
                "macro_vaccine_recall": statistics.mean(r["vaccine_recall"] for r in rr),
                "recognition_fraction": sums["recognized_selected_overlap"] / sums["tested_selected_overlap"]
                if sums["tested_selected_overlap"]
                else None,
                "macro_auc": statistics.mean(r["roc_auc"] for r in ar),
                "within_patient_pair_weighted_auc": credits / pairs,
                "within_patient_pairs": pairs,
                "correct_pair_credits": credits,
            }
        )
    result = {
        "created_at": datetime.datetime.now(datetime.UTC).isoformat(),
        "methodology": METHODOLOGY,
        "baseline_inputs": {str(source / name): hashes[name] for name in input_names},
        "baseline_verified": True,
        "sources": verify_sources(Path(config.data_dir)),
        "inference": inference,
        "coverage": coverage,
        "unique_inference_pairs": len(all_pairs),
        "labeled_targets": len(all_labeled),
        "positives": sum(r["label"] for r in all_labeled),
        "negatives": sum(r["label"] == 0 for r in all_labeled),
        "recall": recall,
        "auc": auc,
        "aggregate": aggregate,
        "score_summaries": summaries,
    }
    overlap_path = output / "audit/overlap.json"
    overlap = json.loads(overlap_path.read_text()) if overlap_path.exists() else {}
    if overlap:
        if overlap["baseline_hashes"] != {name: hashes[name] for name in input_names}:
            raise ValueError("Training-overlap audit used a different baseline")
        result["training_overlap"] = {
            "path": str(overlap_path),
            "sha256": digest(overlap_path),
            "evaluation_summary": overlap["evaluation_summary"],
        }
        if overlap["evaluation_summary"]["with_training_or_validation_pair_overlap"]:
            result["evaluation_status"] = (
                "Disqualified as an independent benchmark due to confirmed training/validation overlap; "
                "retained as descriptive pilot only."
            )
    for name, rows in [
        ("aggregate", aggregate),
        ("recall", recall),
        ("auc", auc),
        ("target-rankings", all_targets),
        ("vaccine-components", all_components),
        ("labeled-targets", all_labeled),
        ("excluded", all_excluded),
        ("top20-shifts", shifts),
        ("coverage", coverage),
        ("score-summaries", summaries),
    ]:
        write_csv(output / f"{name}.csv", rows)
    write_json(output / "analysis.json", result)
    render(output, aggregate, recall, auc, coverage, overlap)
    write_json(
        output / "artifact-hashes.json",
        {p.name: digest(p) for p in output.iterdir() if p.is_file() and p.name != "artifact-hashes.json"},
    )
    print(json.dumps(aggregate, indent=2), flush=True)


if __name__ == "__main__":
    main()
