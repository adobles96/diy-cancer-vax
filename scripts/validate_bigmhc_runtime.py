"""Validate MPS pilot scores and mutation ranking against official CPU inference."""

import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path

from neovax.assay_auc import roc_metrics
from neovax.bigmhc import predict_immunogenicity, read_predictions
from neovax.io import digest, write_json
from neovax.ranking import fractional_ranks


def read_csv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="data/outputs/bigmhc-pilot")
    args = parser.parse_args()
    output = Path(args.output_dir).resolve()
    analysis = json.loads((output / "analysis.json").read_text())
    inference = analysis["inference"]
    targets = read_csv(output / "target-rankings.csv")
    labeled = read_csv(output / "labeled-targets.csv")
    rows = read_csv(output / "inference/input.csv")
    pairs = [(r["pep"], r["mhc"]) for r in rows]
    gpu = read_predictions(output / "inference/predictions.csv", pairs)
    groups = defaultdict(list)
    for pair in pairs:
        groups[(pair[1], len(pair[0]))].append(pair)
    rng = random.Random(20261028)
    validation_pairs = {(r["bigmhc_peptide"], r["bigmhc_actual_allele"]) for r in targets}
    for group in groups.values():
        validation_pairs.update(rng.sample(group, min(32, len(group))))
    cpu, cpu_metadata = predict_immunogenicity(
        validation_pairs,
        inference["model"]["root"],
        inference["command"][0],
        output / "cpu-validation",
        batch_size=512,
        threads=8,
        device="cpu",
    )
    differences = [abs(cpu[pair] - gpu[pair]) for pair in sorted(validation_pairs)]
    comparisons = []
    for patient in range(1, 7):
        patient_targets = [r for r in targets if int(r["patient"]) == patient]
        scores = [cpu[(r["bigmhc_peptide"], r["bigmhc_actual_allele"])] for r in patient_targets]
        costs = fractional_ranks(scores, descending=True)
        corrected = {}
        for row, cost in zip(patient_targets, costs):
            corrected[row["mutation_id"]] = {
                "bigmhc": cost,
                "bigmhc_ensemble": (cost + float(row["netmhcpan_cost"])) / 2,
            }
        labels = [r for r in labeled if int(r["patient"]) == patient]
        for method in ["bigmhc", "bigmhc_ensemble"]:
            original_top = {r["mutation_id"] for r in patient_targets if int(r[f"{method}_rank"]) <= 20}
            corrected_top = set(
                sorted(corrected, key=lambda mutation: (corrected[mutation][method], mutation))[:20]
            )
            old_auc = roc_metrics(
                [int(r["label"]) for r in labels], [float(r[f"{method}_cost"]) for r in labels]
            )["roc_auc"]
            new_auc = roc_metrics(
                [int(r["label"]) for r in labels], [corrected[r["mutation_id"]][method] for r in labels]
            )["roc_auc"]
            comparisons.append(
                {
                    "patient": patient,
                    "method": method,
                    "top20_identical": original_top == corrected_top,
                    "original_auc": old_auc,
                    "cpu_winner_auc": new_auc,
                    "changed_target_costs": sum(
                        corrected[r["mutation_id"]][method] != float(r[f"{method}_cost"])
                        for r in patient_targets
                    ),
                }
            )
    result = {
        "validation_pairs": len(validation_pairs),
        "mutation_winner_pairs": len({(r["bigmhc_peptide"], r["bigmhc_actual_allele"]) for r in targets}),
        "allele_length_strata": len(groups),
        "max_absolute_difference": max(differences),
        "mean_absolute_difference": sum(differences) / len(differences),
        "absolute_tolerance": 0.00001,
        "all_within_tolerance": max(differences) <= 0.00001,
        "ranking_comparisons": comparisons,
        "cpu_inference": cpu_metadata,
        "mps_predictions_sha256": digest(output / "inference/predictions.csv"),
        "scope": "All MPS mutation-winning pairs plus 32 random pairs per actual allele/length stratum were rescored with official CPU predict.py. Mutation maxima were replaced with their CPU-rescored MPS winner scores for rank and AUC comparison; alternative nonwinning pairs were not all CPU-rescored, so this is a numerical agreement check, not full CPU replay.",
    }
    write_json(output / "runtime-validation.json", result)
    print(json.dumps({k: v for k, v in result.items() if k != "cpu_inference"}, indent=2))
    if not result["all_within_tolerance"] or any(
        not r["top20_identical"] or r["original_auc"] != r["cpu_winner_auc"] for r in comparisons
    ):
        raise ValueError("CPU/MPS score or decision agreement requires investigation")


if __name__ == "__main__":
    main()
