"""Audit exact sequence overlap with the released BigMHC immunogenicity datasets."""

import argparse
import csv
import io
import json
import zipfile
from collections import defaultdict
from pathlib import Path

from neovax.assay_auc import assay_labels
from neovax.config import Config
from neovax.data import load_reference, patient_data
from neovax.io import digest, write_json
from neovax.peptides import enumerate_epitopes, reconstruct
from neovax.reference_filter import filter_reference_peptides

DATASET_SHA256 = "0d152a452756cf2e0014ffced6afc25118e7c11cf1f626b26e49f50f79edffaa"
DATASET_URL = "https://data.mendeley.com/datasets/dvmz6pkzvb/4"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", default="data/outputs/bigmhc-pilot/audit/datasets.zip")
    parser.add_argument("--baseline-dir", default="data/outputs/cohort-recall")
    parser.add_argument("--output", default="data/outputs/bigmhc-pilot/audit/overlap.json")
    args = parser.parse_args()
    archive = Path(args.archive)
    if digest(archive) != DATASET_SHA256:
        raise ValueError("BigMHC dataset archive does not match its published SHA256")
    baseline = Path(args.baseline_dir)
    hashes = json.loads((baseline / "artifact-hashes.json").read_text())
    for name in ["analysis.json", "vaccine-component-audit.csv", "target-rankings.csv"]:
        if digest(baseline / name) != hashes[name]:
            raise ValueError(f"Baseline artifact changed: {name}")
    metadata = json.loads((baseline / "analysis.json").read_text())
    with (baseline / "vaccine-component-audit.csv").open() as handle:
        components = list(csv.DictReader(handle))
    with (baseline / "target-rankings.csv").open() as handle:
        targets = list(csv.DictReader(handle))
    datasets = {}
    with zipfile.ZipFile(archive) as source:
        for split in ["train", "val", "test"]:
            filename = next(n for n in source.namelist() if n.endswith(f"im_{split}.csv"))
            datasets[split] = list(csv.DictReader(io.StringIO(source.read(filename).decode())))
    results, evaluation_overlap = [], []
    for patient in range(1, 7):
        saved = next(p for p in metadata["patients"] if p["patient"] == patient)
        config = Config(**saved["config"])
        variants, alleles, _, observations = patient_data(config.data_dir, patient)
        reference = load_reference(config.data_dir)
        proteins, _ = reconstruct(variants, reference)
        epitopes, _, _ = filter_reference_peptides(enumerate_epitopes(proteins, config.lengths), reference)
        by_sequence = defaultdict(set)
        for epitope in epitopes:
            by_sequence[epitope["peptide"]].add(epitope["mutation_id"])
        expected = {r["mutation_id"] for r in targets if int(r["patient"]) == patient}
        assert set().union(*by_sequence.values()) == expected
        vaccine = [r for r in components if int(r["patient"]) == patient]
        labeled, _ = assay_labels(vaccine, [r for r in targets if int(r["patient"]) == patient], observations)
        for split, rows in datasets.items():
            matches, payload_matches = [], []
            for row in rows:
                peptide, allele = row["pep"], row["mhc"]
                mutations = by_sequence.get(peptide, set())
                if mutations:
                    matches.append({
                        "peptide": peptide, "hla": allele, "dataset_label": row["tgt"],
                        "exact_patient_hla": allele in alleles,
                        "mutation_ids": sorted(mutations),
                        "vaccine_imp_ids": sorted(r["imp_id"] for r in vaccine if r["mutation_id"] in mutations),
                    })
                payloads = [r["imp_id"] for r in vaccine if peptide in r["immunizing_sequence"]]
                if payloads:
                    payload_matches.append({
                        "peptide": peptide, "hla": allele, "dataset_label": row["tgt"],
                        "exact_patient_hla": allele in alleles, "vaccine_imp_ids": sorted(payloads),
                    })
            results.append({
                "patient": patient, "split": split, "dataset_rows": len(rows),
                "candidate_sequence_matches": len({r["peptide"] for r in matches}),
                "candidate_pair_matches": len({(r["peptide"], r["hla"]) for r in matches if r["exact_patient_hla"]}),
                "candidate_mutations_with_sequence_overlap": sorted({m for r in matches for m in r["mutation_ids"]}),
                "candidate_mutations_with_pair_overlap": sorted({m for r in matches if r["exact_patient_hla"] for m in r["mutation_ids"]}),
                "matches": matches, "vaccine_payload_matches": payload_matches,
            })
        patient_splits = {r["split"]: r for r in results if r["patient"] == patient}
        exposed = {
            split: set(patient_splits[split]["candidate_mutations_with_pair_overlap"])
            for split in ["train", "val"]
        }
        exposed_union = exposed["train"] | exposed["val"]
        for row in labeled:
            evaluation_overlap.append({
                "patient": patient, "mutation_id": row["mutation_id"], "gene": row["gene"],
                "protein_change": row["protein_change"], "label": row["label"],
                "imp_ids": row["imp_ids"],
                "train_pair_overlap": row["mutation_id"] in exposed["train"],
                "validation_pair_overlap": row["mutation_id"] in exposed["val"],
                "train_or_validation_pair_overlap": row["mutation_id"] in exposed_union,
            })
    result = {
        "dataset_url": DATASET_URL, "archive_sha256": DATASET_SHA256,
        "baseline_hashes": {name: hashes[name] for name in ["analysis.json", "vaccine-component-audit.csv", "target-rankings.csv"]},
        "scope": "Exact sequences and exact peptide-HLA pairs against released IM train/validation/test splits. Training and validation overlap are potential benchmark contamination; test overlap alone is not evidence of training exposure. Payload matches include unchanged flanks and lengths outside our 9/10-mer inputs, so they do not necessarily identify the same mutant epitope. Related sequences, other training/development data, and EL pretraining overlap are not audited.",
        "patients": results, "evaluation_targets": evaluation_overlap,
        "evaluation_summary": {
            "targets": len(evaluation_overlap),
            "with_training_pair_overlap": sum(r["train_pair_overlap"] for r in evaluation_overlap),
            "with_training_or_validation_pair_overlap": sum(r["train_or_validation_pair_overlap"] for r in evaluation_overlap),
            "remaining_positive": sum(r["label"] == 1 and not r["train_or_validation_pair_overlap"] for r in evaluation_overlap),
            "remaining_negative": sum(r["label"] == 0 and not r["train_or_validation_pair_overlap"] for r in evaluation_overlap),
        },
    }
    write_json(args.output, result)
    for row in results:
        print(f"P{row['patient']} {row['split']}: {row['candidate_sequence_matches']} sequences / {row['candidate_pair_matches']} peptide-HLA pairs")


if __name__ == "__main__":
    main()
