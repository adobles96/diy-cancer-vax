"""Audit exact candidate overlap against the official NetMHCpan-4.1 BA/EL release."""

import argparse
import csv
import hashlib
import json
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

from neovax.config import Config
from neovax.data import load_reference, patient_data
from neovax.io import digest, write_csv, write_json
from neovax.peptides import enumerate_epitopes, reconstruct
from neovax.reference_filter import filter_reference_peptides

SOURCE_URL = "https://services.healthtech.dtu.dk/suppl/immunology/NAR_NetMHCpan_NetMHCIIpan/"
ARCHIVE_URL = f"{SOURCE_URL}NetMHCpan_train.tar.gz"
ARCHIVE_SHA256 = "06f2c9f20bb959238bf5d601fca0489a0ed3f17648b952f30640205afca8f9b4"


def normalize_hla(allele):
    """Normalize the literal asterisk difference; preserve allele resolution."""
    return allele.replace("*", "")


def read_csv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def candidate_context(baseline, pair_file, labels_file):
    hashes = json.loads((baseline / "artifact-hashes.json").read_text())
    for name in ["analysis.json", "vaccine-component-audit.csv", "target-rankings.csv"]:
        if digest(baseline / name) != hashes[name]:
            raise ValueError(f"Baseline artifact changed: {name}")
    metadata = json.loads((baseline / "analysis.json").read_text())
    targets = read_csv(baseline / "target-rankings.csv")
    components = read_csv(baseline / "vaccine-component-audit.csv")
    labels = read_csv(labels_file)
    actual_pairs = {(r["peptide"], normalize_hla(r["allele"])) for r in read_csv(pair_file)}
    reference = load_reference(metadata["patients"][0]["config"]["data_dir"])
    patients, reconstructed_pairs = {}, set()
    by_sequence = defaultdict(list)
    payloads = defaultdict(list)
    for saved in metadata["patients"]:
        patient = saved["patient"]
        config = Config(**saved["config"])
        variants, alleles, _, _ = patient_data(config.data_dir, patient)
        proteins, _ = reconstruct(variants, reference)
        epitopes, _, _ = filter_reference_peptides(enumerate_epitopes(proteins, config.lengths), reference)
        mutations = defaultdict(set)
        for epitope in epitopes:
            mutations[epitope["peptide"]].add(epitope["mutation_id"])
        expected = {r["mutation_id"] for r in targets if int(r["patient"]) == patient}
        if set().union(*mutations.values()) != expected:
            raise ValueError(f"Reconstructed mutations disagree with baseline for P{patient}")
        normalized = {normalize_hla(a) for a in alleles}
        pairs = {(p, a) for p in mutations for a in normalized}
        reconstructed_pairs.update(pairs)
        patients[patient] = {
            "patient": patient, "candidate_sequences": len(mutations), "candidate_pairs": len(pairs),
            "candidate_mutations": len(expected), "alleles": sorted(alleles),
        }
        for peptide, ids in mutations.items():
            by_sequence[peptide].append((patient, normalized, ids))
        for row in components:
            if int(row["patient"]) != patient:
                continue
            sequence = row["immunizing_sequence"]
            for length in range(1, len(sequence) + 1):
                for start in range(len(sequence) - length + 1):
                    payloads[sequence[start:start + length]].append((patient, row["imp_id"], normalized))
    if reconstructed_pairs != actual_pairs:
        raise ValueError("Reconstructed actual peptide-HLA pairs differ from saved candidate inputs")
    target_ids = {(int(r["patient"]), r["mutation_id"]) for r in targets}
    if any((int(r["patient"]), r["mutation_id"]) not in target_ids for r in labels):
        raise ValueError("An evaluation label does not map to a reconstructed candidate target")
    inputs = {str(pair_file): digest(pair_file), str(labels_file): digest(labels_file)}
    inputs.update({str(baseline / n): hashes[n] for n in [
        "analysis.json", "vaccine-component-audit.csv", "target-rankings.csv",
    ]})
    return patients, by_sequence, payloads, labels, actual_pairs, inputs, len(components), reference


def scan_archive(archive, by_sequence, payloads):
    matches, payload_matches, inventory = [], [], []
    with tarfile.open(archive, "r:gz") as source:
        members = {Path(m.name).name: m for m in source.getmembers() if m.isfile()}
        expected = {f"c00{i}_{kind}" for i in range(5) for kind in ["ba", "el"]}
        if set(members) != expected | {"allelelist", "MHC_pseudo.dat"}:
            raise ValueError(f"Unexpected archive members: {sorted(members)}")
        allele_bytes = source.extractfile(members["allelelist"]).read()
        allele_map = {}
        for line in allele_bytes.decode().splitlines():
            cell, alleles = line.split()
            allele_map[cell] = {normalize_hla(a) for a in alleles.split(",")}
        pseudo_bytes = source.extractfile(members["MHC_pseudo.dat"]).read()
        known_alleles = {normalize_hla(line.split()[0]) for line in pseudo_bytes.decode().splitlines()}
        if not set().union(*allele_map.values()) <= known_alleles:
            raise ValueError("An allelelist membership has no released MHC pseudosequence")
        for name in sorted(members):
            member = members[name]
            file_hash = hashlib.sha256()
            counts, label_counts = Counter(), Counter()
            with source.extractfile(member) as handle:
                for line_number, line in enumerate(handle, 1):
                    file_hash.update(line)
                    if name not in expected:
                        continue
                    peptide, label, identifier = line.decode().split()
                    counts["rows"] += 1
                    is_ma = normalize_hla(identifier) not in known_alleles and identifier in allele_map
                    counts["ma_rows" if is_ma else "sa_rows"] += 1
                    if name.endswith("_el"):
                        label_counts[label] += 1
                    hla = allele_map[identifier] if is_ma else {normalize_hla(identifier)}
                    if not is_ma and normalize_hla(identifier) not in known_alleles:
                        raise ValueError(f"Unknown allele or unmapped cell line: {identifier}")
                    base = {
                        "file": name, "line": line_number, "data_type": name[-2:].upper(),
                        "peptide": peptide, "dataset_target": label, "training_identifier": identifier,
                        "annotation": "MA cell-line membership" if is_ma else "single-allele annotation",
                        "training_alleles": sorted(hla),
                    }
                    for patient, actual_hla, mutations in by_sequence.get(peptide, []):
                        matches.append({
                            **base, "patient": patient, "mutation_ids": sorted(mutations),
                            "matching_patient_alleles": sorted(hla & actual_hla),
                        })
                    seen_payloads = set()
                    for patient, imp_id, actual_hla in payloads.get(peptide, []):
                        if (patient, imp_id) in seen_payloads:
                            continue
                        seen_payloads.add((patient, imp_id))
                        payload_matches.append({
                            **base, "patient": patient, "imp_id": imp_id,
                            "matching_patient_alleles": sorted(hla & actual_hla),
                        })
            inventory.append({
                "file": name, "size": member.size, "sha256": file_hash.hexdigest(),
                **dict(counts), "dataset_target_counts": dict(label_counts),
            })
            print(f"Scanned {name}: {counts['rows']} training rows", flush=True)
    return matches, payload_matches, inventory, allele_map


def summarize(matches):
    sequences = {r["peptide"] for r in matches}
    pairs = {(r["peptide"], a) for r in matches for a in r["matching_patient_alleles"]}
    mutations = {(r["patient"], m) for r in matches for m in r["mutation_ids"]}
    pair_mutations = {(r["patient"], m) for r in matches if r["matching_patient_alleles"] for m in r["mutation_ids"]}
    return {
        "matching_training_rows_by_patient": len(matches), "candidate_sequences": len(sequences),
        "candidate_pairs_with_sa_or_ma_membership": len(pairs), "mutations_with_sequence_overlap": len(mutations),
        "mutations_with_sa_or_ma_pair_overlap": len(pair_mutations),
        "sa_pairs": len({(r["peptide"], a) for r in matches if r["annotation"] == "single-allele annotation"
                         for a in r["matching_patient_alleles"]}),
        "ma_membership_pairs": len({(r["peptide"], a) for r in matches if r["annotation"] == "MA cell-line membership"
                                    for a in r["matching_patient_alleles"]}),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=Path("data/outputs/training-overlap/netmhcpan/NetMHCpan_train.tar.gz"))
    parser.add_argument("--baseline-dir", type=Path, default=Path("data/outputs/cohort-recall"))
    parser.add_argument("--candidate-pairs", type=Path, default=Path("data/outputs/bigmhc-pilot/candidate-pairs.csv"))
    parser.add_argument("--labels", type=Path, default=Path("data/outputs/cohort-auc/labeled-targets.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/outputs/training-overlap/netmhcpan"))
    args = parser.parse_args()
    archive_hash = digest(args.archive)
    if archive_hash != ARCHIVE_SHA256:
        raise ValueError(f"NetMHCpan-4.1 archive SHA256 mismatch: {archive_hash}")
    patients, by_sequence, payloads, labels, pairs, inputs, components, reference = candidate_context(
        args.baseline_dir, args.candidate_pairs, args.labels,
    )
    matches, payload_matches, inventory, allele_map = scan_archive(args.archive, by_sequence, payloads)
    payload_sequences = {r["peptide"] for r in payload_matches}
    _, reference_matches, _ = filter_reference_peptides(
        [{"peptide": p, "mutation_id": p} for p in sorted(payload_sequences)], reference,
    )
    reference_sequences = {r["peptide"] for r in reference_matches}
    for row in payload_matches:
        row["exact_sequence_in_pinned_reference"] = row["peptide"] in reference_sequences
    types = {kind: [r for r in matches if r["data_type"] == kind] for kind in ["BA", "EL"]}
    exposed = {
        (r["patient"], m, kind): True
        for kind, rows in types.items() for r in rows for m in r["mutation_ids"]
    }
    pair_exposed = {
        (r["patient"], m, kind): True for kind, rows in types.items()
        for r in rows if r["matching_patient_alleles"] for m in r["mutation_ids"]
    }
    evaluated = []
    for row in labels:
        patient, mutation = int(row["patient"]), row["mutation_id"]
        evaluated.append({
            "patient": patient, "mutation_id": mutation, "gene": row["gene"],
            "protein_change": row["protein_change"], "label": int(row["label"]), "imp_ids": row["imp_ids"],
            **{f"{kind.lower()}_sequence_overlap": (patient, mutation, kind) in exposed for kind in types},
            **{f"{kind.lower()}_sa_or_ma_pair_overlap": (patient, mutation, kind) in pair_exposed for kind in types},
        })
    per_patient = []
    for patient, context in patients.items():
        patient_matches = [r for r in matches if r["patient"] == patient]
        per_patient.append({**context, "combined": summarize(patient_matches), **{
            kind: summarize([r for r in patient_matches if r["data_type"] == kind]) for kind in types
        }})
    result = {
        "official_source": SOURCE_URL, "archive_url": ARCHIVE_URL, "archive_sha256": archive_hash,
        "archive_bytes": args.archive.stat().st_size, "input_sha256": inputs,
        "scope": (
            "Exact full peptide sequences against all five BA and all five EL partitions in the official "
            "2020 NetMHCpan-4.1 training release. Actual patient HLA pairs use only literal asterisk normalization. "
            "Single-allele annotations and MA cell-line HLA membership are distinguished: membership is not a "
            "confirmed experimental restriction or final deconvoluted assignment. BA and EL labels concern "
            "binding/presentation, not the vaccine immunogenicity outcome. Homology, older or additional "
            "development data, and the precise weights deployed behind the IEDB version name are not audited. "
            "Vaccine payload substring matches include unchanged flanks and lengths beyond the actual 9/10mer inputs."
        ),
        "candidate_scope": {
            "patients": len(patients), "unique_sequences": len(by_sequence), "unique_actual_hla_pairs": len(pairs),
            "patient_mutations": sum(p["candidate_mutations"] for p in patients.values()),
            "vaccine_components": components,
        },
        "training_inventory": inventory,
        "ma_cell_line_count": len(allele_map),
        "combined": summarize(matches), **{kind: summarize(rows) for kind, rows in types.items()},
        "patients": per_patient,
        "evaluation_summary": {
            "targets": len(evaluated), "positive": sum(r["label"] for r in evaluated),
            "negative": sum(r["label"] == 0 for r in evaluated),
            **{f"{kind.lower()}_sequence_overlap": sum(r[f"{kind.lower()}_sequence_overlap"] for r in evaluated)
               for kind in types},
            **{f"{kind.lower()}_sa_or_ma_pair_overlap": sum(r[f"{kind.lower()}_sa_or_ma_pair_overlap"] for r in evaluated)
               for kind in types},
            "any_sequence_overlap": sum(any(r[f"{k.lower()}_sequence_overlap"] for k in types) for r in evaluated),
            "any_sa_or_ma_pair_overlap": sum(any(r[f"{k.lower()}_sa_or_ma_pair_overlap"] for k in types) for r in evaluated),
        },
        "evaluation_targets": evaluated, "candidate_matches": matches,
        "vaccine_payload_summary": {
            "matching_sequences": len(payload_sequences),
            "matching_sequences_in_pinned_reference": len(reference_sequences),
            "components_with_substring_overlap": len({(r["patient"], r["imp_id"]) for r in payload_matches}),
            "components_with_sa_or_ma_hla_membership": len({
                (r["patient"], r["imp_id"]) for r in payload_matches if r["matching_patient_alleles"]
            }),
        },
        "vaccine_payload_substring_matches": payload_matches,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_json(args.output_dir / "overlap.json", result)
    write_csv(args.output_dir / "candidate-matches.csv", matches)
    write_csv(args.output_dir / "evaluation-targets.csv", evaluated)
    write_csv(args.output_dir / "vaccine-payload-matches.csv", payload_matches)
    print(json.dumps({k: result[k] for k in ["candidate_scope", "combined", "BA", "EL", "evaluation_summary"]}, indent=2))


if __name__ == "__main__":
    main()
