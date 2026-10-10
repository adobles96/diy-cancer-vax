from collections import Counter, defaultdict

from neovax.io import stable_id
from neovax.peptides import AMINO_ACIDS
from neovax.ranking import epitope_sort_key


def long_bounds(start, end, protein_length, requested_length):
    if requested_length < end - start:
        raise ValueError("Long peptide must contain the entire epitope")
    length = min(requested_length, protein_length)
    if length < 15:
        raise ValueError("Protein is shorter than the minimum 15 residues")
    left = max(0, min((start + end - length) // 2, protein_length - length))
    return left, left + length


def select_vaccine(ranked, mutation_scores, proteins, config):
    """Select from discovery records only, before joining any historical evidence."""
    protein_map = {r["variant_id"]: r for r in proteins}
    by_mutation, by_protein = defaultdict(list), defaultdict(list)
    for row in ranked:
        by_protein[row["variant_id"]].append(row)
        if row["eligible"]:
            by_mutation[row["mutation_id"]].append(row)
    if config.selection_unit == "mutation":
        candidates = [
            r
            for mutation in mutation_scores
            for r in sorted(by_mutation[mutation["mutation_id"]], key=epitope_sort_key)
        ]
    else:
        candidates = [r for r in ranked if r["eligible"]]
    costs = {r["mutation_id"]: r["ranking_cost"] for r in mutation_scores}
    counts, seen, selected, exclusions = Counter(), set(), [], []
    for epitope in candidates:
        mutation_id = epitope["mutation_id"]
        if counts[mutation_id] >= config.per_mutation_cap:
            continue
        protein = protein_map[epitope["variant_id"]]
        try:
            start, end = long_bounds(
                epitope["start"], epitope["end"], len(protein["mutant_protein"]), config.long_length
            )
        except ValueError as error:
            exclusions.append({"epitope_id": epitope["epitope_id"], "reason": str(error)})
            continue
        sequence = protein["mutant_protein"][start:end]
        if set(sequence) - AMINO_ACIDS:
            exclusions.append(
                {"epitope_id": epitope["epitope_id"], "reason": "ambiguous residues in long context"}
            )
            continue
        if sequence in seen:
            continue
        if not start <= epitope["start"] < epitope["end"] <= end or not start <= protein["position"] < end:
            raise ValueError("Long peptide containment invariant failed")
        covered = [
            r["epitope_id"]
            for r in by_protein[epitope["variant_id"]]
            if start <= r["start"] and r["end"] <= end and r["all_alleles_scored"]
        ]
        selected.append(
            {
                "component_id": stable_id(mutation_id, sequence),
                "selection_rank": len(selected) + 1,
                "mutation_id": mutation_id,
                "variant_id": epitope["variant_id"],
                "gene": epitope["gene"],
                "protein_change": epitope["protein_change"],
                "transcript": epitope["transcript"],
                "sequence": sequence,
                "reference_sequence": protein["reference_protein"][start:end],
                "start": start,
                "end": end,
                "mutation_offset": protein["position"] - start,
                "epitope_offset": epitope["start"] - start,
                "driving_epitope": epitope["peptide"],
                "epitope_id": epitope["epitope_id"],
                "best_allele": epitope["best_allele"],
                "ic50_nm": epitope["ic50_nm"],
                "deepimmuno_score": epitope["deepimmuno_score"],
                "ranking_cost": costs[mutation_id]
                if config.selection_unit == "mutation"
                else epitope["ranking_cost"],
                "score_origin": config.selection_unit,
                "covered_epitope_ids": covered,
                "length": len(sequence),
                "cysteine_count": sequence.count("C"),
                "hydrophobic_fraction": sum(a in "AILMFWVY" for a in sequence) / len(sequence),
            }
        )
        seen.add(sequence)
        counts[mutation_id] += 1
        if len(selected) >= config.budget:
            break
    return selected, {
        "requested": config.budget,
        "selected": len(selected),
        "shortfall": max(0, config.budget - len(selected)),
        "excluded": exclusions,
    }
