import statistics
from collections import defaultdict

from neovax.io import stable_id


def fractional_ranks(values, descending=False):
    ordered = sorted(enumerate(values), key=lambda item: item[1], reverse=descending)
    result = [0.0] * len(values)
    i = 0
    while i < len(ordered):
        j = i + 1
        while j < len(ordered) and ordered[j][1] == ordered[i][1]:
            j += 1
        score = ((i + j - 1) / 2) / (len(values) - 1) if len(values) > 1 else 0.0
        for index, _ in ordered[i:j]:
            result[index] = score
        i = j
    return result


def rank_epitopes(epitopes, alleles, predictions, config, immunogenicity=None):
    """No observation or historical membership may enter this function."""
    immunogenicity = immunogenicity or {}
    unique = sorted({(e["peptide"], allele) for e in epitopes for allele in alleles})
    pairs = []
    for peptide, allele in unique:
        prediction = predictions.get((peptide, allele))
        score = immunogenicity.get((peptide, allele))
        eligible = prediction is not None and (
            config.threshold is None or prediction["ic50"] < config.threshold
        )
        if config.ranking != "binding" and score is None:
            raise ValueError(
                f"Missing DeepImmuno result for {peptide} / {allele}; combined ranking unavailable"
            )
        pairs.append(
            {
                "pair_id": stable_id(peptide, allele),
                "peptide": peptide,
                "allele": allele,
                "ic50_nm": prediction["ic50"] if prediction else None,
                "percentile_rank": prediction["percentile_rank"] if prediction else None,
                "request_id": prediction["request_id"] if prediction else None,
                "deepimmuno_score": score,
                "eligible": eligible,
                "ranking_cost": None,
                "binding_fractional_rank": None,
                "immunogenicity_fractional_rank": None,
                "status": "scored" if prediction else "missing binding score",
            }
        )
    eligible_pairs = [p for p in pairs if p["eligible"]]
    binding_ranks = fractional_ranks([p["ic50_nm"] for p in eligible_pairs])
    immuno_ranks = (
        fractional_ranks([p["deepimmuno_score"] for p in eligible_pairs], True)
        if config.ranking != "binding"
        else []
    )
    for index, pair in enumerate(eligible_pairs):
        pair["binding_fractional_rank"] = binding_ranks[index]
        if config.ranking == "binding":
            pair["ranking_cost"] = pair["ic50_nm"]
        else:
            pair["immunogenicity_fractional_rank"] = immuno_ranks[index]
            pair["ranking_cost"] = (
                1 - pair["deepimmuno_score"]
                if config.ranking == "immunogenicity"
                else config.alpha * binding_ranks[index] + (1 - config.alpha) * immuno_ranks[index]
            )
        pair["normalization_population"] = len(eligible_pairs)
    by_peptide = defaultdict(list)
    for pair in pairs:
        by_peptide[pair["peptide"]].append(pair)
    ranked = []
    for epitope in epitopes:
        options = [p for p in by_peptide[epitope["peptide"]] if p["eligible"]]
        best = (
            min(
                options,
                key=lambda p: (p["ranking_cost"], p["ic50_nm"], p["peptide"], p["allele"], p["pair_id"]),
            )
            if options
            else None
        )
        complete = all(p["ic50_nm"] is not None for p in by_peptide[epitope["peptide"]])
        reference = predictions.get((epitope["reference_peptide"], best["allele"])) if best else None
        ranked.append(
            {
                **epitope,
                "ranking_cost": best["ranking_cost"] if best and complete else None,
                "eligible": best is not None and complete,
                "all_alleles_scored": complete,
                "best_allele": best["allele"] if best else None,
                "ic50_nm": best["ic50_nm"] if best else None,
                "deepimmuno_score": best["deepimmuno_score"] if best else None,
                "reference_ic50_nm": reference["ic50"] if reference else None,
                "reference_over_mutant_ratio": reference["ic50"] / best["ic50_nm"] if reference else None,
            }
        )
    ranked.sort(key=epitope_sort_key)
    # Sequence rank is shared across transcript/position associations.
    sequence_rank = {}
    for row in ranked:
        if row["eligible"]:
            sequence_rank.setdefault(row["peptide"], len(sequence_rank) + 1)
        row["epitope_rank"] = sequence_rank.get(row["peptide"])
    return pairs, ranked


def epitope_sort_key(row):
    return (
        row["ranking_cost"] if row["ranking_cost"] is not None else float("inf"),
        row["ic50_nm"] if row["ic50_nm"] is not None else float("inf"),
        row["peptide"],
        row["best_allele"] or "",
        row["epitope_id"],
    )


def summarize_mutations(ranked, aggregation):
    groups = defaultdict(dict)
    for epitope in ranked:
        groups[epitope["mutation_id"]].setdefault(epitope["peptide"], epitope)
    aggregate = {"min": min, "mean": statistics.mean, "median": statistics.median}[aggregation]
    rows = []
    for mutation_id, distinct in groups.items():
        costs = [r["ranking_cost"] for r in distinct.values() if r["eligible"]]
        best = min(distinct.values(), key=epitope_sort_key)
        rows.append(
            {
                "mutation_id": mutation_id,
                "gene": best["gene"],
                "protein_change": best["protein_change"],
                "aggregation": aggregation,
                "ranking_cost": aggregate(costs) if costs else None,
                "eligible_epitopes": len(costs),
                "total_epitopes": len(distinct),
                "best_epitope_id": best["epitope_id"] if costs else None,
            }
        )
    rows.sort(
        key=lambda r: (r["ranking_cost"] if r["ranking_cost"] is not None else float("inf"), r["mutation_id"])
    )
    return rows
