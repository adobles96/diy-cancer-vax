import math
from collections import Counter, defaultdict

from scipy.stats import spearmanr


def annotate_membership(rows, vaccine, published, assays=()):
    by_target = defaultdict(list)
    for component in vaccine:
        if component["mutation_id"]:
            by_target[component["mutation_id"]].append(component)
    assayed = {
        (p["imp_id"], p["sequence"])
        for p in assays
        if p["cell_type"] == "CD8"
        and p["assay_type"] == "peptide"
        and p["result"] in {"positive", "negative"}
    }
    result = []
    for row in rows:
        matches = by_target[row["mutation_id"]]
        peptide = row.get("driving_epitope", row.get("peptide"))
        containment = [v for v in matches if peptide in v["immunizing_sequence"]]
        result.append(
            {
                **row,
                "historical_target": bool(matches),
                "historical_imp_ids": [v["imp_id"] for v in matches],
                "short_contained_in_historical": bool(containment),
                "exact_short_assayed": any((v["imp_id"], peptide) in assayed for v in matches),
                "exact_long_match": any(row.get("sequence") == v["immunizing_sequence"] for v in matches),
            }
        )
    return result


def assay_matrix(assays):
    rows = {}
    for observation in assays:
        if observation["assay_type"] not in {"peptide", "minigene", "tumor"}:
            continue
        key = (
            observation["sequence"],
            observation["cell_type"],
            observation["timepoint"],
            observation["pre_stimulation"],
        )
        if key not in rows:
            rows[key] = {
                k: observation[k]
                for k in [
                    "patient",
                    "gene",
                    "sequence",
                    "imp_id",
                    "cell_type",
                    "timepoint",
                    "pre_stimulation",
                    "treatment_context",
                ]
            }
            rows[key].update(
                peptide="not reported", minigene="not reported", tumor="not reported", sources=[]
            )
        rows[key][observation["assay_type"]] = observation["result"]
        rows[key]["sources"].extend(observation["sources"])
    return list(rows.values())


def assay_metrics(ranked, assays, budget):
    best = {}
    for row in ranked:
        best.setdefault(row["peptide"], row)
    observed = [
        r
        for r in assays
        if r["cell_type"] == "CD8"
        and r["assay_type"] == "peptide"
        and r["timepoint"] == "week 16"
        and r["result"] in {"positive", "negative"}
    ]
    tested = {r["sequence"]: r for r in observed}
    scored, missing = [], []
    for sequence, assay in tested.items():
        prediction = best.get(sequence)
        if (
            not prediction
            or not prediction["all_alleles_scored"]
            or not prediction.get("ranking_available", True)
        ):
            missing.append(
                {
                    "sequence": sequence,
                    "result": assay["result"],
                    "reason": "not independently reconstructed"
                    if not prediction
                    else "requested ranking unavailable"
                    if not prediction.get("ranking_available", True)
                    else "incomplete binding scores",
                }
            )
        else:
            scored.append(
                {
                    "sequence": sequence,
                    "result": assay["result"],
                    "positive": assay["result"] == "positive",
                    "cost": prediction["ranking_cost"],
                    "discovery_rank": prediction["epitope_rank"],
                    "binding_gate_failed": not prediction["eligible"],
                }
            )
    scored.sort(key=lambda r: (r["cost"] if r["cost"] is not None else float("inf"), r["sequence"]))
    # Threshold-group AP is invariant to ordering of tied scores, including tied gate failures.
    groups = defaultdict(list)
    for row in scored:
        groups[row["cost"]].append(row)
    positives = sum(r["positive"] for r in scored)
    total, tp, ap, curve = 0, 0, 0.0, []
    for group in groups.values():
        group_tp = sum(r["positive"] for r in group)
        first_rank = total + 1
        total += len(group)
        tp += group_tp
        average_rank = (first_rank + total) / 2
        for row in group:
            row["evaluation_rank"] = average_rank
        if positives:
            ap += (group_tp / positives) * tp / total
        curve.append(
            {"recall": tp / positives if positives else None, "precision": tp / total, "count": total}
        )
    # A boundary tie contributes its expected positive fraction to precision at K.
    k = min(budget, len(scored))
    remaining, expected_tp = k, 0.0
    for group in groups.values():
        included = min(remaining, len(group))
        expected_tp += included * sum(r["positive"] for r in group) / len(group)
        remaining -= included
    return {
        "full_tested_denominator": len(tested),
        "intersection_count": len(scored),
        "positives": positives,
        "negatives": len(scored) - positives,
        "full_tested_class_counts": dict(Counter(r["result"] for r in tested.values())),
        "average_precision": ap if positives else None,
        "precision_at_k": expected_tp / k if k else None,
        "k": k,
        "positive_ranks": [r for r in scored if r["positive"]],
        "tested_ranking": scored,
        "missing_tested": missing,
        "precision_recall": curve,
        "rule": "Best eligible allele per peptide from discovery scores; binding-gate failures tied last. AP groups ties; precision at K averages boundary ties. Untested candidates are not negatives.",
    }


def evaluate(selected, ranked, vaccine, published, assays, predictions, budget):
    historical_targets = {r["mutation_id"] for r in vaccine if r["mutation_id"]}
    reconstructable = {r["mutation_id"] for r in vaccine if r["reconstructable"]}
    selected_targets = {r["mutation_id"] for r in selected}
    comparisons, seen = [], set()
    for row in published:
        for state, sequence_field, score_field in [
            ("mutant", "mutant_epitope", "historical_mutant_ic50_nm"),
            ("reference", "wildtype_epitope", "historical_wildtype_ic50_nm"),
        ]:
            sequence, allele = row[sequence_field], row["hla"]
            key = (sequence, allele, state)
            new = predictions.get((sequence, allele))
            old = row[score_field]
            if key in seen or not new or not isinstance(old, (int, float)):
                continue
            seen.add(key)
            comparisons.append(
                {
                    "sequence": sequence,
                    "allele": allele,
                    "state": state,
                    "gene": row["gene"],
                    "historical_2_4_ic50_nm": old,
                    "new_4_1_ic50_nm": new["ic50"],
                    "source": row["source"],
                }
            )
    correlations = {}
    for state in ["mutant", "reference"]:
        group = [r for r in comparisons if r["state"] == state]
        rho = (
            float(
                spearmanr(
                    [r["historical_2_4_ic50_nm"] for r in group], [r["new_4_1_ic50_nm"] for r in group]
                ).statistic
            )
            if len(group) > 1
            else None
        )
        correlations[state] = {
            "n": len(group),
            "spearman_rho": rho if rho is not None and math.isfinite(rho) else None,
        }
    tumor = [r for r in assays if r["assay_type"] == "tumor"]
    return {
        "interpretation": "Exploratory retrospective comparison; model training overlap is not fully audited",
        "overlap": {
            "budget": budget,
            "selected_components": len(selected),
            "historical_components": len(vaccine),
            "audited_historical_mutations": len(historical_targets),
            "reconstructable_historical_mutations": len(reconstructable),
            "target_overlap": len(selected_targets & historical_targets),
            "short_containment_count": sum(r["short_contained_in_historical"] for r in selected),
            "exact_long_matches": sum(r["exact_long_match"] for r in selected),
            "excluded_historical_components": [r for r in vaccine if not r["reconstructable"]],
        },
        "cd8": assay_metrics(ranked, assays, budget),
        "historical_scores": comparisons,
        "score_correlations": correlations,
        "tumor_evidence": tumor,
        "tumor_tested_count": sum(r["result"] in {"positive", "negative"} for r in tumor),
        "tumor_caveat": "Selected reactive cell lines, not a complete candidate negative-label set. IFN-γ ELISpot measures secretion, not direct killing.",
    }


def rescore_historical(vaccine, predictions, alleles, lengths):
    rows = []
    for component in vaccine:
        sequence = component["immunizing_sequence"]
        pairs = {
            (sequence[start : start + length], allele)
            for length in lengths
            for start in range(len(sequence) - length + 1)
            for allele in alleles
        }
        scored = [predictions[key] for key in pairs if key in predictions]
        best = min(scored, key=lambda p: (p["ic50"], p["peptide"], p["allele"])) if scored else None
        rows.append(
            {
                "imp_id": component["imp_id"],
                "gene": component["gene"],
                "protein_change": component["protein_change"],
                "scope": "Historical sequence rescore only; all windows, not discovery candidates",
                "expected_pairs": len(set(pairs)),
                "scored_pairs": len(scored),
                "best_epitope": best["peptide"] if best else None,
                "ic50_nm": best["ic50"] if best else None,
                "allele": best["allele"] if best else None,
            }
        )
    return rows
