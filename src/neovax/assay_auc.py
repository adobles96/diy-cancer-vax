"""Threshold-free discrimination within the observed, reconstructable vaccine set."""

import math
from collections import defaultdict


def assay_labels(components, targets, observations):
    """One row per mutation; positives are permissive, negatives require measured evidence."""
    by_target = {r["mutation_id"]: r for r in targets}
    assays = defaultdict(list)
    for observation in observations:
        if (
            observation["timepoint"] == "week 16"
            and observation["assay_type"] == "peptide"
            and observation["cell_type"] in {"CD4", "CD8"}
        ):
            assays[observation["imp_id"]].append(observation)
    groups, excluded = defaultdict(list), []
    for component in components:
        mutation = component["mutation_id"]
        if mutation not in by_target:
            excluded.append({**component, "exclusion_reason": "No reconstructed, scored mutation target"})
            continue
        groups[mutation].append(component)
    labeled = []
    for mutation, group in sorted(groups.items()):
        evidence = [o for c in group for o in assays[c["imp_id"]]]
        positive = any(o["result"] == "positive" for o in evidence)
        all_negative = all(
            {o["cell_type"] for o in assays[c["imp_id"]]} == {"CD4", "CD8"}
            and all(o["result"] == "negative" for o in assays[c["imp_id"]])
            for c in group
        )
        row = {
            **by_target[mutation],
            "imp_ids": sorted(c["imp_id"] for c in group),
            "label": 1 if positive else 0 if all_negative else None,
            "week16_observations": evidence,
        }
        if row["label"] is None:
            excluded.append(
                {**row, "exclusion_reason": "No positive; incomplete or ambiguous negative evidence"}
            )
        else:
            labeled.append(row)
    return labeled, excluded


def roc_metrics(labels, costs):
    """Lower cost predicts positive; ROC groups ties and AUC gives tied pairs half credit."""
    if len(labels) != len(costs) or any(label not in {0, 1} for label in labels):
        raise ValueError("Expected equally sized binary labels and costs")
    if any(not math.isfinite(cost) for cost in costs):
        raise ValueError("ROC costs must be finite")
    positives = [cost for label, cost in zip(labels, costs) if label == 1]
    negatives = [cost for label, cost in zip(labels, costs) if label == 0]
    n_positive, n_negative = len(positives), len(negatives)
    denominator = n_positive * n_negative
    wins = sum(a < b for a in positives for b in negatives)
    ties = sum(a == b for a in positives for b in negatives)
    curve = []
    if denominator:
        grouped = defaultdict(list)
        for label, cost in zip(labels, costs):
            grouped[cost].append(label)
        tp, fp = 0, 0
        curve = [{"cost_threshold": None, "fpr": 0.0, "tpr": 0.0}]
        for cost, group in sorted(grouped.items()):
            tp += sum(group)
            fp += len(group) - sum(group)
            curve.append({"cost_threshold": cost, "fpr": fp / n_negative, "tpr": tp / n_positive})
    return {
        "n": len(labels),
        "positives": n_positive,
        "negatives": n_negative,
        "positive_negative_pairs": denominator,
        "correctly_ordered_pairs": wins,
        "tied_pairs": ties,
        "roc_auc": (wins + 0.5 * ties) / denominator if denominator else None,
        "undefined_reason": None if denominator else "Both positive and negative targets are required",
        "roc_curve": curve,
    }
