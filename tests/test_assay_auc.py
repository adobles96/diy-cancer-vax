import math

import pytest
from scipy.stats import rankdata

from neovax.assay_auc import assay_labels, roc_metrics


def test_auc_direction_ties_and_undefined():
    assert roc_metrics([1, 1, 0, 0], [0, 1, 2, 3])["roc_auc"] == 1
    assert roc_metrics([1, 1, 0, 0], [3, 2, 1, 0])["roc_auc"] == 0
    result = roc_metrics([1, 1, 0, 0], [1, 1, 1, 1])
    assert result["roc_auc"] == 0.5
    assert result["tied_pairs"] == 4
    assert len(result["roc_curve"]) == 2
    assert roc_metrics([1, 1], [0, 1])["roc_auc"] is None
    assert roc_metrics([], [])["roc_auc"] is None
    with pytest.raises(ValueError):
        roc_metrics([1], [math.nan])
    with pytest.raises(ValueError):
        roc_metrics([1, 0], [1])


def test_auc_matches_independent_rank_sum_and_curve_area():
    labels = [1, 0, 1, 0, 0, 1, 1]
    costs = [5, 2, 2, 1, 4, 3, 2]
    result = roc_metrics(labels, costs)
    ranks = rankdata([-cost for cost in costs])
    n_positive, n_negative = sum(labels), len(labels) - sum(labels)
    expected = (
        sum(rank for rank, label in zip(ranks, labels) if label) - n_positive * (n_positive + 1) / 2
    ) / (n_positive * n_negative)
    assert result["roc_auc"] == pytest.approx(expected)
    points = result["roc_curve"]
    area = sum(
        (right["fpr"] - left["fpr"]) * (right["tpr"] + left["tpr"]) / 2
        for left, right in zip(points, points[1:])
    )
    assert area == pytest.approx(expected)


def observation(imp_id, cell, result, **changes):
    return {
        "imp_id": imp_id,
        "cell_type": cell,
        "result": result,
        "assay_type": "peptide",
        "timepoint": "week 16",
        **changes,
    }


def test_labels_deduplicate_mutations_and_keep_unknowns_out():
    components = [
        {"mutation_id": mutation, "imp_id": imp}
        for mutation, imp in [("m1", "a"), ("m1", "b"), ("m2", "c"), ("m3", "d"), (None, "e")]
    ]
    targets = [{"mutation_id": mutation} for mutation in ["m1", "m2", "m3"]]
    observations = [
        observation("a", "CD8", "negative"),
        observation("b", "CD4", "positive"),
        observation("c", "CD4", "negative"),
        observation("c", "CD8", "negative"),
        observation("c", "CD4", "positive", timepoint="year 3–4.5"),
        observation("c", "CD8", "positive", assay_type="tumor"),
        observation("d", "CD4", "negative"),
        observation("d", "CD8", "not determined"),
        observation("e", "CD4", "positive"),
    ]
    labeled, excluded = assay_labels(components, targets, observations)
    assert [(r["mutation_id"], r["label"]) for r in labeled] == [("m1", 1), ("m2", 0)]
    assert labeled[0]["imp_ids"] == ["a", "b"]
    assert len(excluded) == 2
    assert {r["mutation_id"] for r in excluded} == {None, "m3"}
    assert len(assay_labels(components, targets, [])[0]) == 0
