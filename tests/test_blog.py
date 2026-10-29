"""Protect benchmark denominators and the reproducible six-patient workflow."""

import csv
from pathlib import Path

import pytest

from neovax.blog import pooled_auc_metrics, predict_retained_binding, recognition_tables, reproduce_blog
from neovax.config import Config
from neovax.iedb import Cancelled, InferenceError


def test_binding_requests_only_retained_exact_length_candidates_and_actual_hla(monkeypatch):
    calls = []
    alleles = ["HLA-A*01:01", "HLA-B*56:01"]
    epitopes = [{"peptide": "A" * 9}, {"peptide": "V" * 10}, {"peptide": "A" * 9}]

    def predict(cache, sequences, actual, lengths, **kwargs):
        calls.append((sequences, actual, lengths, kwargs["cached_only"]))
        scores = {(sequence, allele): {"ic50": 10} for sequence in sequences for allele in actual}
        scores[("UNRELATED", actual[0])] = {"ic50": 1}
        return scores, {"remaining_requests": 0}, "complete", []

    monkeypatch.setattr("neovax.blog.predict", predict)
    binding, statistics = predict_retained_binding(
        None, epitopes, alleles, Config(), lambda _: None, lambda: False
    )
    assert calls == [(["A" * 9], alleles, (9,), False), (["V" * 10], alleles, (10,), False)]
    assert set(binding) == {(row["peptide"], allele) for row in epitopes for allele in alleles}
    assert [row["length"] for row in statistics] == [9, 10]


def test_partial_binding_cannot_produce_a_misleading_cohort(monkeypatch):
    monkeypatch.setattr(
        "neovax.blog.predict", lambda *args, **kwargs: ({}, {"remaining_requests": 3}, "partial", [])
    )
    with pytest.raises(InferenceError, match="3 binding requests are not cached"):
        predict_retained_binding(
            None,
            [{"peptide": "A" * 9}],
            ["HLA-A*01:01"],
            Config(cached_only=True),
            lambda _: None,
            lambda: False,
        )
    with pytest.raises(Cancelled):
        predict_retained_binding(None, [], [], Config(), lambda _: None, lambda: True)


def test_tumor_denominator_uses_tested_components_and_week16_only():
    components = [
        {"patient": 1, "imp_id": "positive", "immunogenic": True, "netmhcpan_rank": 1},
        {"patient": 1, "imp_id": "negative", "immunogenic": False, "netmhcpan_rank": 2},
        {"patient": 1, "imp_id": "missing", "immunogenic": False, "netmhcpan_rank": 3},
    ]
    base = {"patient": 1, "assay_type": "tumor", "timepoint": "week 16", "cell_type": "CD8"}
    observations = [
        {**base, "imp_id": "positive", "result": "positive"},
        {**base, "imp_id": "positive", "result": "positive", "cell_type": "CD4"},
        {**base, "imp_id": "negative", "result": "negative"},
        {**base, "imp_id": "missing", "result": "not determined"},
        {**base, "imp_id": "missing", "result": "positive", "timepoint": "year 3–4.5"},
    ]
    peptide, tumor, _ = recognition_tables(components, observations)
    total = next(row for row in tumor if row["patient"] == "Total" and row["method"] == "netmhcpan")
    assert total == {
        "patient": "Total",
        "method": "netmhcpan",
        "positive": 1,
        "tested": 2,
        "included": 3,
        "untested": 1,
        "fraction": 0.5,
    }
    empty = next(row for row in tumor if row["patient"] == "Total" and row["method"] == "deepimmuno")
    assert empty["fraction"] is None
    assert (
        next(row for row in peptide if row["patient"] == "Total" and row["method"] == "neovax")["included"]
        == 3
    )


@pytest.mark.integration
def test_cached_real_blog_reproduces_saved_results(tmp_path):
    manifest = Path("data/outputs/models/deepimmuno/manifest.json")
    if not manifest.exists() or not Path("data/raw/ott2017_table5.xlsx").exists():
        pytest.skip("Downloaded sources and pinned model are required")
    if not Path("data/outputs/cache/deepimmuno").exists():
        pytest.skip("The six-patient DeepImmuno cache is required")
    output = tmp_path / "blog"
    summary = reproduce_blog(Config(cached_only=True, deepimmuno_manifest=str(manifest)), output)
    assert [row["hits"] for row in summary["vaccine_recall"]] == [16, 6, 19]
    assert [row["denominator"] for row in summary["vaccine_recall"]] == [98, 98, 98]
    assert [(row["positive"], row["included"]) for row in summary["peptide_recognition"]] == [
        (63, 98),
        (14, 16),
        (5, 6),
        (11, 19),
    ]
    assert [(row["positive"], row["tested"]) for row in summary["tumor_recognition"]] == [
        (3, 31),
        (2, 8),
        (0, 2),
        (0, 7),
    ]
    assert (summary["labeled_targets"], summary["positive_labels"], summary["negative_labels"]) == (
        72,
        45,
        27,
    )
    assert [row["pooled_auc"] for row in summary["auc"]] == pytest.approx(
        [0.6477366255144033, 0.5012345679012346, 0.5539094650205761]
    )
    assert all(row["positive_negative_pairs"] == 1215 for row in summary["auc"])
    with (output / "roc-auc-table.csv").open() as handle:
        auc_rows = list(csv.DictReader(handle))
    assert [row["Patient"] for row in auc_rows] == ["1", "2", "3", "4", "5", "6", "Overall (pooled)"]
    for path in [
        output / "report.html",
        output / "cohort-auc/report.html",
        output / "peptide-recognition-roc-auc.svg",
    ]:
        assert "pair-weighted" not in path.read_text()
    assert summary["allow_hla_rescue"]
    assert any(summary["rescued_alleles"].values())
    for name in ("vaccine-recall-at-20", "peptide-recognition-roc-auc"):
        for extension in ("png", "svg", "pdf"):
            assert (output / f"{name}.{extension}").stat().st_size > 1000
    for name in ("peptide-recognition-table.csv", "tumor-recognition-tested-table.csv", "roc-auc-table.csv"):
        assert name in (output / "report.html").read_text()
    baseline = Path("data/outputs/blog/peptide-recognition-table.csv")
    if baseline.exists():
        with baseline.open() as handle:
            saved = list(csv.DictReader(handle))
        with (output / "peptide-recognition-table.csv").open() as handle:
            regenerated = list(csv.DictReader(handle))
        assert regenerated == saved


def test_pooled_auc_uses_raw_scores_and_cross_patient_pairs():
    # Both patients rank their positive first, but their raw score scales differ.
    # All four positive-negative comparisons must count, including across patients.
    rows = [
        {"patient": 1, "label": 1, "binding_ic50_nm": 10, "deepimmuno_score": 0.4, "mean_ensemble_cost": 0.1},
        {"patient": 1, "label": 0, "binding_ic50_nm": 20, "deepimmuno_score": 0.3, "mean_ensemble_cost": 0.4},
        {
            "patient": 2,
            "label": 1,
            "binding_ic50_nm": 100,
            "deepimmuno_score": 0.9,
            "mean_ensemble_cost": 0.4,
        },
        {
            "patient": 2,
            "label": 0,
            "binding_ic50_nm": 200,
            "deepimmuno_score": 0.8,
            "mean_ensemble_cost": 0.8,
        },
    ]
    aggregate = pooled_auc_metrics(rows)
    assert [r["pooled_auc"] for r in aggregate] == [0.75, 0.75, 0.875]
    assert all(r["positive_negative_pairs"] == 4 for r in aggregate)
    assert aggregate[2]["tied_pairs"] == 1
