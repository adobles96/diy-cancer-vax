import datetime
from dataclasses import replace

import pytest

from neovax.config import Config
from neovax.data import hla_value, result_state
from neovax.evaluation import assay_metrics
from neovax.peptides import enumerate_epitopes, reconstruct
from neovax.ranking import fractional_ranks, rank_epitopes, summarize_mutations
from neovax.selection import long_bounds, select_vaccine


def variant(identifier="v1", mutation="m1", change="p.A15V", transcript="tx.1"):
    return {
        "variant_id": identifier,
        "mutation_id": mutation,
        "patient": 1,
        "gene": "GENE",
        "protein_change": change,
        "transcript": transcript,
        "consequence": "Missense_Mutation",
    }


def protein_records():
    return reconstruct([variant()], {"tx.1": "A" * 40})[0]


def scores(epitopes, alleles=("HLA-A*02:01",)):
    return {
        (r[field], a): {"ic50": float(index + 10), "percentile_rank": 0.5, "request_id": "fixture"}
        for index, r in enumerate(epitopes)
        for a in alleles
        for field in ["peptide", "reference_peptide"]
    }


@pytest.mark.parametrize("position,length", [(0, 9), (1, 10), (14, 9), (14, 10), (38, 10), (39, 9)])
def test_all_mutation_containing_windows(position, length):
    proteins, skipped = reconstruct([variant(change=f"p.A{position + 1}V")], {"tx.1": "A" * 40})
    rows = enumerate_epitopes(proteins, [length])
    expected = {start for start in range(41 - length) if start <= position < start + length}
    assert {r["start"] for r in rows} == expected
    assert not skipped
    for row in rows:
        assert row["peptide"][row["mutation_offset"]] == "V"
        assert row["reference_peptide"] == "A" * length


def test_exact_transcript_and_reference_checks():
    assert (
        reconstruct([variant(transcript="tx.2")], {"tx.1": "A" * 40})[1][0]["reason"]
        == "exact transcript version missing"
    )
    assert (
        reconstruct([variant(change="p.C15V")], {"tx.1": "A" * 40})[1][0]["reason"]
        == "reference amino acid mismatch"
    )
    assert (
        reconstruct([variant(change="p.A0V")], {"tx.1": "A" * 40})[1][0]["reason"]
        == "protein position out of bounds"
    )


def test_remote_ambiguous_residue_does_not_remove_valid_windows():
    proteins, skipped = reconstruct([variant()], {"tx.1": "A" * 39 + "X"})
    assert len(proteins) == 1 and not skipped
    assert len(enumerate_epitopes(proteins, (9, 10))) == 19
    proteins, _ = reconstruct([variant()], {"tx.1": "A" * 18 + "X" + "A" * 21})
    assert all("X" not in r["peptide"] for r in enumerate_epitopes(proteins, (9, 10)))


@pytest.mark.parametrize(
    "value,expected",
    [
        (datetime.time(2, 1), "02:01"),
        (datetime.timedelta(hours=44, minutes=2), "44:02"),
        (datetime.timedelta(hours=66, minutes=1), "66:01"),
        ("24:02", "24:02"),
    ],
)
def test_hla_time_conversion(value, expected):
    assert hla_value(value) == expected


@pytest.mark.parametrize(
    "value,state",
    [
        (None, "not reported"),
        ("n.d.", "not determined"),
        (0, "negative"),
        (1, "positive"),
        ("?", "ambiguous"),
    ],
)
def test_result_states(value, state):
    assert result_state(value) == state


def test_ranks_ties_and_score_directions():
    assert fractional_ranks([10, 10, 30]) == [0.25, 0.25, 1]
    assert fractional_ranks([1]) == [0]
    assert fractional_ranks([0.1, 0.8, 0.8], True) == [1, 0.25, 0.25]
    epitopes = enumerate_epitopes(protein_records(), [9])[:2]
    predictions = scores(epitopes)
    immuno = {(r["peptide"], "HLA-A*02:01"): 0.1 + index * 0.8 for index, r in enumerate(epitopes)}
    config = Config(ranking="immunogenicity", deepimmuno_manifest="fixture")
    _, ranked = rank_epitopes(epitopes, ["HLA-A*02:01"], predictions, config, immuno)
    assert ranked[0]["peptide"] == epitopes[1]["peptide"]
    pairs, _ = rank_epitopes(
        epitopes, ["HLA-A*02:01"], predictions, replace(config, ranking="combined"), immuno
    )
    assert [p["ranking_cost"] for p in pairs] == [0.5, 0.5]
    with pytest.raises(ValueError, match="Missing DeepImmuno"):
        rank_epitopes(epitopes, ["HLA-A*02:01"], predictions, config, {})


def test_missing_binding_and_strict_gate():
    ep = enumerate_epitopes(protein_records(), [9])[:1]
    pred = scores(ep)
    key = (ep[0]["peptide"], "HLA-A*02:01")
    pred[key]["ic50"] = 500
    assert not rank_epitopes(ep, ["HLA-A*02:01"], pred, Config())[1][0]["eligible"]
    assert rank_epitopes(ep, ["HLA-A*02:01"], pred, Config(threshold=None))[1][0]["eligible"]
    assert not rank_epitopes(ep, ["HLA-A*02:01", "HLA-B*07:02"], pred, Config(threshold=None))[1][0][
        "eligible"
    ]


def test_aggregation_deduplicates_and_only_changes_mutation_selection():
    base = enumerate_epitopes(protein_records(), [9])[:3]
    records = []
    for index, (mutation, cost) in enumerate([("m1", 1), ("m1", 99), ("m2", 30)]):
        records.append(
            {
                **base[index],
                "mutation_id": mutation,
                "ranking_cost": cost,
                "ic50_nm": cost,
                "eligible": True,
                "best_allele": "HLA-A*02:01",
                "all_alleles_scored": True,
                "deepimmuno_score": None,
            }
        )
    duplicate = {**records[0], "epitope_id": "duplicate"}
    minimum = summarize_mutations([*records, duplicate], "min")
    mean = summarize_mutations([*records, duplicate], "mean")
    median = summarize_mutations([*records, duplicate], "median")
    assert minimum[0]["mutation_id"] == "m1"
    assert mean[0]["mutation_id"] == median[0]["mutation_id"] == "m2"
    assert mean[1]["ranking_cost"] == 50 and mean[1]["total_epitopes"] == 2
    proteins = protein_records()
    assert select_vaccine(records, minimum, proteins, Config(budget=1))[0][0]["mutation_id"] == "m1"
    assert select_vaccine(records, mean, proteins, Config(budget=1))[0][0]["mutation_id"] == "m1"
    assert (
        select_vaccine(records, mean, proteins, Config(budget=1, selection_unit="mutation"))[0][0][
            "mutation_id"
        ]
        == "m2"
    )


@pytest.mark.parametrize(
    "start,end,n,width", [(0, 9, 40, 25), (31, 40, 40, 25), (10, 20, 40, 25), (2, 11, 18, 30)]
)
def test_long_boundaries(start, end, n, width):
    left, right = long_bounds(start, end, n, width)
    assert 0 <= left <= start < end <= right <= n
    assert right - left == min(n, width)


def test_short_proteins_and_width_fail():
    with pytest.raises(ValueError):
        long_bounds(0, 9, 14, 25)
    with pytest.raises(ValueError):
        long_bounds(0, 10, 40, 9)


def test_determinism_and_valid_long_containment():
    proteins = protein_records()
    epitopes = enumerate_epitopes(proteins, [9, 10])
    prediction = scores(epitopes)
    _, a = rank_epitopes(epitopes, ["HLA-A*02:01"], prediction, Config())
    _, b = rank_epitopes(list(reversed(epitopes)), ["HLA-A*02:01"], prediction, Config())
    assert a == b
    selected, summary = select_vaccine(a, summarize_mutations(a, "min"), proteins, Config())
    assert summary["shortfall"] == 19
    assert len(selected) == 1 and selected[0]["driving_epitope"] in selected[0]["sequence"]
    assert selected[0]["sequence"][selected[0]["mutation_offset"]] == "V"


def test_worst_gate_ties_do_not_inflate_average_precision():
    ranked = [
        {"peptide": "A", "all_alleles_scored": True, "ranking_cost": 1, "eligible": True, "epitope_rank": 1},
        {
            "peptide": "B",
            "all_alleles_scored": True,
            "ranking_cost": None,
            "eligible": False,
            "epitope_rank": None,
        },
        {
            "peptide": "C",
            "all_alleles_scored": True,
            "ranking_cost": None,
            "eligible": False,
            "epitope_rank": None,
        },
    ]
    assays = [
        {
            "sequence": seq,
            "cell_type": "CD8",
            "assay_type": "peptide",
            "timepoint": "week 16",
            "result": state,
        }
        for seq, state in [("A", "negative"), ("B", "positive"), ("C", "negative")]
    ]
    result = assay_metrics(ranked, assays, 2)
    assert result["average_precision"] == pytest.approx(1 / 3)
    assert result["precision_at_k"] == 0.25
    assert result["positive_ranks"][0]["evaluation_rank"] == 2.5
