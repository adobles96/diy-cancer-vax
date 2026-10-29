import random

import pytest

from neovax.deepimmuno import cached_prediction_subset
from neovax.io import write_json
from neovax.reference_filter import filter_reference_peptides


def candidate(peptide, mutation="m1", transcript="mutated.1"):
    return {"peptide": peptide, "mutation_id": mutation, "transcript": transcript}


def test_exact_full_length_and_cross_protein_boundary():
    reference = {"other.1": "XXACDEFGHIKM", "left.1": "YYQRST", "right.1": "VWYAA"}
    rows = [candidate(p) for p in ["ACDEFGHIK", "ACDEFGHIKL", "ACDEYGHIK", "QRSTVWYAA"]]
    kept, removed, summary = filter_reference_peptides(rows, reference)
    assert [r["peptide"] for r in kept] == ["ACDEFGHIKL", "ACDEYGHIK", "QRSTVWYAA"]
    assert removed[0]["reference_matches"] == [{"transcript": "other.1", "start": 2, "end": 11}]
    assert summary["unique_peptides_removed"] == 1
    assert rows[0] == candidate("ACDEFGHIK")  # Input records remain unchanged.


def test_all_associations_isoforms_and_source_protein_are_filtered():
    reference = {"mutated.1": "ACDEFGHIKACDEFGHIK", "isoform.2": "ACDEFGHIK"}
    rows = [candidate("ACDEFGHIK"), candidate("ACDEFGHIK", "m2"), candidate("ACDEFGHIL", "m2")]
    kept, removed, summary = filter_reference_peptides(rows, reference)
    assert len(kept) == 1 and len(removed) == 2
    assert summary["candidate_associations_removed"] == 2
    assert summary["unique_peptides_removed"] == 1
    assert summary["mutations_affected"] == 2
    assert summary["mutations_fully_removed"] == ["m1"]
    assert len(removed[0]["reference_matches"]) == 3
    assert removed[0]["reference_matches"] == removed[1]["reference_matches"]


@pytest.mark.parametrize("length", [8, 9, 10, 11])
def test_supported_lengths_and_overlaps(length):
    peptide = "A" * length
    kept, removed, _ = filter_reference_peptides([candidate(peptide)], {"tx": "A" * (length + 1)})
    assert not kept
    assert [m["start"] for m in removed[0]["reference_matches"]] == [0, 1]


def test_matcher_matches_independent_naive_search():
    rng = random.Random(472)
    reference = {f"tx{i}": "".join(rng.choices("ACDEFG", k=40)) for i in range(40)}
    peptides = {seq[i : i + 9] for seq in list(reference.values())[:5] for i in range(3)}
    peptides.update("".join(rng.choices("ACDEFGHIK", k=9)) for _ in range(40))
    rows = [candidate(p, p) for p in sorted(peptides)]
    _, removed, _ = filter_reference_peptides(rows, reference)
    actual = {r["peptide"]: r["reference_matches"] for r in removed}
    expected = {
        p: [
            {"transcript": tx, "start": i, "end": i + len(p)}
            for tx, seq in sorted(reference.items())
            for i in range(len(seq) - len(p) + 1)
            if seq[i : i + len(p)] == p
        ]
        for p in peptides
    }
    assert actual == {p: matches for p, matches in expected.items() if matches}
    assert filter_reference_peptides([], reference)[2]["unique_peptides_removed"] == 0


def test_subset_cache_reuse_rejects_different_models_and_rescues(tmp_path):
    manifest = {"commit": "pinned", "files": [{"sha256": "weights"}]}
    record = {"peptide": "ACDEFGHIK", "allele": "HLA-A*66:01", "score": 0.8}
    wanted = [(record["peptide"], record["allele"])]
    representations = {"HLA-A*66:01": "HLA-A*68:01"}
    saved = {"model": manifest, "predictions": [record, {**record, "peptide": "EXTRA"}]}
    write_json(tmp_path / "cache.json", saved)
    assert cached_prediction_subset(tmp_path, manifest, wanted, representations) is None
    saved["allele_representations"] = representations
    write_json(tmp_path / "cache.json", saved)
    assert cached_prediction_subset(tmp_path, manifest, wanted, representations) == [record]
    assert (
        cached_prediction_subset(tmp_path, {**manifest, "commit": "other"}, wanted, representations) is None
    )
    assert cached_prediction_subset(tmp_path, manifest, [], representations) == []
