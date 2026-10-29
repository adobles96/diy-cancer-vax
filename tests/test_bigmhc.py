import csv
import importlib.util
import json
from pathlib import Path

import pytest

from neovax.bigmhc import COMMIT, exact_alleles, predict_immunogenicity, read_predictions

SPEC = importlib.util.spec_from_file_location(
    "bigmhc_pilot", Path(__file__).parents[1] / "scripts/bigmhc_pilot.py"
)
PILOT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PILOT)


def write_predictions(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mhc", "pep", "BigMHC_IM"])
        writer.writerows(rows)


def test_predictions_join_by_pair_not_output_order(tmp_path):
    path = tmp_path / "predictions.csv"
    pairs = [("AAAAAAAAA", "HLA-A*02:01"), ("CCCCCCCCC", "HLA-B*07:02")]
    write_predictions(path, [(pairs[1][1], pairs[1][0], 0.8), (pairs[0][1], pairs[0][0], 0.2)])
    assert read_predictions(path, pairs) == {pairs[0]: 0.2, pairs[1]: 0.8}


@pytest.mark.parametrize(
    "rows",
    [
        [],
        [("HLA-A*02:01", "AAAAAAAAA", 0.2), ("HLA-A*02:01", "AAAAAAAAA", 0.2)],
        [("HLA-A*02:02", "AAAAAAAAA", 0.2)],
        [("HLA-A*02:01", "AAAAAAAAA", "nan")],
        [("HLA-A*02:01", "AAAAAAAAA", 1.01)],
    ],
)
def test_predictions_reject_incomplete_duplicate_remapped_and_invalid(tmp_path, rows):
    path = tmp_path / "predictions.csv"
    write_predictions(path, rows)
    with pytest.raises(ValueError):
        read_predictions(path, [("AAAAAAAAA", "HLA-A*02:01")])


def test_exact_hla_disallows_upstream_fuzzy_rescue(tmp_path):
    (tmp_path / "data").mkdir()
    (tmp_path / "data/pseudoseqs.csv").write_text("mhc,1_A\nHLA-A*02:01,1\n")
    assert exact_alleles(tmp_path, ["HLA-A*02:01"]) == {"HLA-A*02:01": "HLA-A*02:01"}
    with pytest.raises(ValueError, match="No exact"):
        exact_alleles(tmp_path, ["HLA-A*02:02"])


def test_runtime_package_change_invalidates_cached_predictions(tmp_path, monkeypatch):
    root, output = tmp_path / "model", tmp_path / "result"
    (root / "data").mkdir(parents=True)
    (root / "src").mkdir()
    (root / "src/predict.py").write_text("# synthetic runner\n")
    (root / "data/pseudoseqs.csv").write_text("mhc,1_A\nHLA-A*02:01,1\n")
    (root / "manifest.json").write_text(json.dumps({"commit": COMMIT, "model": "BigMHC_IM", "files": []}))
    runtime = {"packages": {"torch": "2.2.2"}}
    runs = []

    def fake_run(command, **kwargs):
        runs.append(command)
        write_predictions(Path(command[command.index("-o") + 1]), [("HLA-A*02:01", "AAAAAAAAA", 0.5)])

    monkeypatch.setattr("neovax.bigmhc.subprocess.check_output", lambda *args, **kwargs: json.dumps(runtime))
    monkeypatch.setattr("neovax.bigmhc.subprocess.run", fake_run)
    pairs = [("AAAAAAAAA", "HLA-A*02:01")]
    predict_immunogenicity(pairs, root, "/synthetic/python", output)
    predict_immunogenicity(pairs, root, "/synthetic/python", output)
    assert len(runs) == 1
    runtime["packages"]["torch"] = "2.3.0"
    predict_immunogenicity(pairs, root, "/synthetic/python", output)
    assert len(runs) == 2


def test_pilot_aggregation_max_across_all_pairs_and_average_rank_ties():
    alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    epitopes = [
        {"mutation_id": m, "peptide": p, "gene": m, "protein_change": "p.A1C"}
        for m, p in [("m1", "AAAAAAAAA"), ("m1", "CCCCCCCCC"), ("m2", "DDDDDDDDD"), ("m3", "EEEEEEEEE")]
    ]
    pairs = {(e["peptide"], a) for e in epitopes for a in alleles}
    binding = {p: {"ic50": 100.0} for p in pairs}
    binding[("AAAAAAAAA", alleles[0])] = {"ic50": 1.0}
    binding[("DDDDDDDDD", alleles[0])] = {"ic50": 10.0}
    predictions = {p: 0.1 for p in pairs}
    predictions[("CCCCCCCCC", alleles[1])] = 0.9
    predictions[("DDDDDDDDD", alleles[0])] = 0.9
    rows = PILOT.bigmhc_targets(epitopes, alleles, binding, predictions)
    by_id = {r["mutation_id"]: r for r in rows}
    assert by_id["m1"]["binding_peptide"] == "AAAAAAAAA"
    assert by_id["m1"]["bigmhc_peptide"] == "CCCCCCCCC"
    assert by_id["m1"]["bigmhc_actual_allele"] == alleles[1]
    assert by_id["m1"]["bigmhc_score"] == 0.9
    assert by_id["m1"]["bigmhc_cost"] == by_id["m2"]["bigmhc_cost"] == 0.25
    assert by_id["m1"]["bigmhc_ensemble_cost"] == 0.125
    assert by_id["m2"]["bigmhc_ensemble_cost"] == 0.375
    assert by_id["m3"]["bigmhc_ensemble_rank"] == 3
    del predictions[("AAAAAAAAA", alleles[0])]
    with pytest.raises(ValueError, match="Incomplete model scores"):
        PILOT.bigmhc_targets(epitopes, alleles, binding, predictions)
