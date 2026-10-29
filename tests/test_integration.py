import copy
import csv
import json
from pathlib import Path

import pytest

from neovax import pipeline
from neovax.config import Config
from neovax.data import actual_vaccine, load_reference, patient_data
from neovax.deepimmuno import validate_model
from neovax.io import digest, write_json
from neovax.peptides import reconstruct
from neovax.report import render_report

DATA = Path("data/raw")
CACHE = Path("data/outputs/cache")
pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not (DATA / "hu2021_datasets.xlsx").exists(), reason="Fetch published source workbooks first"
    ),
]


@pytest.mark.parametrize(
    "patient,counts",
    [
        (1, (179, 84, 79, 5)),
        (2, (1887, 960, 899, 61)),
        (3, (820, 395, 356, 39)),
        (4, (1995, 1036, 960, 76)),
        (5, (728, 350, 325, 25)),
        (6, (3920, 1864, 1723, 141)),
    ],
)
def test_all_six_reconstruction_counts(patient, counts):
    variants, alleles, published, assays = patient_data(str(DATA), patient)
    reconstructed, skipped = reconstruct(variants, load_reference(str(DATA)))
    assert (
        len(variants),
        sum(v["consequence"] == "Missense_Mutation" for v in variants),
        len(reconstructed),
        sum(v["reason"] == "exact transcript version missing" for v in skipped),
    ) == counts
    assert len(alleles) == len(set(alleles)) and all(a.startswith(("HLA-A", "HLA-B")) for a in alleles)
    assert all(a["sources"] and all(s["cell"] and s["sha256"] for s in a["sources"]) for a in assays)
    assert len({a["observation_id"] for a in assays}) == len(assays)


def test_patient1_assay_regression():
    variants, _, published, assays = patient_data(str(DATA), 1)
    rows = [
        r
        for r in assays
        if r["cell_type"] == "CD8"
        and r["assay_type"] == "peptide"
        and r["timepoint"] == "week 16"
        and r["result"] in {"positive", "negative"}
    ]
    assert len(rows) == 23
    assert sum(r["result"] == "positive" for r in rows) == 2
    for sequence in ["KLKFVTLVF", "RFLEYLPLRF"]:
        response = {
            r["assay_type"]: r["result"]
            for r in assays
            if r["sequence"] == sequence and r["cell_type"] == "CD8" and r["timepoint"] == "week 16"
        }
        assert response == {"peptide": "positive", "minigene": "positive", "tumor": "negative"}
    assert all(
        r["result"] == "not determined"
        for r in assays
        if r["gene"] in {"CASP5a", "CASP5b"} and r["cell_type"] == "CD8"
    )
    proteins, _ = reconstruct(variants, load_reference(str(DATA)))
    vaccine = actual_vaccine(published, variants, proteins)
    assert sum(v["reconstructable"] for v in vaccine) == 11


@pytest.fixture
def legacy_run(tmp_path):
    if not (CACHE / "binding.sqlite").exists():
        pytest.skip("Import the verified Patient 1 cache first")
    config = Config(
        patient=1,
        threshold=None,
        budget=11,
        selection_unit="mutation",
        cached_only=True,
        output_dir=str(tmp_path),
    )
    return pipeline.run_pipeline(config, make_report=False)


def test_legacy_patient1(legacy_run):
    assert legacy_run["summary"]["status"] == "complete"
    assert legacy_run["summary"]["overlap"]["target_overlap"] == 4
    assert legacy_run["summary"]["coverage"]["reconstructed_rows"] == 79
    assert legacy_run["evaluations"]["score_correlations"]["mutant"]["n"] == 25


def test_label_leakage(legacy_run, tmp_path, monkeypatch):
    actual_reader = pipeline.patient_data

    def corrupted(directory, patient):
        variants, alleles, published, assays = copy.deepcopy(actual_reader(directory, patient))
        for row in published:
            row["gene_tpm"] = 999999
            row["protein_change"] = "p.A999999V"
        for row in assays:
            row["result"] = "negative" if row["result"] == "positive" else "positive"
        return variants, alleles, published, assays

    monkeypatch.setattr(pipeline, "patient_data", corrupted)
    config = Config(**{**legacy_run["config"], "output_dir": str(tmp_path)})
    changed = pipeline.run_pipeline(config, make_report=False)
    assert changed["pairs"] == legacy_run["pairs"]
    assert [r["component_id"] for r in changed["proposed"]] == [
        r["component_id"] for r in legacy_run["proposed"]
    ]
    assert changed["summary"]["overlap"]["target_overlap"] == 0


def test_deepimmuno_exact_allele_coverage():
    path = Path("data/outputs/models/deepimmuno/manifest.json")
    if not path.exists():
        pytest.skip("Optional model not fetched")
    for patient in range(1, 7):
        alleles = patient_data(str(DATA), patient)[1]
        if patient == 1:
            validate_model(path, alleles, (9, 10))
        else:
            with pytest.raises(ValueError, match="Unsupported alleles"):
                validate_model(path, alleles, (9, 10))


def test_failed_optional_model_preserves_binding_without_combined_metrics(tmp_path):
    path = Path("data/outputs/models/deepimmuno/manifest.json")
    if not path.exists():
        pytest.skip("Optional model not fetched")
    result = pipeline.run_pipeline(
        Config(
            patient=2,
            ranking="combined",
            deepimmuno_manifest=str(path),
            cached_only=True,
            output_dir=str(tmp_path),
        ),
        make_report=False,
    )
    assert result["summary"]["status"] == "partial"
    assert not result["proposed"] and result["pairs"][0]["ic50_nm"] is not None
    assert all(row["ranking_cost"] is None for row in result["pairs"])
    assert result["evaluations"]["cd8"]["intersection_count"] == 0
    assert result["evaluations"]["cd8"]["average_precision"] is None
    root = Path(result["summary"]["run_directory"])
    with (root / "proposed_vaccine.csv").open() as handle:
        shortlist = csv.DictReader(handle)
        assert "sequence" in shortlist.fieldnames and "ranking_cost" in shortlist.fieldnames
        assert list(shortlist) == []
    assert {p.name for p in root.iterdir()} == {"proposed_vaccine.csv", "run.json"}


def test_patient2_pipeline_uses_audited_publication_suffixes(tmp_path, monkeypatch):
    forwarded = []

    def unavailable_model(epitopes, alleles, config, cancel, progress, *, allow_hla_rescue):
        forwarded.append(allow_hla_rescue)
        raise ValueError("Test model unavailable")

    monkeypatch.setattr(pipeline, "predict_immunogenicity", unavailable_model)
    result = pipeline.run_pipeline(
        Config(
            patient=2,
            ranking="combined",
            deepimmuno_manifest="fixture",
            rescue_unknown_hla=True,
            cached_only=True,
            output_dir=str(tmp_path),
        ),
        make_report=False,
    )
    assert forwarded == [True]
    assert result["summary"]["overlap"]["reconstructable_historical_mutations"] == 10
    repaired = [r for r in result["actual"] if r["gene"] in {"C14orf101a", "C14orf101b"}]
    assert len(repaired) == 2
    assert all(r["reconstructable"] and r["mutation_id"] for r in repaired)
    assert all("full mutant long-peptide sequence verified" in r["target_join"] for r in repaired)


def test_legacy_report_regeneration(legacy_run):
    root = Path(legacy_run["summary"]["run_directory"])
    (root / "run.json").rename(root / "run_data.json")
    write_json(root / "manifest.json", legacy_run["manifest"])
    report = render_report(root)
    saved = json.loads((root / "manifest.json").read_text())
    assert saved["exports"]["report.html"] == digest(report)
    assert not (root / "figures").exists()


def test_failed_run_retains_one_diagnostic_record(tmp_path, monkeypatch):
    def fail(*args):
        raise ValueError("Source data is invalid")

    monkeypatch.setattr(pipeline, "execute", fail)
    with pytest.raises(ValueError, match="Source data"):
        pipeline.run_pipeline(Config(output_dir=str(tmp_path)))
    root = next(tmp_path.iterdir())
    assert [p.name for p in root.iterdir()] == ["run.json"]
    saved = json.loads((root / "run.json").read_text())
    assert saved["manifest"]["status"] == "failed"
    assert saved["config"]["patient"] == 1
    with pytest.raises(ValueError, match="did not finish"):
        render_report(root)
