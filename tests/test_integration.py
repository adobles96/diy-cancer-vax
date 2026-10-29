import copy
import json
from pathlib import Path

import pytest

from neovax import pipeline
from neovax.config import Config
from neovax.data import actual_vaccine, load_reference, patient_data
from neovax.deepimmuno import validate_model
from neovax.peptides import reconstruct

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


@pytest.mark.parametrize("patient", range(1, 7))
def test_six_complete_export_runs(patient):
    summary_path = Path("data/outputs/six-patient-validation.json")
    if not summary_path.exists():
        pytest.skip("Six-patient validation run has not been created")
    summaries = json.loads(summary_path.read_text())
    summary = next((s for s in summaries if s["patient"] == patient), None)
    assert summary is not None and summary["status"] == "complete"
    assert summary["coverage"]["missing_unique_predictions"] == 0
    root = Path(summary["run_directory"])
    for filename in [
        "config.json",
        "manifest.json",
        "summary.json",
        "variants.csv",
        "skipped_variants.csv",
        "epitopes.csv",
        "pair_predictions.csv",
        "mutation_scores.csv",
        "proposed_vaccine.csv",
        "proposed_vaccine.fasta",
        "actual_vaccine.csv",
        "assays.csv",
        "evaluations.json",
        "report.html",
        "figures/shortlist.svg",
        "figures/shortlist.png",
    ]:
        assert (root / filename).exists()
