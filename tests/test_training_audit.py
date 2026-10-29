"""The NetMHCpan audit must work independently of any discarded model pilot."""

import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "audit_netmhcpan_overlap", Path(__file__).parents[1] / "scripts/audit_netmhcpan_overlap.py"
)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


@pytest.mark.integration
def test_candidate_context_from_current_blog_without_external_pairs():
    baseline = Path("data/outputs/blog/cohort-recall")
    labels_file = Path("data/outputs/blog/cohort-auc/labeled-targets.csv")
    if not labels_file.exists() or not (baseline / "analysis.json").exists():
        pytest.skip("Reproduced blog data is required")
    if not Path("data/raw/hu2021_datasets.xlsx").exists():
        pytest.skip("Downloaded source data is required")
    patients, sequences, _, labels, pairs, inputs, components, _ = AUDIT.candidate_context(
        baseline, labels_file
    )
    assert len(patients) == 6
    assert len(sequences) == 80175
    assert len(pairs) == 272705
    assert sum(p["candidate_mutations"] for p in patients.values()) == 4335
    assert components == 98 and len(labels) == 72
    assert sum(int(row["label"]) for row in labels) == 45
    assert set(inputs) == {
        str(labels_file),
        *(
            str(baseline / name)
            for name in ["analysis.json", "target-rankings.csv", "vaccine-component-audit.csv"]
        ),
    }
    assert all(allele.startswith(("HLA-A", "HLA-B")) for _, allele in pairs)
