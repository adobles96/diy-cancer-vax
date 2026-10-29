from pathlib import Path

import pytest

from neovax.cohort import audit_components, rank_targets, recall_at_k
from neovax.data import actual_vaccine, load_reference, patient_data
from neovax.deepimmuno import allele_representations, validate_model
from neovax.peptides import enumerate_epitopes, reconstruct


def example_ranking(scale=1):
    epitopes = [
        {"mutation_id": f"m{i}", "peptide": f"pep{i}", "gene": "gene", "protein_change": "p.A1V"}
        for i in range(3)
    ]
    binding = {(f"pep{i}", "HLA-A*01:01"): {"ic50": value * scale} for i, value in enumerate([1, 10, 1000])}
    immuno = {(f"pep{i}", "HLA-A*01:01"): value for i, value in enumerate([0.1, 0.9, 0.5])}
    return epitopes, binding, immuno


def test_target_rank_ensemble_and_no_binding_gate():
    epitopes, binding, immuno = example_ranking()
    mapping = {"HLA-A*01:01": "HLA-A*01:01"}
    result = rank_targets(epitopes + epitopes, list(mapping), binding, immuno, mapping)
    assert len(result) == 3
    assert all(r["unique_pairs"] == 1 for r in result)
    assert [r["mean_ensemble_cost"] for r in result] == [0.5, 0.25, 0.75]
    assert result[1]["mean_ensemble_rank"] == 1
    assert result[2]["deepimmuno_rank"] == 2  # IC50 > 500 remains in this standalone ranking.
    _, rescaled, _ = example_ranking(100)
    again = rank_targets(epitopes, list(mapping), rescaled, immuno, mapping)
    assert [r["mean_ensemble_rank"] for r in result] == [r["mean_ensemble_rank"] for r in again]
    del immuno[("pep2", "HLA-A*01:01")]
    with pytest.raises(ValueError, match="Incomplete"):
        rank_targets(epitopes, list(mapping), binding, immuno, mapping)


def test_full_and_conditional_recall_keep_unavailable_and_unknown_separate():
    epitopes, binding, immuno = example_ranking()
    mapping = {"HLA-A*01:01": "HLA-A*01:01"}
    targets = rank_targets(epitopes, list(mapping), binding, immuno, mapping)
    components = [
        {"imp_id": "a", "mutation_id": "m1", "immunogenic": True},
        {"imp_id": "b", "mutation_id": "m2", "immunogenic": False},
        {"imp_id": "c", "mutation_id": None, "immunogenic": True},
    ]
    rows = recall_at_k(targets, components, 1)
    vaccine = next(r for r in rows if r["method"] == "deepimmuno" and r["endpoint"] == "vaccine")
    positive = next(r for r in rows if r["method"] == "deepimmuno" and r["endpoint"] == "immunogenic")
    assert vaccine["recall"] == 1 / 3
    assert vaccine["reconstructable_recall"] == 0.5
    assert positive["recall"] == 0.5
    assert positive["reconstructable_recall"] == 1
    assert positive["unavailable_imp_ids"] == ["c"]
    assert positive["selected_targets"] == 1


def test_broad_immunogenicity_and_verified_gene_suffix():
    vaccine = [
        {
            "imp_id": "v1",
            "gene": "GENEa",
            "protein_change": "p.A1V",
            "mutation_id": None,
            "immunizing_sequence": "VPEPTIDE",
            "reconstructable": False,
        }
    ]
    proteins = [
        {"gene": "GENE", "protein_change": "p.A1V", "mutant_protein": "VPEPTIDEXYZ", "mutation_id": "m1"}
    ]
    base = {
        "imp_id": "v1",
        "timepoint": "week 16",
        "assay_type": "peptide",
        "cell_type": "CD4",
        "result": "positive",
        "observation_id": "o1",
        "sources": [{"cell": "K5"}],
    }
    audited = audit_components(vaccine, proteins, [base])[0]
    assert audited["immunogenic"] and audited["reconstructable"] and audited["mutation_id"] == "m1"
    assert vaccine[0]["mutation_id"] is None
    for change in [{"timepoint": "year 3–4.5"}, {"assay_type": "tumor"}, {"result": "not determined"}]:
        assert not audit_components(vaccine, proteins, [{**base, **change}])[0]["immunogenic"]
    mismatched = [{**proteins[0], "mutant_protein": "DIFFERENT"}]
    assert audit_components(vaccine, mismatched, [base])[0]["mutation_id"] is None
    ambiguous = [*proteins, {**proteins[0], "mutation_id": "m2"}]
    assert audit_components(vaccine, ambiguous, [base])[0]["mutation_id"] is None


def test_numerical_rescue_keeps_upstream_table_order(tmp_path):
    table = tmp_path / "data/hla2paratopeTable_aligned.txt"
    table.parent.mkdir()
    table.write_text("HLA\tpseudo\nHLA-B*5706\tAAA\nHLA-B*5701\tBBB\nHLA-B*2703\tCCC\nHLA-B*2705\tDDD\n")
    result = allele_representations(tmp_path, ["HLA-B*56:01", "HLA-B*27:02", "HLA-B*57:01"])
    assert result == {
        "HLA-B*56:01": "HLA-B*57:06",
        "HLA-B*27:02": "HLA-B*27:03",
        "HLA-B*57:01": "HLA-B*57:01",
    }


@pytest.mark.integration
def test_real_cohort_denominators_and_rescue():
    if not Path("data/raw/ott2017_table5.xlsx").exists():
        pytest.skip("Downloaded sources unavailable")
    counts, available, positive = [], [], []
    for patient in range(1, 7):
        variants, alleles, published, assays = patient_data("data/raw", patient)
        proteins, _ = reconstruct(variants, load_reference("data/raw"))
        components = audit_components(actual_vaccine(published, variants, proteins), proteins, assays)
        candidates = {r["mutation_id"] for r in enumerate_epitopes(proteins, (9, 10))}
        counts.append(len(components))
        available.append(sum(c["mutation_id"] in candidates for c in components))
        positive.append(sum(c["immunogenic"] for c in components))
    assert counts == [13, 17, 14, 14, 20, 20]
    assert available == [11, 10, 12, 9, 15, 15]
    assert positive == [9, 11, 9, 9, 11, 14]
    manifest = "data/outputs/models/deepimmuno/manifest.json"
    if Path(manifest).exists():
        with pytest.raises(ValueError, match="exact allele"):
            validate_model(manifest, alleles, (9, 10))
        validate_model(manifest, alleles, (9, 10), allow_hla_rescue=True)
        with pytest.raises(ValueError):
            validate_model(manifest, alleles, (8,), allow_hla_rescue=True)
