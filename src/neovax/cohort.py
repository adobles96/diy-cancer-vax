"""Retrospective recall of vaccine components using independently ranked mutation targets."""

from collections import defaultdict

from neovax.ranking import fractional_ranks

METHODS = ("netmhcpan", "deepimmuno", "mean_ensemble")
DEFAULT_RECALL_K = 20


def rank_targets(epitopes, alleles, binding, immunogenicity, representations):
    """Best pair per mutation for each model; mean of target-level fractional ranks."""
    groups = defaultdict(set)
    annotations = {}
    for epitope in epitopes:
        mutation = epitope["mutation_id"]
        annotations[mutation] = epitope
        groups[mutation].update((epitope["peptide"], allele) for allele in alleles)
    rows = []
    for mutation, pairs in sorted(groups.items()):
        if any(pair not in binding or pair not in immunogenicity for pair in pairs):
            raise ValueError(f"Incomplete model scores for {mutation}")
        best_binding = min(pairs, key=lambda pair: (binding[pair]["ic50"], pair))
        best_immuno = min(pairs, key=lambda pair: (-immunogenicity[pair], pair))
        rows.append(
            {
                "mutation_id": mutation,
                "gene": annotations[mutation]["gene"],
                "protein_change": annotations[mutation]["protein_change"],
                "unique_pairs": len(pairs),
                "binding_ic50_nm": binding[best_binding]["ic50"],
                "binding_peptide": best_binding[0],
                "binding_allele": best_binding[1],
                "deepimmuno_score": immunogenicity[best_immuno],
                "deepimmuno_peptide": best_immuno[0],
                "deepimmuno_actual_allele": best_immuno[1],
                "deepimmuno_model_allele": representations[best_immuno[1]],
                "deepimmuno_winner_rescued": best_immuno[1] != representations[best_immuno[1]],
            }
        )
    binding_ranks = fractional_ranks([r["binding_ic50_nm"] for r in rows])
    immuno_ranks = fractional_ranks([r["deepimmuno_score"] for r in rows], descending=True)
    for row, binding_rank, immuno_rank in zip(rows, binding_ranks, immuno_ranks):
        row.update(
            netmhcpan_cost=binding_rank,
            deepimmuno_cost=immuno_rank,
            mean_ensemble_cost=(binding_rank + immuno_rank) / 2,
        )
    for method in METHODS:
        for rank, row in enumerate(sorted(rows, key=lambda r: (r[f"{method}_cost"], r["mutation_id"])), 1):
            row[f"{method}_rank"] = rank
    return rows


def audit_components(vaccine, proteins, assays):
    """Allow publication gene suffixes only with exact consequence and full peptide identity."""
    result = []
    for original in vaccine:
        component = dict(original)
        if component["mutation_id"] is None:
            matches = {
                protein["mutation_id"]
                for protein in proteins
                if component["gene"][:-1] == protein["gene"]
                and component["gene"][-1:] in "abcde"
                and component["protein_change"] == protein["protein_change"]
                and component["immunizing_sequence"] in protein["mutant_protein"]
            }
            if len(matches) == 1:
                component.update(
                    mutation_id=next(iter(matches)),
                    reconstructable=True,
                    target_join="Publication a–e suffix removed; exact protein change, unique genomic mutation, and full mutant long-peptide sequence verified",
                )
        positives = [
            observation
            for observation in assays
            if observation["imp_id"] == component["imp_id"]
            and observation["timepoint"] == "week 16"
            and observation["assay_type"] == "peptide"
            and observation["cell_type"] in {"CD4", "CD8"}
            and observation["result"] == "positive"
        ]
        component.update(
            immunogenic=bool(positives),
            positive_cell_types=sorted({o["cell_type"] for o in positives}),
            positive_observation_ids=[o["observation_id"] for o in positives],
            positive_sources=[source for o in positives for source in o["sources"]],
        )
        result.append(component)
    return result


def recall_at_k(targets, components, k=DEFAULT_RECALL_K):
    candidate_ids = {r["mutation_id"] for r in targets}
    results = []
    for method in METHODS:
        ordered = sorted(targets, key=lambda r: r[f"{method}_rank"])
        selected = ordered[:k]
        selected_ids = {r["mutation_id"] for r in selected}
        cutoff = selected[-1][f"{method}_cost"] if selected else None
        tied = [r for r in ordered if r[f"{method}_cost"] == cutoff]
        for endpoint in ("vaccine", "immunogenic"):
            relevant = [r for r in components if endpoint == "vaccine" or r["immunogenic"]]
            available = [r for r in relevant if r["mutation_id"] in candidate_ids]
            hits = [r for r in relevant if r["mutation_id"] in selected_ids]
            results.append(
                {
                    "method": method,
                    "endpoint": endpoint,
                    "k": k,
                    "candidate_targets": len(targets),
                    "selected_targets": len(selected),
                    "hits": len(hits),
                    "denominator": len(relevant),
                    "recall": len(hits) / len(relevant) if relevant else None,
                    "reconstructable_denominator": len(available),
                    "reconstructable_recall": len(hits) / len(available) if available else None,
                    "coverage_ceiling": len(available) / len(relevant) if relevant else None,
                    "hit_imp_ids": [r["imp_id"] for r in hits],
                    "unavailable_imp_ids": [
                        r["imp_id"] for r in relevant if r["mutation_id"] not in candidate_ids
                    ],
                    "selected_with_rescued_immuno_winner": sum(
                        r["deepimmuno_winner_rescued"] for r in selected
                    ),
                    "boundary_tie_size": len(tied),
                    "boundary_tie_selected": sum(r["mutation_id"] in selected_ids for r in tied),
                }
            )
    return results
