"""Exclude mutant short peptides already present in the pinned reference proteome."""

from collections import defaultdict

import ahocorasick

from neovax.data import source_manifest


def filter_reference_peptides(epitopes, reference):
    """Exact whole-peptide matches, scanning each reference protein separately."""
    unique = {r["peptide"] for r in epitopes}
    matches = defaultdict(list)
    if unique:
        automaton = ahocorasick.Automaton()
        for peptide in sorted(unique):
            automaton.add_word(peptide, peptide)
        automaton.make_automaton()
        for transcript, protein in sorted(reference.items()):
            for last, peptide in automaton.iter(protein):
                matches[peptide].append(
                    {
                        "transcript": transcript,
                        "start": last + 1 - len(peptide),
                        "end": last + 1,
                    }
                )
    retained, removed = [], []
    for row in epitopes:
        if row["peptide"] in matches:
            removed.append(
                {
                    **row,
                    "reason": "Exact peptide match in reference proteome",
                    "reference_matches": matches[row["peptide"]],
                }
            )
        else:
            retained.append(row)
    before_targets = {r["mutation_id"] for r in epitopes}
    after_targets = {r["mutation_id"] for r in retained}
    summary = {
        "mode": "exact whole-peptide amino-acid sequence match",
        "reference_protein_records": len(reference),
        "reference_residues": sum(map(len, reference.values())),
        "candidate_associations_before": len(epitopes),
        "candidate_associations_removed": len(removed),
        "candidate_associations_after": len(retained),
        "unique_peptides_before": len(unique),
        "unique_peptides_removed": len(matches),
        "unique_peptides_after": len(unique) - len(matches),
        "mutations_affected": len({r["mutation_id"] for r in removed}),
        "mutations_fully_removed": sorted(before_targets - after_targets),
        "coordinates": "reference protein: zero-based half-open",
        "scope": "All pinned reference protein records, including isoforms and the source protein. No approximate matches, patient germline variants, expression or processing assumptions. Absence is not proof of tumor specificity.",
    }
    return retained, removed, summary


def reference_filter_source():
    return next(r for r in source_manifest() if r["file"] == "hg19_knownGeneTxPep.txt.gz")
