import re
from collections import Counter

from neovax.io import stable_id

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def reconstruct(variants, proteins):
    reconstructed, skipped = [], []
    for variant in variants:
        row = dict(variant)
        protein = proteins.get(row["transcript"])
        match = re.fullmatch(r"p\.([A-Z])(\d+)([A-Z])", row["protein_change"] or "")
        reason = None
        if row["consequence"] != "Missense_Mutation":
            reason = "unsupported consequence"
        elif protein is None:
            reason = "exact transcript version missing"
        elif not match:
            reason = "unsupported protein change"
        else:
            reference, position, mutant = match[1], int(match[2]) - 1, match[3]
            if position < 0 or position >= len(protein):
                reason = "protein position out of bounds"
            elif protein[position] != reference:
                reason = "reference amino acid mismatch"
            elif mutant not in AMINO_ACIDS or reference not in AMINO_ACIDS:
                reason = "ambiguous residues"
            elif reference == mutant:
                reason = "unchanged residue"
            else:
                row.update(
                    position=position,
                    ambiguous_protein_positions=[i for i, aa in enumerate(protein) if aa not in AMINO_ACIDS],
                    reference_protein=protein,
                    mutant_protein=f"{protein[:position]}{mutant}{protein[position + 1 :]}",
                )
        if reason:
            skipped.append({**variant, "reason": reason})
        else:
            reconstructed.append(row)
    return reconstructed, skipped


def enumerate_epitopes(proteins, lengths):
    rows = []
    for protein in proteins:
        p, n = protein["position"], len(protein["mutant_protein"])
        for length in lengths:
            for start in range(max(0, p - length + 1), min(p, n - length) + 1):
                end = start + length
                if set(protein["mutant_protein"][start:end]) - AMINO_ACIDS:
                    continue
                rows.append(
                    {
                        "epitope_id": stable_id(protein["variant_id"], start, end),
                        "variant_id": protein["variant_id"],
                        "mutation_id": protein["mutation_id"],
                        "patient": protein["patient"],
                        "gene": protein["gene"],
                        "protein_change": protein["protein_change"],
                        "transcript": protein["transcript"],
                        "start": start,
                        "end": end,
                        "mutation_offset": p - start,
                        "length": length,
                        "peptide": protein["mutant_protein"][start:end],
                        "reference_peptide": protein["reference_protein"][start:end],
                    }
                )
    return rows


def coverage(variants, proteins, skipped, epitopes):
    return {
        "variant_rows": len(variants),
        "missense_rows": sum(v["consequence"] == "Missense_Mutation" for v in variants),
        "reconstructed_rows": len(proteins),
        "genomic_mutations": len({r["mutation_id"] for r in proteins}),
        "epitope_associations": len(epitopes),
        "unique_mutant_epitopes": len({r["peptide"] for r in epitopes}),
        "skipped_reasons": dict(Counter(r["reason"] for r in skipped)),
    }
