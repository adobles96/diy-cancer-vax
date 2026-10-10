import datetime
import json
import uuid
from collections import Counter
from dataclasses import replace
from importlib.resources import files
from pathlib import Path

from neovax import __version__
from neovax.data import actual_vaccine, load_reference, patient_data, verify_sources
from neovax.deepimmuno import predict_immunogenicity
from neovax.evaluation import annotate_membership, assay_matrix, evaluate, rescore_historical
from neovax.iedb import METHOD, BindingCache, Cancelled, InferenceError, contexts_for, predict
from neovax.io import digest, write_csv, write_json
from neovax.peptides import coverage, enumerate_epitopes, reconstruct
from neovax.ranking import rank_epitopes, summarize_mutations
from neovax.report import render_report
from neovax.selection import select_vaccine


def run_pipeline(config, progress=lambda message: None, cancel=lambda: False, make_report=True):
    now = datetime.datetime.now(datetime.UTC)
    run_dir = (
        Path(config.output_dir).resolve()
        / f"patient-{config.patient}-{now.strftime('%Y%m%dT%H%M%S')}-{uuid.uuid4().hex[:8]}"
    )
    run_dir.mkdir(parents=True)
    write_json(run_dir / "config.json", config.to_dict())
    manifest = {
        "status": "running",
        "created_at": now.isoformat(),
        "software_version": __version__,
        "binding_model": METHOD,
        "historical_model": "NetMHCpan 2.4",
        "coordinate_system": "protein: zero-based half-open; genomic: published GRCh37 coordinates",
        "run_directory": str(run_dir),
    }
    write_json(run_dir / "manifest.json", manifest)
    try:
        result = execute(config, run_dir, manifest, progress, cancel)
    except (OSError, ValueError, KeyError) as error:
        manifest.update(status="failed", errors=[str(error)])
        write_json(run_dir / "manifest.json", manifest)
        raise
    if make_report:
        render_report(run_dir, result)
    manifest["exports"] = {
        p.relative_to(run_dir).as_posix(): digest(p)
        for p in run_dir.rglob("*")
        if p.is_file() and p.name != "manifest.json"
    }
    write_json(run_dir / "manifest.json", manifest)
    progress(f"{manifest['status'].capitalize()} · {run_dir}")
    return result


def execute(config, run_dir, manifest, progress, cancel):
    progress(f"Patient {config.patient}: normalizing published variants and observations")
    variants, alleles, published, assays = patient_data(config.data_dir, config.patient)
    manifest["sources"] = verify_sources(Path(config.data_dir))
    proteins, skipped = reconstruct(variants, load_reference(config.data_dir))
    epitopes = enumerate_epitopes(proteins, config.lengths)
    counts = coverage(variants, proteins, skipped, epitopes)
    progress(f"Reconstructed {len(proteins)} missense rows · {len(epitopes)} epitope associations")
    write_csv(run_dir / "variants.csv", variants)
    write_csv(run_dir / "skipped_variants.csv", skipped)
    write_csv(run_dir / "reconstructed_proteins.csv", proteins)
    cache = BindingCache(config.cache_dir)
    try:
        predictions, workload, status, errors = predict(
            cache,
            contexts_for(proteins, config.lengths),
            alleles,
            config.lengths,
            config.batch_size,
            config.cached_only,
            cancel,
            progress,
        )
    finally:
        cache.close()
    expected = {
        (e[field], a) for e in epitopes for field in ["peptide", "reference_peptide"] for a in alleles
    }
    missing = expected - predictions.keys()
    if missing and status == "complete":
        status = "partial"
    if cancel():
        status = "cancelled"
    immuno, immuno_metadata, ranking_available = {}, None, True
    if config.ranking != "binding" and status != "cancelled":
        try:
            immuno, immuno_metadata = predict_immunogenicity(epitopes, alleles, config, cancel, progress)
        except (InferenceError, ValueError, OSError) as error:
            status = "cancelled" if isinstance(error, Cancelled) else "partial"
            errors.append(str(error))
            ranking_available = False
    elif config.ranking != "binding":
        ranking_available = False
    ranking_config = config if ranking_available else replace(config, ranking="binding")
    pairs, ranked = rank_epitopes(epitopes, alleles, predictions, ranking_config, immuno)
    for row in ranked:
        row["ranking_available"] = ranking_available
    if not ranking_available:
        for row in [*pairs, *ranked]:
            row.update(ranking_cost=None, eligible=False)
            if "epitope_rank" in row:
                row["epitope_rank"] = None
    mutation_scores = summarize_mutations(ranked, config.aggregation)
    selected, selection_summary = select_vaccine(ranked, mutation_scores, proteins, config)
    # Historical labels enter only after every discovery prediction, rank and selection is fixed.
    vaccine = actual_vaccine(published, variants, proteins)
    if status != "cancelled":
        progress("Comparison: rescoring published long sequences separately from discovery")
        historical_sequences = {r["immunizing_sequence"] for r in vaccine}
        historical_sequences.update(r["wildtype_epitope"] for r in published if r["wildtype_epitope"])
        historical_sequences.update(r["mutant_epitope"] for r in published if r["mutant_epitope"])
        cache = BindingCache(config.cache_dir)
        try:
            historical_predictions, historical_workload, historical_status, historical_errors = predict(
                cache,
                sorted(historical_sequences),
                alleles,
                config.lengths,
                config.batch_size,
                config.cached_only,
                cancel,
                progress,
            )
        finally:
            cache.close()
        predictions.update(historical_predictions)
        manifest["historical_inference"] = {
            "status": historical_status,
            "workload": historical_workload,
            "errors": historical_errors,
        }
        if historical_status != "complete" and status == "complete":
            status = historical_status
        errors.extend(historical_errors)
    ranked = annotate_membership(ranked, vaccine, published, assays)
    selected = annotate_membership(selected, vaccine, published, assays)
    protein_map = {r["variant_id"]: r for r in proteins}
    for component in selected:
        protein = protein_map[component["variant_id"]]
        contexts = []
        for actual in vaccine:
            if actual["mutation_id"] != component["mutation_id"]:
                continue
            start = protein["mutant_protein"].find(actual["immunizing_sequence"])
            contexts.append(
                {
                    "imp_id": actual["imp_id"],
                    "start": start if start >= 0 else None,
                    "end": start + len(actual["immunizing_sequence"]) if start >= 0 else None,
                    "sequence": actual["immunizing_sequence"],
                }
            )
        component["historical_boundaries"] = contexts
    evaluations = evaluate(selected, ranked, vaccine, published, assays, predictions, config.budget)
    historical_rescore = rescore_historical(vaccine, predictions, alleles, config.lengths)
    supporting = json.loads(files("neovax").joinpath("resources/observations.json").read_text())
    supporting = [r for r in supporting if config.patient in r["patients"]]
    counts.update(
        expected_unique_predictions=len(expected),
        completed_unique_predictions=len(expected) - len(missing),
        missing_unique_predictions=len(missing),
        eligible_unique_pairs=sum(p["eligible"] for p in pairs),
        proposed_long_peptides=len(selected),
    )
    manifest.update(
        status=status,
        errors=errors,
        workload=workload,
        immunogenicity=immuno_metadata,
        requested_ranking_available=ranking_available,
        request_ids=sorted({predictions[k]["request_id"] for k in expected if k in predictions}),
    )
    summary = {
        "patient": config.patient,
        "status": status,
        "alleles": alleles,
        "coverage": counts,
        "selection": selection_summary,
        "overlap": evaluations["overlap"],
        "assay_states": dict(Counter(r["result"] for r in assays)),
        "errors": errors,
        "interpretation": evaluations["interpretation"],
        "run_directory": str(run_dir),
    }
    result = {
        "config": config.to_dict(),
        "manifest": manifest,
        "summary": summary,
        "proposed": selected,
        "epitopes": ranked,
        "pairs": pairs,
        "mutations": mutation_scores,
        "actual": vaccine,
        "assays": assays,
        "assay_matrix": assay_matrix(assays),
        "supporting": supporting,
        "evaluations": evaluations,
        "historical_rescore": historical_rescore,
        "skipped": skipped,
    }
    for name, rows in [
        ("epitopes", ranked),
        ("pair_predictions", pairs),
        ("mutation_scores", mutation_scores),
        ("proposed_vaccine", selected),
        ("actual_vaccine", vaccine),
        ("assays", assays),
        ("historical_rescore", historical_rescore),
    ]:
        write_csv(run_dir / f"{name}.csv", rows)
    # Aligned reference scores for every association and allele, not just the winning allele.
    write_csv(
        run_dir / "reference_predictions.csv",
        [
            {
                "epitope_id": e["epitope_id"],
                "reference_peptide": e["reference_peptide"],
                "allele": a,
                "ic50_nm": predictions[(e["reference_peptide"], a)]["ic50"]
                if (e["reference_peptide"], a) in predictions
                else None,
            }
            for e in epitopes
            for a in alleles
        ],
    )
    (run_dir / "proposed_vaccine.fasta").write_text(
        "".join(
            f">{r['component_id']} {r['gene']} {r['protein_change']}\n{r['sequence']}\n" for r in selected
        )
    )
    for name, value in [
        ("manifest", manifest),
        ("summary", summary),
        ("evaluations", evaluations),
        ("run_data", result),
    ]:
        write_json(run_dir / f"{name}.json", value)
    return result
