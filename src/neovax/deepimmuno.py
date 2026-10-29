import csv
import json
import math
import subprocess
import tempfile
import time
from importlib.resources import files
from pathlib import Path

import httpx

from neovax.iedb import Cancelled, InferenceError
from neovax.io import digest, stable_id, write_json


def setup_model(directory, python):
    directory = Path(directory).resolve()
    specification = json.loads(files("neovax").joinpath("resources/deepimmuno_sources.json").read_text())
    for record in specification["files"]:
        path = directory / record["file"]
        if not path.exists() or digest(path) != record["sha256"]:
            response = httpx.get(record["url"], follow_redirects=True, timeout=120)
            response.raise_for_status()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(response.content)
        if digest(path) != record["sha256"]:
            raise ValueError(f"DeepImmuno hash mismatch: {path}")
    manifest = {
        **specification,
        "root": str(directory),
        "python": str(Path(python).absolute()),
        "runtime": "TensorFlow 2.15.1 / Python 3.11; isolated from the binding environment",
    }
    write_json(directory / "manifest.json", manifest)
    return directory / "manifest.json"


def allele_representations(root, alleles):
    """Mirror the pinned upstream numerical rescue, retaining the actual patient HLA."""
    with (Path(root) / "data/hla2paratopeTable_aligned.txt").open() as handle:
        supported = [row["HLA"] for row in csv.DictReader(handle, delimiter="\t")]
    inventory = {}
    for allele in supported:
        inventory.setdefault(allele[4], {}).setdefault(allele[6:8], []).append(allele[8:])
    result = {}
    for allele in alleles:
        compact = allele.replace(":", "")
        model = compact
        if compact not in supported:
            groups = inventory[compact[4]]
            group, subtype = compact[6:8], compact[8:]
            if group in groups:
                subtype = min(groups[group], key=lambda value: abs(int(subtype) - int(value)))
            else:
                group = min(groups, key=lambda value: abs(int(group) - int(value)))
                subtype = groups[group][0]
            model = f"HLA-{compact[4]}*{group}{subtype}"
        result[allele] = f"{model[:-2]}:{model[-2:]}"
    return result


def validate_model(manifest_path, alleles, lengths, allow_hla_rescue=False):
    manifest = json.loads(Path(manifest_path).read_text())
    pinned = json.loads(files("neovax").joinpath("resources/deepimmuno_sources.json").read_text())
    if manifest["commit"] != pinned["commit"] or manifest["files"] != pinned["files"]:
        raise ValueError("DeepImmuno manifest differs from the pinned model specification")
    root = Path(manifest["root"])
    for record in manifest["files"]:
        if digest(root / record["file"]) != record["sha256"]:
            raise ValueError(f"DeepImmuno file hash mismatch: {record['file']}")
    with (root / "data/hla2paratopeTable_aligned.txt").open() as handle:
        supported = {row["HLA"] for row in csv.DictReader(handle, delimiter="\t")}
    missing = [a for a in alleles if a.replace(":", "") not in supported]
    if (missing and not allow_hla_rescue) or not set(lengths) <= {9, 10}:
        raise ValueError(
            f"DeepImmuno requires exact allele support and 9/10-mers. Unsupported alleles: {missing}. Use binding-only mode."
        )
    return manifest


def training_overlap(manifest, sequences):
    path = Path(manifest["root"]) / "data/remove_low_negative/remove0123_sample100.csv"
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    training = {r["peptide"] for r in rows}
    matches = sorted(set(sequences) & training)
    return {
        "training_file": str(path.name),
        "sha256": digest(path),
        "exact_sequence_matches": matches,
        "query_unique_sequences": len(set(sequences)),
        "training_unique_sequences": len(training),
        "scope": "Exact sequence audit against one released training CSV only. Related sequences, other inputs, development overlap and NetMHCpan overlap remain unaudited. Exploratory retrospective comparison.",
    }


def cached_prediction_subset(directory, manifest, pairs, representations):
    """Reuse a validated model/allele-compatible batch after candidate filtering."""
    wanted = set(pairs)
    if not wanted:
        return []
    for path in sorted(Path(directory).glob("*.json")):
        saved = json.loads(path.read_text())
        model = saved["model"]
        if model["commit"] != manifest["commit"] or model["files"] != manifest["files"]:
            continue
        mapping = saved.get("allele_representations", {})
        if any(mapping.get(allele, allele) != target for allele, target in representations.items()):
            continue
        records = [r for r in saved["predictions"] if (r["peptide"], r["allele"]) in wanted]
        if {(r["peptide"], r["allele"]) for r in records} == wanted:
            return records
    return None


def predict_immunogenicity(
    epitopes, alleles, config, cancel=lambda: False, progress=lambda message: None, *, allow_hla_rescue=False
):
    manifest = validate_model(config.deepimmuno_manifest, alleles, config.lengths, allow_hla_rescue)
    pairs = sorted({(r["peptide"], a) for r in epitopes for a in alleles})
    root = Path(manifest["root"])
    representations = allele_representations(root, alleles)
    rescued = {a: model for a, model in representations.items() if a != model}
    cache_identity = [manifest["commit"], manifest["files"], pairs]
    if rescued:
        cache_identity.append({"rescue_unknown_hla": rescued})
    cache_path = Path(config.cache_dir) / "deepimmuno" / f"{stable_id(*cache_identity)}.json"
    if cache_path.exists():
        records = json.loads(cache_path.read_text())["predictions"]
    else:
        records = cached_prediction_subset(cache_path.parent, manifest, pairs, representations)
    if records is None:
        if config.cached_only:
            raise InferenceError("DeepImmuno scores are not cached; binding results remain available")
        progress(f"DeepImmuno: {len(pairs)} peptide–HLA pairs in an isolated process")
        with tempfile.TemporaryDirectory(prefix="neovax-deepimmuno-") as temporary:
            temp = Path(temporary)
            input_path = temp / "input.csv"
            with input_path.open("w", newline="") as handle:
                csv.writer(handle).writerows((peptide, allele.replace(":", "")) for peptide, allele in pairs)
            command = [
                manifest["python"],
                str(root / "deepimmuno-cnn.py"),
                "--mode",
                "multiple",
                "--intdir",
                str(input_path),
                "--outdir",
                str(temp),
            ]
            with (temp / "inference.log").open("w+") as log:
                process = subprocess.Popen(command, cwd=root, stdout=log, stderr=log)
                started = time.monotonic()
                while process.poll() is None:
                    if cancel() or time.monotonic() - started > 1800:
                        process.terminate()
                        process.wait(timeout=30)
                        raise Cancelled("DeepImmuno inference cancelled or exceeded 30 minutes")
                    time.sleep(0.1)
                if process.returncode:
                    log.seek(0)
                    raise InferenceError(f"DeepImmuno process failed: {log.read()[-1500:]}")
            with (temp / "deepimmuno-cnn-result.txt").open() as handle:
                records = [
                    {
                        "peptide": r["peptide"],
                        "allele": f"{r['HLA'][:-2]}:{r['HLA'][-2:]}",
                        "score": float(r["immunogenicity"]),
                    }
                    for r in csv.DictReader(handle, delimiter="\t")
                ]
    result = {(r["peptide"], r["allele"]): r["score"] for r in records}
    if (
        len(records) != len(pairs)
        or set(result) != set(pairs)
        or any(not math.isfinite(v) or not 0 <= v <= 1 for v in result.values())
    ):
        raise InferenceError("DeepImmuno returned incomplete, duplicate or invalid scores")
    if not cache_path.exists():
        write_json(
            cache_path, {"model": manifest, "allele_representations": representations, "predictions": records}
        )
    return result, {
        "model": manifest,
        "allele_representations": representations,
        "rescued_alleles": rescued,
        "training_overlap": training_overlap(manifest, [r["peptide"] for r in epitopes]),
    }
