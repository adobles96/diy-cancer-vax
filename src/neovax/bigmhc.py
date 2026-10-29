"""Strict subprocess adapter for the optional, pinned BigMHC IM pilot."""

import csv
import json
import math
import os
import subprocess
import time
from pathlib import Path

from neovax.io import digest, stable_id, write_json

COMMIT = "199e890fa69a63646730bcecea0a520e08ca4f02"
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def exact_alleles(root, alleles):
    """Require literal patient HLA support; upstream's fuzzy fallback is forbidden."""
    with (Path(root) / "data/pseudoseqs.csv").open() as handle:
        supported = {row[0] for row in list(csv.reader(handle))[1:]}
    missing = sorted(set(alleles) - supported)
    if missing:
        raise ValueError(f"No exact BigMHC pseudosequence for: {missing}")
    return {allele: allele for allele in sorted(set(alleles))}


def read_predictions(path, pairs):
    expected = set(pairs)
    result = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle):
            pair = (row["pep"], row["mhc"])
            value = float(row["BigMHC_IM"])
            if pair in result or pair not in expected or not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"Duplicate, unexpected, or invalid BigMHC prediction: {pair}")
            result[pair] = value
    if set(result) != expected:
        raise ValueError(f"BigMHC coverage mismatch: {len(result)} of {len(expected)} pairs")
    return result


def predict_immunogenicity(pairs, root, python, output, batch_size=256, threads=8, device="cpu"):
    """Run official -m im CLI and retain original input, output, log, and provenance."""
    root, output = Path(root).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    pairs = sorted(set(pairs))
    if not pairs or any(len(p) not in {9, 10} or not set(p) <= AMINO_ACIDS for p, _ in pairs):
        raise ValueError("Pilot expects nonempty canonical mutation-containing 9/10-mer pairs")
    manifest = json.loads((root / "manifest.json").read_text())
    if manifest["commit"] != COMMIT or manifest["model"] != "BigMHC_IM":
        raise ValueError("BigMHC manifest differs from pinned IM model")
    for record in manifest["files"]:
        if digest(root / record["file"]) != record["sha256"]:
            raise ValueError(f"BigMHC file checksum mismatch: {record['file']}")
    representations = exact_alleles(root, [a for _, a in pairs])
    if device not in {"cpu", "mps"}:
        raise ValueError(f"Unsupported pilot device: {device}")
    runner = (
        root / "src/predict.py"
        if device == "cpu"
        else Path(__file__).resolve().parents[2] / "scripts/run_bigmhc_mps.py"
    )
    runtime = json.loads(
        subprocess.check_output(
            [
                str(Path(python).absolute()),
                "-c",
                "import importlib.metadata, json, platform, torch; "
                "print(json.dumps(dict(python=platform.python_version(), platform=platform.platform(), "
                "torch=torch.__version__, mps_available=torch.backends.mps.is_available(), "
                "packages=dict(sorted((d.metadata['Name'], d.version) for d in importlib.metadata.distributions())))))",
            ],
            text=True,
        )
    )
    runner_hash = digest(runner)
    identity = stable_id(manifest, pairs, runtime, runner_hash, device, batch_size, threads)
    input_path, result_path = output / "input.csv", output / "predictions.csv"
    metadata_path = output / "inference.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        if metadata["identity"] == identity and metadata["output_sha256"] == digest(result_path):
            return read_predictions(result_path, pairs), metadata
    with input_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mhc", "pep"])
        writer.writerows((allele, peptide) for peptide, allele in pairs)
    command = [
        str(Path(python).absolute()),
        str(runner),
        "-i",
        str(input_path),
        "-o",
        str(result_path),
        "-m",
        "im",
        "-d",
        "cpu",
        "-a",
        "0",
        "-p",
        "1",
        "-c",
        "1",
        "-b",
        str(batch_size),
        "-j",
        "1",
        "-f",
        "2",
    ]
    environment = {
        **os.environ,
        "OMP_NUM_THREADS": str(threads),
        "VECLIB_MAXIMUM_THREADS": str(threads),
        "MKL_NUM_THREADS": str(threads),
        "PYTHONUNBUFFERED": "1",
        "PYTHONPATH": str(root / "src"),
        "PYTORCH_ENABLE_MPS_FALLBACK": "1",
    }
    started = time.monotonic()
    with (output / "inference.log").open("w") as log:
        subprocess.run(command, cwd=root, env=environment, stdout=log, stderr=log, check=True)
    predictions = read_predictions(result_path, pairs)
    metadata = {
        "identity": identity,
        "model": manifest,
        "command": command,
        "device": device,
        "threads": threads,
        "batch_size": batch_size,
        "elapsed_seconds": time.monotonic() - started,
        "pairs": len(pairs),
        "allele_representations": representations,
        "input_sha256": digest(input_path),
        "output_sha256": digest(result_path),
        "runtime": runtime,
        "runner_sha256": runner_hash,
    }
    write_json(metadata_path, metadata)
    return predictions, metadata
