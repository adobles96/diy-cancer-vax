"""Reproduce all six patient runs, keeping every output outside version control."""

import argparse
from pathlib import Path

from neovax.config import Config
from neovax.io import write_json
from neovax.pipeline import run_pipeline


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--live", action="store_true", help="Allow missing prediction batches to contact IEDB"
    )
    args = parser.parse_args()
    summaries = []
    for patient in range(1, 7):
        result = run_pipeline(Config(patient=patient, cached_only=not args.live), progress=print)
        summaries.append(result["summary"])
        write_json(Path("data/outputs/six-patient-validation.json"), summaries)
    if any(r["status"] != "complete" for r in summaries):
        raise SystemExit("Some runs are incomplete; inspect manifests and rerun with --live to resume")


if __name__ == "__main__":
    main()
