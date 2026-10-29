"""Run: .venv/bin/python scripts/cohort_auc.py (uses the saved cohort scores)."""

import argparse
import json
from pathlib import Path

from neovax.blog import run_cohort_auc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="data/outputs/cohort-recall")
    parser.add_argument("--output-dir", default="data/outputs/cohort-auc")
    parser.add_argument("--data-dir", help="Override the source-data directory in the input provenance")
    args = parser.parse_args()
    result = run_cohort_auc(
        Path(args.input_dir),
        Path(args.output_dir),
        data_dir=args.data_dir,
        progress=lambda message: print(message, flush=True),
    )
    print(json.dumps(result["aggregate"], indent=2))
    print(f"Report: {Path(args.output_dir).resolve() / 'report.html'}")


if __name__ == "__main__":
    main()
