"""Run: .venv/bin/python scripts/cohort_recall.py --rescue-unknown-hla"""

import argparse
import json
from pathlib import Path

from neovax.blog import run_cohort_recall
from neovax.config import Config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rescue-unknown-hla", action="store_true")
    parser.add_argument("--cached-only", action="store_true")
    parser.add_argument("--data-dir", default="data/raw")
    parser.add_argument("--cache-dir", default="data/outputs/cache")
    parser.add_argument("--deepimmuno-manifest", default="data/outputs/models/deepimmuno/manifest.json")
    parser.add_argument("--output-dir", default="data/outputs/cohort-recall")
    args = parser.parse_args()
    config = Config(
        data_dir=args.data_dir,
        cache_dir=args.cache_dir,
        cached_only=args.cached_only,
        deepimmuno_manifest=args.deepimmuno_manifest,
    )
    result = run_cohort_recall(
        config,
        Path(args.output_dir),
        progress=lambda message: print(message, flush=True),
        allow_hla_rescue=args.rescue_unknown_hla,
    )
    print(json.dumps(result["pooled"], indent=2))
    print(f"Report: {Path(args.output_dir).resolve() / 'report.html'}")


if __name__ == "__main__":
    main()
