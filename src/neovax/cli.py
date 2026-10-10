import argparse
import signal
import threading
from dataclasses import replace
from pathlib import Path

from rich.console import Console
from rich.table import Table

from neovax.config import Config
from neovax.data import fetch_data
from neovax.deepimmuno import setup_model
from neovax.iedb import BindingCache
from neovax.pipeline import run_pipeline
from neovax.report import render_report
from neovax.tui import NeoVaxApp


def print_summary(result, console=None):
    console = console or Console()
    summary = result["summary"]
    console.print(f"[bold]Patient {summary['patient']} · {summary['status'].upper()}[/bold]")
    table = Table("Gene", "Proposed long peptide", "Driving epitope", "Score", "Historical target")
    for row in result["proposed"]:
        table.add_row(
            row["gene"],
            row["sequence"],
            row["driving_epitope"],
            f"{row['ranking_cost']:.4g}",
            "yes" if row["historical_target"] else "no",
        )
    console.print(table)
    console.print(
        f"Target overlap: {summary['overlap']['target_overlap']} / {summary['overlap']['reconstructable_historical_mutations']} reconstructable historical mutations; {summary['overlap']['historical_components']} total administered components."
    )
    console.print(f"Reconstruction gaps: {summary['coverage']['skipped_reasons']}")
    console.print(f"Assay observations: {summary['assay_states']}")
    for error in summary["errors"]:
        console.print(error, style="yellow", markup=False)
    root = Path(summary["run_directory"])
    for name in [
        "report.html",
        "proposed_vaccine.csv",
        "proposed_vaccine.fasta",
        "assays.csv",
        "manifest.json",
    ]:
        console.print(str(root / name), markup=False)


def main():
    parser = argparse.ArgumentParser(description="NeoVax Explorer: published candidate and assay exploration")
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ["run", "tui"]:
        command = commands.add_parser(name)
        command.add_argument("--config", type=Path)
        command.add_argument("--patient", type=int, choices=range(1, 7))
        command.add_argument("--cached-only", action="store_true", default=None)
    fetch = commands.add_parser("fetch-data")
    fetch.add_argument("--data-dir", default="outputs/data/raw")
    fetch.add_argument("--from-directory", type=Path)
    report = commands.add_parser("report")
    report.add_argument("run_directory", type=Path)
    legacy = commands.add_parser("import-cache")
    legacy.add_argument("directory", type=Path)
    legacy.add_argument("--cache-dir", default="outputs/cache")
    deep = commands.add_parser("setup-deepimmuno")
    deep.add_argument("--directory", default="outputs/models/deepimmuno")
    deep.add_argument(
        "--python", required=True, help="Python executable in a separate TensorFlow 2.15.1 environment"
    )
    args = parser.parse_args()
    console = Console()
    if args.command == "fetch-data":
        fetch_data(args.data_dir, args.from_directory)
        console.print(f"Pinned sources verified in {Path(args.data_dir).resolve()}")
    elif args.command == "report":
        console.print(str(render_report(args.run_directory)))
    elif args.command == "import-cache":
        cache = BindingCache(args.cache_dir)
        try:
            console.print(f"Imported {cache.import_legacy(args.directory)} validated prediction rows")
        finally:
            cache.close()
    elif args.command == "setup-deepimmuno":
        console.print(str(setup_model(args.directory, args.python)))
    else:
        config = Config.load(args.config) if args.config else Config()
        overrides = {
            key: getattr(args, key) for key in ["patient", "cached_only"] if getattr(args, key) is not None
        }
        config = replace(config, **overrides)
        if args.command == "tui":
            app = NeoVaxApp(config)
            app.run()
            if app.result:
                print_summary(app.result, console)
        else:
            event = threading.Event()
            previous = signal.signal(signal.SIGINT, lambda signum, frame: event.set())
            try:
                result = run_pipeline(
                    config, progress=lambda message: console.print(message, markup=False), cancel=event.is_set
                )
            finally:
                signal.signal(signal.SIGINT, previous)
            print_summary(result, console)
            if result["summary"]["status"] != "complete":
                raise SystemExit(2)


if __name__ == "__main__":
    main()
