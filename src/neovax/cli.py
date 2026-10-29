"""A small guided command line, with explicit commands for unattended runs."""

import argparse
import json
import signal
import subprocess
import sys
import threading
from dataclasses import replace
from pathlib import Path

import httpx
from rich import box
from rich.console import Console
from rich.prompt import Confirm, IntPrompt, Prompt
from rich.table import Table

from neovax.blog import LABELS, reproduce_blog
from neovax.config import Config
from neovax.data import fetch_data, source_manifest, verify_sources
from neovax.deepimmuno import setup_model
from neovax.iedb import BindingCache, Cancelled, InferenceError
from neovax.io import write_json
from neovax.pipeline import run_pipeline
from neovax.report import render_report
from neovax.runtime import ensure_deepimmuno

METHODS = {"binding": "NetMHCpan", "immunogenicity": "DeepImmuno", "combined": "Ensemble"}


def prepare_inputs(config, data_root, explicit_config, console):
    settings_path = Path(".neovax.json")
    remember = data_root is not None
    if data_root is None and not explicit_config and settings_path.exists():
        data_root = Path(json.loads(settings_path.read_text())["data_root"])
    if data_root is None and not explicit_config:
        data_root = Path("data")
        missing = any(not (data_root / "raw" / source["file"]).exists() for source in source_manifest())
        if missing:
            console.print("[yellow]Source data not found.[/yellow] The download is approximately 255 MB.")
            if not sys.stdin.isatty():
                raise ValueError(
                    "Choose a data directory with --data-root PATH when running non-interactively."
                )
            answer = console.input(
                "Data directory (contains raw/ and outputs/) [data]: ", markup=False
            ).strip()
            data_root = Path(answer) if answer else data_root
            remember = True
    if data_root is not None:
        data_root = data_root.expanduser().resolve()
        config = replace(
            config,
            data_dir=str(data_root / "raw"),
            cache_dir=str(data_root / "outputs/cache"),
            output_dir=str(data_root / "outputs/runs"),
        )
    console.print(f"Checking source data: {config.data_dir}", style="cyan", markup=False)
    if config.cached_only:
        verify_sources(Path(config.data_dir))
        console.print("Source data verified. Offline mode: using cached predictions only.")
    else:
        fetch_data(
            config.data_dir, progress=lambda message: console.print(message, markup=False, soft_wrap=True)
        )
        console.print("Source data ready.")
    if remember:
        write_json(settings_path, {"data_root": str(data_root)})
    return config


def patient_settings(console):
    patient = int(Prompt.ask("Patient", choices=[str(p) for p in range(1, 7)], default="1", console=console))
    console.print(
        "\n[bold]Default settings[/bold]\n"
        "NetMHCpan · 9/10-amino-acid epitopes · binding <500 nM\n"
        "Up to 20 mutation targets · one 25-amino-acid peptide per target"
    )
    if Confirm.ask("Use defaults?", default=True, console=console):
        return Config(patient=patient), "binding"
    console.print("\n1  NetMHCpan\n2  DeepImmuno\n3  Equal-weight ensemble")
    method = {"1": "binding", "2": "immunogenicity", "3": "combined"}[
        Prompt.ask("Method", choices=["1", "2", "3"], default="1", console=console)
    ]
    lengths = Prompt.ask("Epitope lengths", choices=["9,10", "9", "10"], default="9,10", console=console)
    budget = ask_integer(console, "Maximum mutation targets", 20, 1, 100)
    long_length = ask_integer(console, "Long-peptide length", 25, 15, 30)
    gate = Prompt.ask("Binding cutoff (nM)", choices=["500", "50", "none"], default="500", console=console)
    rescue = False
    if method != "binding":
        console.print(
            "DeepImmuno needs a separate Python 3.11 runtime; model setup is automatic.\n"
            "Some HLA alleles lack exact support. Rescue substitutes a related allele (an approximation)."
        )
        rescue = Confirm.ask("Allow HLA rescue, as in the blog analysis?", default=True, console=console)
    # Attach the actual model manifest after data/model setup, before activating the chosen ranking.
    return Config(
        patient=patient,
        lengths=tuple(int(n) for n in lengths.split(",")),
        budget=budget,
        long_length=long_length,
        threshold=None if gate == "none" else float(gate),
        rescue_unknown_hla=rescue,
    ), method


def ask_integer(console, label, default, minimum, maximum):
    while True:
        value = IntPrompt.ask(f"{label} ({minimum}–{maximum})", default=default, console=console)
        if minimum <= value <= maximum:
            return value
        console.print(f"Enter a number from {minimum} to {maximum}.", style="yellow")


def print_summary(result, console=None):
    console = console or Console()
    summary = result["summary"]
    console.print(f"\n[bold]Patient {summary['patient']} · {summary['status'].upper()}[/bold]")
    binding = result["config"]["ranking"] == "binding"
    score_label = "IC50 (nM)" if binding else "Rank cost"
    table = Table(
        "Gene / variant",
        "Long peptide",
        "Driving epitope",
        score_label,
        "Historical target",
        box=box.SIMPLE,
        show_lines=True,
    )
    for column in table.columns:
        column.overflow = "fold"
    for row in result["proposed"]:
        table.add_row(
            f"{row['gene']}\n{row['protein_change']}",
            row["sequence"],
            row["driving_epitope"],
            f"{row['ranking_cost']:.4g}",
            "yes" if row["historical_target"] else "no",
        )
    console.print(table)
    console.print("Score: lower is better. The driving epitope is the short peptide used for ranking.")
    console.print(
        "Historical target = the same mutation was targeted by NeoVax; the long sequence may differ."
    )
    overlap = summary["overlap"]
    console.print(
        f"Recovered {overlap['target_overlap']} of {overlap['historical_components']} vaccine components "
        f"({overlap['reconstructable_historical_mutations']} reconstructable mutation targets)."
    )
    for error in summary["errors"]:
        console.print(error, style="yellow", markup=False)
    root = Path(summary["run_directory"])
    console.print(f"\nResults: {root}", style="green", markup=False, soft_wrap=True)
    console.print(f"Open: {root / 'report.html'}", markup=False, soft_wrap=True)
    console.print(
        "Shortlist: proposed_vaccine.csv and proposed_vaccine.fasta · All scores and evidence saved alongside."
    )


def print_blog_summary(summary, console):
    console.print("\n[bold green]Blog results reproduced · all six patients[/bold green]")
    table = Table("Method", "Vaccine recall @20", "Peptide recognition¹", "ROC-AUC²", box=box.SIMPLE)
    for row in summary["vaccine_recall"]:
        method = row["method"]
        recognition = next(r for r in summary["peptide_recognition"] if r["method"] == method)
        auc = next(r for r in summary["auc"] if r["method"] == method)
        table.add_row(
            LABELS[method],
            f"{row['hits']}/{row['denominator']} ({row['micro_recall']:.1%})",
            f"{recognition['positive']}/{recognition['included']} ({recognition['fraction']:.1%})",
            f"{auc['within_patient_pair_weighted_auc']:.3f}",
        )
    console.print(table)
    console.print(
        "Top 20 means distinct mutations, each ranked from its retained 9/10-amino-acid epitopes.\n"
        "¹ Any positive week-16 CD4/CD8 peptide-response assay; only selected, tested vaccine targets.\n"
        f"² Within-patient pair-weighted AUC on {summary['labeled_targets']} reconstructable tested targets.\n"
        "Tumor-recognition results use only tested targets and are saved in a separate table."
    )
    console.print(f"\nOpen: {summary['report_path']}", style="green", markup=False, soft_wrap=True)
    console.print(f"Charts, tables and provenance: {summary['output_dir']}", markup=False, soft_wrap=True)


def execute_job(config, command, output_dir, console):
    event = threading.Event()

    def stop(signum, frame):
        if event.is_set():
            raise KeyboardInterrupt
        event.set()
        console.print(
            "Cancellation requested; saving completed work. Press Ctrl+C again to exit immediately.",
            style="yellow",
        )

    previous = signal.signal(signal.SIGINT, stop)
    try:

        def progress(message):
            console.print(message, markup=False, soft_wrap=True)

        if command == "blog":
            summary = reproduce_blog(config, output_dir, progress=progress, cancel=event.is_set)
            print_blog_summary(summary, console)
        else:
            console.print(
                f"\nRunning Patient {config.patient} · {METHODS[config.ranking]}", style="bold cyan"
            )
            result = run_pipeline(config, progress=progress, cancel=event.is_set)
            print_summary(result, console)
            if result["summary"]["status"] != "complete":
                raise SystemExit(2)
    finally:
        signal.signal(signal.SIGINT, previous)


def main():
    parser = argparse.ArgumentParser(
        description="DIY Cancer Vaccine — run a patient or reproduce the blog results"
    )
    commands = parser.add_subparsers(dest="command")
    for name, help_text in [("run", "Run one patient"), ("blog", "Reproduce the six-patient blog results")]:
        command = commands.add_parser(name, help=help_text)
        command.add_argument("--cached-only", action="store_true", default=None)
        command.add_argument(
            "--data-root", type=Path, help="Parent of raw/ and outputs/; remembered for future runs"
        )
        if name == "run":
            command.add_argument("--config", type=Path)
            command.add_argument("--patient", type=int, choices=range(1, 7))
        else:
            command.add_argument("--output-dir", type=Path, help="Default: <data-root>/outputs/blog")
    fetch = commands.add_parser("fetch-data", help="Download or verify source files")
    fetch.add_argument("--data-dir", default="data/raw")
    fetch.add_argument("--from-directory", type=Path)
    report = commands.add_parser("report", help="Regenerate an existing patient HTML report")
    report.add_argument("run_directory", type=Path)
    legacy = commands.add_parser("import-cache")
    legacy.add_argument("directory", type=Path)
    legacy.add_argument("--cache-dir", default="data/outputs/cache")
    deep = commands.add_parser("setup-deepimmuno")
    deep.add_argument("--directory", default="data/outputs/models/deepimmuno")
    deep.add_argument("--python", required=True)
    args = parser.parse_args()
    console = Console()
    try:
        if args.command == "fetch-data":
            fetch_data(
                args.data_dir,
                args.from_directory,
                progress=lambda message: console.print(message, markup=False, soft_wrap=True),
            )
            console.print(f"Pinned sources verified in {Path(args.data_dir).resolve()}")
            return
        if args.command == "report":
            console.print(str(render_report(args.run_directory)))
            return
        if args.command == "import-cache":
            cache = BindingCache(args.cache_dir)
            try:
                console.print(f"Imported {cache.import_legacy(args.directory)} validated prediction rows")
            finally:
                cache.close()
            return
        if args.command == "setup-deepimmuno":
            console.print(str(setup_model(args.directory, args.python)))
            return
        command = args.command
        config_path = getattr(args, "config", None)
        config = Config.load(config_path) if config_path else Config()
        method = config.ranking
        if command is None:
            if not sys.stdin.isatty():
                parser.error(
                    "The menu needs a terminal. Use ./run --patient 1 --data-root PATH or ./run blog --data-root PATH."
                )
            console.print("\n[bold cyan]DIY Cancer Vaccine[/bold cyan]\n")
            console.print("1  Run an individual patient\n2  Reproduce the blog results\n")
            choice = Prompt.ask("Choose", choices=["1", "2"], default="1", console=console)
            command = "run" if choice == "1" else "blog"
            if command == "run":
                config, method = patient_settings(console)
        overrides = {
            key: getattr(args, key)
            for key in ["patient", "cached_only"]
            if getattr(args, key, None) is not None
        }
        config = replace(config, **overrides)
        if command == "blog":
            console.print(
                "\nReproduce all six patients with NetMHCpan, DeepImmuno and their equal-weight ensemble."
            )
            console.print(
                "Creates recall@20 and ROC-AUC charts, peptide/tumor recognition tables, and an HTML report."
            )
            console.print(
                "First run downloads source data and models; uncached predictions can take a while."
            )
        config = prepare_inputs(config, getattr(args, "data_root", None), config_path is not None, console)
        data_root = Path(config.data_dir).parent
        output_dir = getattr(args, "output_dir", None) or data_root / "outputs/blog"
        if command == "blog":
            console.print(
                f"\nBlog results will be saved in: {output_dir.resolve()}", style="cyan", markup=False
            )
            console.print("Uses the blog's fixed settings, including approximate HLA rescue where needed.")
        if command == "blog" or method != "binding":
            manifest = config.deepimmuno_manifest or ensure_deepimmuno(
                data_root,
                progress=lambda message: console.print(message, markup=False, soft_wrap=True),
                cached_only=config.cached_only,
            )
            config = replace(config, deepimmuno_manifest=manifest, ranking=method)
        execute_job(config, command, output_dir, console)
    except (KeyboardInterrupt, EOFError, Cancelled):
        console.print("\nStopped. Completed predictions remain cached.", style="yellow")
        raise SystemExit(130) from None
    except (
        OSError,
        ValueError,
        RuntimeError,
        httpx.HTTPError,
        subprocess.CalledProcessError,
        InferenceError,
    ) as error:
        console.print(f"\nCould not finish: {error}", style="red", markup=False)
        console.print("Completed predictions are cached. Fix the issue and run the same command to resume.")
        raise SystemExit(2) from error


if __name__ == "__main__":
    main()
