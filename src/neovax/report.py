import json
from importlib.resources import files
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from plotly.offline import get_plotlyjs

from neovax.evaluation import assay_matrix
from neovax.io import digest, write_json


def static_figures(directory, result):
    directory.mkdir(exist_ok=True)
    selected = result["proposed"]
    fig = Figure(figsize=(10, 5), layout="constrained")
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    if selected:
        ax.barh(
            [f"{r['selection_rank']}. {r['gene']}" for r in reversed(selected)],
            [r["ranking_cost"] for r in reversed(selected)],
            color="#147d78",
        )
    else:
        ax.text(0.5, 0.5, "No eligible complete predictions", ha="center", transform=ax.transAxes)
    ax.set_xlabel(
        "Driving score (nM)"
        if result["config"]["ranking"] == "binding"
        else "Driving ranking cost (dimensionless)"
    )
    ax.set_title(f"Patient {result['config']['patient']} · proposed targets")
    for suffix in ["svg", "png"]:
        fig.savefig(directory / f"shortlist.{suffix}", dpi=160)
    rows = result["evaluations"]["historical_scores"]
    fig = Figure(figsize=(7, 5), layout="constrained")
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    for state, color in [("mutant", "#147d78"), ("reference", "#b97542")]:
        group = [r for r in rows if r["state"] == state]
        if group:
            ax.scatter(
                [r["historical_2_4_ic50_nm"] for r in group],
                [r["new_4_1_ic50_nm"] for r in group],
                label=f"{state} (n={len(group)})",
                color=color,
                alpha=0.7,
            )
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="Historical NetMHCpan 2.4 IC50 (nM)",
        ylabel="NetMHCpan 4.1 BA IC50 (nM)",
    )
    if rows:
        ax.legend()
    for suffix in ["svg", "png"]:
        fig.savefig(directory / f"historical_scores.{suffix}", dpi=160)


def render_report(run_directory, result=None):
    run_directory = Path(run_directory)
    if result is None:
        result = json.loads((run_directory / "run_data.json").read_text())
    result = {**result, "assay_matrix": assay_matrix(result["assays"])}
    static_figures(run_directory / "figures", result)
    # Large tables use compact rows in JSON and render only the current page.
    environment = Environment(
        loader=FileSystemLoader(str(files("neovax").joinpath("templates"))),
        autoescape=select_autoescape(["html"]),
    )
    compact = {k: v for k, v in result.items() if k != "pairs"}
    chosen = {r["driving_epitope"] for r in result["proposed"]}
    compact["selected_pairs"] = [r for r in result["pairs"] if r["peptide"] in chosen]
    compact["scatter_pairs"] = [
        {
            k: r[k]
            for k in [
                "peptide",
                "allele",
                "ic50_nm",
                "deepimmuno_score",
                "binding_fractional_rank",
                "immunogenicity_fractional_rank",
                "ranking_cost",
            ]
        }
        for r in result["pairs"]
        if r["deepimmuno_score"] is not None
    ]
    html = environment.get_template("report.html").render(result=compact, plotly_js=get_plotlyjs())
    path = run_directory / "report.html"
    path.write_text(html)
    manifest_path = run_directory / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if "exports" in manifest:
            manifest["exports"]["report.html"] = digest(path)
            for figure in (run_directory / "figures").iterdir():
                manifest["exports"][figure.relative_to(run_directory).as_posix()] = digest(figure)
            write_json(manifest_path, manifest)
    return path
