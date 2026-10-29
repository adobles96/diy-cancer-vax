import json
from importlib.resources import files
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape
from plotly.offline import get_plotlyjs

from neovax.evaluation import assay_matrix
from neovax.io import digest, write_json


def render_report(run_directory, result=None):
    run_directory = Path(run_directory)
    if result is None:
        record_path = run_directory / "run.json"
        if not record_path.exists():
            record_path = run_directory / "run_data.json"  # Previously saved runs remain readable.
        result = json.loads(record_path.read_text())
        if "summary" not in result:
            raise ValueError("This run did not finish producing results; rerun the pipeline first.")
    result = {**result, "assay_matrix": assay_matrix(result["assays"])}
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
    record_path = run_directory / "run.json"
    if record_path.exists():
        saved = json.loads(record_path.read_text())
        if "exports" in saved["manifest"]:
            saved["manifest"]["exports"]["report.html"] = digest(path)
            write_json(record_path, saved)
    else:
        manifest_path = run_directory / "manifest.json"
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            if "exports" in manifest:
                manifest["exports"]["report.html"] = digest(path)
                write_json(manifest_path, manifest)
    return path
