"""Run: .venv/bin/python scripts/cohort_auc.py (uses the saved cohort scores)."""

import argparse
import csv
import datetime
import html
import json
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from neovax.assay_auc import assay_labels, roc_metrics
from neovax.cohort import METHODS
from neovax.data import patient_data
from neovax.io import digest, write_csv, write_json

LABELS = {"netmhcpan": "NetMHCpan", "deepimmuno": "DeepImmuno", "mean_ensemble": "Mean ensemble"}
METHODOLOGY = """Evaluation is restricted to actual vaccinated mutations with independently reconstructed
and scored candidates and usable week-16 peptide-response labels. Each genomic mutation counts once.
Positive = any positive CD4 or CD8 peptide response, ex vivo or after expansion, from any associated IMP.
Specifically: IFN-γ ELISpot against peptide-pulsed autologous antigen-presenting cells, using CD8
cells after pre-stimulation, CD4 responses directly ex vivo, or CD4 cells after pre-stimulation.
We use the authors' positive/negative calls; we do not impose a new threshold on spot counts.
Negative = both CD4 and CD8 have measurements for each associated IMP, and all eligible measurements
are explicitly negative. Missing, not-determined, and ambiguous observations do not become negatives.
Minigene, tumor, and later-timepoint observations do not enter the label. This is observed response
under the reported assays, not proof of biological non-immunogenicity.
Scores are unchanged from the previous cohort analysis: minimum NetMHCpan IC50 or maximum DeepImmuno
score over all reconstructed mutation-containing 9/10-mer–HLA-A/B pairs for each mutation.
The candidate set inherits the reference-proteome exact-match filter recorded in the input provenance.
The ensemble averages the two within-patient fractional ranks over the full discovery candidate set;
it is not refitted or renormalized on the vaccine evaluation set. No top-K or binding cutoff is used.
These are mutation-level scores; the best short peptide need not be the experimentally tested peptide.
AUC is the fraction of positive–negative target pairs correctly ordered, giving score ties half credit.
ROC curves group equal scores; genomic-ID tie breaks from shortlist generation are not used for AUC.
Across patients, report the equally weighted mean AUC and a pair-weighted within-patient AUC.
The latter pools concordance counts but never compares a target from one patient with another patient's
target. It is not a conventional AUC computed by concatenating all patient scores.
The estimand is discrimination conditional on vaccination and successful reconstruction. Historical
binding-based selection still limits generalization to all mutations. These class-I models do not
directly predict CD4 responses. HLA rescue and training-overlap limitations from the prior run persist.
With only 9–15 evaluated mutations per patient, these descriptive estimates have substantial uncertainty."""


def render(output, metrics, aggregate):
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), layout="constrained")
    colors = ["#3266a8", "#bb6233", "#38816a"]
    body = []
    for patient, ax in zip(range(1, 7), axes.flat):
        group = [r for r in metrics if r["patient"] == patient]
        for row, color in zip(group, colors):
            curve = row["roc_curve"]
            auc_text = f"{row['roc_auc']:.3f}" if row["roc_auc"] is not None else "undefined"
            ax.plot(
                [p["fpr"] for p in curve],
                [p["tpr"] for p in curve],
                color=color,
                label=f"{LABELS[row['method']]} ({auc_text})",
            )
        ax.plot([0, 1], [0, 1], color="#999", linestyle="--", linewidth=1)
        ax.set(
            title=f"P{patient}: {group[0]['positives']} positive, {group[0]['negatives']} negative",
            xlabel="False positive rate",
            ylabel="True positive rate",
            xlim=(0, 1),
            ylim=(0, 1.02),
        )
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(fontsize=8, loc="lower right", frameon=False)
        cells = "".join(
            f"<td>{r['roc_auc']:.3f}</td>" if r["roc_auc"] is not None else "<td>Undefined</td>"
            for r in group
        )
        body.append(
            f"<tr><th>P{patient}</th><td>{group[0]['positives']} / {group[0]['negatives']}</td>{cells}</tr>"
        )
    fig.suptitle("Immunogenicity discrimination within the reconstructable vaccine set", fontsize=15)
    fig.savefig(output / "roc-curves.png", dpi=180)
    fig.savefig(output / "roc-curves.svg")
    plt.close(fig)
    total_positive = sum(r["positives"] for r in metrics if r["method"] == METHODS[0])
    total_negative = sum(r["negatives"] for r in metrics if r["method"] == METHODS[0])
    for key, label in [
        ("macro_auc", "Mean across patients"),
        ("within_patient_pair_weighted_auc", "Within-patient pair-weighted"),
    ]:
        cells = "".join(
            f"<td>{r[key]:.3f}</td>" if r[key] is not None else "<td>Undefined</td>" for r in aggregate
        )
        body.append(f"<tr><th>{label}</th><td>{total_positive} / {total_negative}</td>{cells}</tr>")
    page = f"""<!doctype html><html lang="en"><meta charset="utf-8"><title>NeoVax immunogenicity ROC AUC</title>
<style>body{{font:16px/1.5 system-ui;max-width:1100px;margin:40px auto;padding:0 24px;color:#223042}}table{{border-collapse:collapse;width:100%}}td,th{{text-align:left;padding:12px;border-bottom:1px solid #ddd}}pre{{white-space:pre-wrap;background:#f2f5f8;padding:20px}}img{{width:100%}}</style>
<h1>Immunogenicity ROC AUC</h1><p>Actual vaccinated, reconstructable mutations · Any positive week-16 CD4/CD8 peptide response · No top-K cutoff</p>
<table><thead><tr><th>Patient</th><th>Positive / negative</th><th>NetMHCpan</th><th>DeepImmuno</th><th>Mean ensemble</th></tr></thead><tbody>{"".join(body)}</tbody></table>
<p>AUC 0.5 is chance-level ordering; 1 is perfect ordering. AUC is not the fraction of binary predictions classified correctly.</p>
<h2>What counts as immunogenic?</h2>
<p>A mutation is positive if <strong>any associated peptide has an author-reported positive week-16 IFN-γ ELISpot response</strong> in one of these measurements:</p>
<table><thead><tr><th>Response</th><th>Condition</th><th>Source columns</th></tr></thead><tbody>
<tr><td>CD8</td><td>After peptide pre-stimulation / cell expansion; response to peptide-pulsed autologous antigen-presenting cells</td><td>Ott 2017 supplementary workbook 5, O; Hu 2021 Dataset 4a, O</td></tr>
<tr><td>CD4</td><td>Directly ex vivo, without prior cell expansion; response to added peptide</td><td>Hu 2021 Dataset 4b, J</td></tr>
<tr><td>CD4</td><td>After peptide pre-stimulation / cell expansion; response to peptide-pulsed autologous antigen-presenting cells</td><td>Hu 2021 Dataset 4b, K</td></tr>
</tbody></table>
<p>IFN-γ ELISpot detects cells secreting interferon-gamma in response to peptide exposure. We use the authors' binary calls, without applying a new spot-count threshold. Repeated records of the same assay are deduplicated.</p>
<p>Minigene recognition, autologous tumor recognition, HLA blocking, tetramer staining, degranulation, pooled cytokine assays, and later/post-checkpoint-treatment observations do not enter this endpoint. A positive label does not require mutant-versus-wild-type selectivity or demonstrated tumor-cell killing.</p>
<p>The AUC set contains {total_positive} positive and {total_negative} negative reconstructed vaccine mutations. Negatives require explicit negative measurements in both CD4 and CD8; missing or ambiguous measurements are excluded.</p>
<p>Sources: <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5577644/">Ott et al. 2017</a> and <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8273876/">Hu et al. 2021</a>.</p>
<img src="roc-curves.svg" alt="ROC curves and AUC for each of six patients">
<h2>Methods and interpretation</h2><pre>{html.escape(METHODOLOGY)}</pre>
<p><a href="metrics.csv">Per-patient metrics</a> · <a href="labeled-targets.csv">Target scores and assay evidence</a> · <a href="excluded.csv">Excluded components</a> · <a href="analysis.json">Provenance and full results</a></p></html>"""
    (output / "report.html").write_text(page)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="data/outputs/cohort-recall")
    parser.add_argument("--output-dir", default="data/outputs/cohort-auc")
    args = parser.parse_args()
    source, output = Path(args.input_dir).resolve(), Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    hashes = json.loads((source / "artifact-hashes.json").read_text())
    inputs = ["analysis.json", "target-rankings.csv", "vaccine-component-audit.csv"]
    for name in inputs:
        if digest(source / name) != hashes[name]:
            raise ValueError(f"Input file changed since the cohort run: {name}")
    provenance = json.loads((source / "analysis.json").read_text())
    with (source / "target-rankings.csv").open() as handle:
        targets = list(csv.DictReader(handle))
    with (source / "vaccine-component-audit.csv").open() as handle:
        components = list(csv.DictReader(handle))
    for target in targets:
        for method in METHODS:
            target[f"{method}_cost"] = float(target[f"{method}_cost"])
    metrics, all_labeled, all_excluded = [], [], []
    for patient in range(1, 7):
        config = next(p["config"] for p in provenance["patients"] if p["patient"] == patient)
        _, _, _, observations = patient_data(config["data_dir"], patient)
        labeled, excluded = assay_labels(
            [r for r in components if int(r["patient"]) == patient],
            [r for r in targets if int(r["patient"]) == patient],
            observations,
        )
        all_labeled.extend(labeled)
        all_excluded.extend(excluded)
        for method in METHODS:
            metrics.append(
                {
                    "patient": patient,
                    "method": method,
                    **roc_metrics([r["label"] for r in labeled], [r[f"{method}_cost"] for r in labeled]),
                }
            )
    aggregate = []
    for method in METHODS:
        group = [r for r in metrics if r["method"] == method and r["roc_auc"] is not None]
        pairs = sum(r["positive_negative_pairs"] for r in group)
        credits = sum(r["correctly_ordered_pairs"] + 0.5 * r["tied_pairs"] for r in group)
        aggregate.append(
            {
                "method": method,
                "patients_with_defined_auc": len(group),
                "macro_auc": statistics.mean(r["roc_auc"] for r in group) if group else None,
                "within_patient_pair_weighted_auc": credits / pairs if pairs else None,
                "correct_pair_credits": credits,
                "within_patient_pairs": pairs,
            }
        )
    result = {
        "created_at": datetime.datetime.now(datetime.UTC).isoformat(),
        "methodology": METHODOLOGY,
        "inputs": {str(source / name): hashes[name] for name in inputs},
        "score_provenance": provenance,
        "labeled_targets": len(all_labeled),
        "metrics": metrics,
        "aggregate": aggregate,
        "excluded": all_excluded,
    }
    write_json(output / "analysis.json", result)
    write_csv(output / "metrics.csv", [{k: v for k, v in r.items() if k != "roc_curve"} for r in metrics])
    write_csv(output / "aggregate.csv", aggregate)
    write_csv(output / "labeled-targets.csv", all_labeled)
    write_csv(output / "excluded.csv", all_excluded)
    render(output, metrics, aggregate)
    write_json(
        output / "artifact-hashes.json",
        {p.name: digest(p) for p in output.iterdir() if p.is_file() and p.name != "artifact-hashes.json"},
    )
    print(json.dumps(aggregate, indent=2))
    print(f"Report: {output / 'report.html'}")


if __name__ == "__main__":
    main()
