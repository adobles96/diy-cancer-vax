"""Portable PNG, SVG and PDF charts for the retrospective benchmark."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

from neovax.cohort import METHODS

LABELS = {"netmhcpan": "NetMHCpan", "deepimmuno": "DeepImmuno", "mean_ensemble": "Ensemble"}
COLORS = ("#176E60", "#AC562E", "#64736D")


def grouped_chart(output, name, title, subtitle, ticks, values, labels, note, percentage=False):
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 11, "svg.fonttype": "none"}):
        figure, axis = plt.subplots(figsize=(14, 7), facecolor="#FAF8F3")
        figure.subplots_adjust(left=0.075, right=0.98, bottom=0.25, top=0.72)
        axis.set_facecolor("#FAF8F3")
        figure.text(0.075, 0.91, title, size=23, weight="bold", color="#23332F")
        figure.text(0.075, 0.85, subtitle, size=11, color="#64736D")
        axis.axvspan(5.5, len(ticks) - 0.5, color="#E8EEE9", zorder=0)
        for offset, (method, color) in enumerate(zip(METHODS, COLORS)):
            bars = axis.bar(
                [i + (offset - 1) * 0.28 for i in range(len(ticks))],
                values[method],
                width=0.25,
                label=LABELS[method],
                color=color,
                zorder=3,
            )
            axis.bar_label(bars, labels=labels[method], padding=5, fontsize=9)
        axis.set_xticks(range(len(ticks)), ticks)
        axis.set(
            xlim=(-0.6, len(ticks) - 0.4),
            ylim=(0, 100 if percentage else 1),
            ylabel="Vaccine recall@20" if percentage else "ROC AUC",
        )
        if percentage:
            axis.yaxis.set_major_formatter(PercentFormatter(100, decimals=0))
        else:
            axis.axhline(0.5, color="#88948E", linestyle="--", linewidth=1)
        axis.grid(axis="y", color="#D5DED6", linewidth=0.7, zorder=0)
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.spines["bottom"].set_color("#D5DED6")
        axis.tick_params(length=0, pad=10)
        axis.legend(loc="lower left", bbox_to_anchor=(-0.005, 1.05), ncol=3, frameon=False)
        figure.text(0.075, 0.10, note, size=10, color="#64736D")
        for extension in ("png", "svg", "pdf"):
            figure.savefig(output / f"{name}.{extension}", dpi=180, facecolor=figure.get_facecolor())
        plt.close(figure)


def render_recall(output, metrics, pooled):
    values, labels = {}, {}
    ticks = []
    for method in METHODS:
        rows = [
            next(
                r
                for r in metrics
                if r["patient"] == patient and r["method"] == method and r["endpoint"] == "vaccine"
            )
            for patient in range(1, 7)
        ]
        overall = next(r for r in pooled if r["method"] == method and r["endpoint"] == "vaccine")
        values[method] = [100 * r["recall"] for r in rows] + [100 * overall["micro_recall"]]
        labels[method] = [str(r["hits"]) for r in rows] + [str(overall["hits"])]
        ticks = [f"Patient {r['patient']}\n(n = {r['denominator']})" for r in rows]
        ticks.append(f"Overall\n(n = {overall['denominator']})")
    grouped_chart(
        output,
        "vaccine-recall-at-20",
        "Recovery of actual vaccine targets at top 20",
        "Recall@20 = recovered vaccine components / all actual vaccine components",
        ticks,
        values,
        labels,
        "Labels show recovered component counts. Overall pools counts across patients; overlap is matched by mutation.",
        percentage=True,
    )


def render_auc(output, metrics, aggregate):
    values, labels = {}, {}
    ticks = []
    for method in METHODS:
        rows = [
            next(r for r in metrics if r["patient"] == patient and r["method"] == method)
            for patient in range(1, 7)
        ]
        overall = next(r for r in aggregate if r["method"] == method)
        scores = [r["roc_auc"] for r in rows]
        scores.extend([overall["macro_auc"], overall["within_patient_pair_weighted_auc"]])
        ticks = [f"Patient {r['patient']}\n(n = {r['n']})" for r in rows]
        ticks.extend(["Patient\nmean", "Within-patient\npair-weighted"])
        values[method] = [score if score is not None else 0 for score in scores]
        labels[method] = [f"{score:.3f}" if score is not None else "Undefined" for score in scores]
    grouped_chart(
        output,
        "peptide-recognition-roc-auc",
        "Predicting peptide recognition: ROC AUC",
        "Reconstructed, scored vaccine mutations with week-16 labels; no top-20 cutoff",
        ticks,
        values,
        labels,
        "Patient mean gives equal weight; pair-weighted compares targets only within patients.",
    )
