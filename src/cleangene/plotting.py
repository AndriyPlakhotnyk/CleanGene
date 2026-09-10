from __future__ import annotations
import csv
from pathlib import Path

def plot_presence_absence(matrix_tsv: Path, outdir: Path, organism: str, max_cluster: int = 2000) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    outdir.mkdir(parents=True, exist_ok=True)
    with matrix_tsv.open(newline="") as h:
        reader=csv.reader(h,delimiter="\t"); header=next(reader); isolates=header[1:]; rows=[r for r in reader if r]
    genes=[r[0] for r in rows]; data=np.array([[int(x) for x in r[1:]] for r in rows],dtype=np.uint8) if rows else np.zeros((0,0),dtype=np.uint8)
    if data.size:
        order=np.argsort(-data.mean(axis=1),kind="stable"); data=data[order]; genes=[genes[i] for i in order]
        if data.shape[1] <= max_cluster and data.shape[1] > 1:
            try:
                from scipy.cluster.hierarchy import leaves_list, linkage
                from scipy.spatial.distance import pdist
                col_order=leaves_list(linkage(pdist(data.T,metric="jaccard"),method="average"))
            except Exception:
                col_order=np.argsort(-data.mean(axis=0),kind="stable")
        else:
            col_order=np.argsort(-data.mean(axis=0),kind="stable")
        data=data[:,col_order]
    prevalence=data.mean(axis=1) if data.size else np.array([])
    fig,(ax,curve)=plt.subplots(1,2,figsize=(10,7),gridspec_kw={"width_ratios":[12,2]})
    if data.size:
        ax.imshow(data,aspect="auto",interpolation="nearest",cmap=plt.matplotlib.colors.ListedColormap(["#f7f5ef","#333333"]),rasterized=True)
    else:
        ax.text(0.5,0.5,"No genes",ha="center",va="center",transform=ax.transAxes)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(organism.replace("_"," "),fontstyle="italic")
    curve.plot(prevalence,range(len(prevalence)),color="#333333",linewidth=1.2)
    curve.set_ylim(len(prevalence),0); curve.set_xticks([0,1]); curve.set_yticks([])
    curve.set_xlabel("prevalence")
    fig.text(0.5,0.02,f"{data.shape[1]} isolates | {data.shape[0]} genes",ha="center")
    fig.tight_layout(rect=(0,0.04,1,1))
    fig.savefig(outdir/"pangenome_presence_absence.svg",dpi=200)
    fig.savefig(outdir/"pangenome_presence_absence.png",dpi=200)
    plt.close(fig)


def plot_decision_upset(data_tsv: Path, outdir: Path, group: str, max_columns: int = 30) -> None:
    """UpSet-style decision columns with direction-colored call counts.

    A column is a decision_reason plus its recorded metric set/provenance. This
    preserves reason labels even when several rules use the same metric set.
    Large reports are paginated without dropping rare decisions.
    """
    from collections import Counter
    import math
    import textwrap
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from .util import read_tsv

    counts = Counter()
    for row in read_tsv(data_tsv):
        metrics = tuple(filter(None, row["decision_metrics"].split(";")))
        counts[(row["decision_reason"], metrics, row["metrics_provenance"], row["change"])] += int(row["n_calls"])
    keys = sorted({key[:3] for key in counts}, key=lambda key: (-sum(counts[(*key, direction)] for direction in ("0->1", "1->0")), key))
    outdir.mkdir(parents=True, exist_ok=True)
    # Remove stale extra pages if a resumed report now has fewer decisions.
    for suffix in ("png", "svg"):
        for old in outdir.glob(f"decision_reason_upset_page_*.{suffix}"): old.unlink()
    total = sum(counts.values())
    pages = max(1, math.ceil(len(keys) / max_columns))
    colors = {"0->1": "#0072B2", "1->0": "#D55E00"}
    for page in range(pages):
        subset = keys[page * max_columns:(page + 1) * max_columns]
        stem = "decision_reason_upset" if page == 0 else f"decision_reason_upset_page_{page + 1:02d}"
        if not subset:
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.axis("off")
            ax.text(.5, .55, "No gene calls changed", ha="center", va="center", fontsize=19)
            ax.text(.5, .42, "No 0→1 or 1→0 decisions to display", ha="center", color="#555555")
            fig.suptitle(f"{group} · read-validation decisions")
        else:
            metrics = sorted({metric for _, members, _ in subset for metric in members})
            fig = plt.figure(figsize=(max(11, min(26, 5 + .72 * len(subset))), 9))
            grid = fig.add_gridspec(2, 2, width_ratios=[3.0, max(5, .72 * len(subset))], height_ratios=[3, 2], hspace=.12, wspace=.18)
            bars = fig.add_subplot(grid[0, 1]); matrix = fig.add_subplot(grid[1, 1], sharex=bars)
            sets = fig.add_subplot(grid[1, 0], sharey=matrix)
            legend = fig.add_subplot(grid[0, 0]); legend.axis("off")
            positions = list(range(len(subset))); bottom = [0] * len(subset)
            for direction in ("0->1", "1->0"):
                heights = [counts[(*key, direction)] for key in subset]
                bars.bar(positions, heights, bottom=bottom, color=colors[direction], label=direction.replace("->", "→"), width=.65)
                bottom = [a + b for a, b in zip(bottom, heights)]
            bars.set_ylabel("Changed gene calls"); bars.yaxis.set_major_locator(MaxNLocator(integer=True))
            bars.tick_params(axis="x", bottom=False, labelbottom=False)
            bars.spines[["top", "right"]].set_visible(False)
            for x, count in zip(positions, bottom): bars.text(x, count, str(count), ha="center", va="bottom", fontsize=8)
            bars.margins(y=.18)
            for x, (_, members, _) in enumerate(subset):
                matrix.scatter([x] * len(metrics), range(len(metrics)), color="#e3e7eb", s=22, zorder=1)
                active = [i for i, metric in enumerate(metrics) if metric in members]
                if active:
                    matrix.plot([x, x], [min(active), max(active)], color="#263746", linewidth=1.4)
                    matrix.scatter([x] * len(active), active, color="#263746", s=42, zorder=3)
            matrix.set_ylim(len(metrics) - .5, -.5)
            matrix.set_yticks(range(len(metrics))); matrix.tick_params(axis="y", left=False, labelleft=False)
            labels = [textwrap.fill(reason, 20) + (" *" if provenance != "recorded" else "") for reason, _, provenance in subset]
            matrix.set_xticks(positions, labels, rotation=55, ha="right", fontsize=9)
            matrix.set_xlabel("decision_reason")
            matrix.spines[["top", "right", "left"]].set_visible(False)
            set_counts = [sum(counts[(*key, direction)] for key in subset if metric in key[1] for direction in ("0->1", "1->0")) for metric in metrics]
            sets.barh(range(len(metrics)), set_counts, color="#8399ab", height=.6)
            sets.set_yticks(range(len(metrics)), [metric.replace("_", " ") for metric in metrics], fontsize=9)
            sets.invert_xaxis(); sets.set_xlabel("Calls using metric")
            sets.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=3)); sets.spines[["top", "right", "left"]].set_visible(False)
            handles, labels = bars.get_legend_handles_labels()
            legend.legend(handles, labels, title="Call change", loc="upper left", frameon=False)
            legend.text(0, .5, "Dots: inputs to the deciding rule\n(not thresholds passed).\n\n* Metric set inferred from legacy\nevidence or unavailable.", fontsize=9, va="top", color="#444444")
            fig.suptitle(f"{group} · read-validation decisions\n{total:,} changed calls · page {page + 1}/{pages}", fontsize=15)
            fig.subplots_adjust(left=.18, right=.98, bottom=.29, top=.89)
        for suffix in ("png", "svg"):
            fig.savefig(outdir / f"{stem}.{suffix}", dpi=180, bbox_inches="tight")
        plt.close(fig)
