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
