"""Figures from the simulation study.

Run after scripts/simulation_study.py:
    python scripts/make_figures.py
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
RESULTS = ROOT / "results"


def _need(name):
    path = RESULTS / name
    if not path.exists():
        raise SystemExit(f"missing {path}; run scripts/simulation_study.py first")
    return pd.read_csv(path)


def fig_null():
    """False positives when nothing changes."""
    df = _need("sim_null.csv")
    n_genes = int(df["n_genes"].iloc[0])

    fig, ax = plt.subplots(figsize=(5.5, 4))
    x = np.arange(len(df))
    ax.bar(x - 0.2, df["per_cell"], width=0.4, label="per-cell Wilcoxon",
           color="C3", edgecolor="black", linewidth=0.5)
    ax.bar(x + 0.2, df["pseudobulk"], width=0.4, label="pseudobulk",
           color="C0", edgecolor="black", linewidth=0.5)
    ax.axhline(0.05 * n_genes, ls="--", lw=1, color="grey",
               label="5% of genes tested")
    ax.set(xlabel="simulation replicate", ylabel="genes called significant",
           title=f"False positives on null data ({n_genes} genes, nothing changed)")
    ax.title.set_fontsize(10)
    ax.set_xticks(x)
    ax.set_xticklabels(df["rep"] + 1)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "null_false_positives.png", dpi=150)
    plt.close(fig)


def fig_power():
    """Power is not the problem; the false discovery proportion is."""
    df = _need("sim_power.csv")
    agg = df.groupby("method")[["tp", "fp", "power", "fdp"]].mean()

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))

    methods = list(agg.index)
    x = np.arange(len(methods))
    axes[0].bar(x - 0.2, agg["tp"], width=0.4, label="true positives",
                color="C0", edgecolor="black", linewidth=0.5)
    axes[0].bar(x + 0.2, agg["fp"], width=0.4, label="false positives",
                color="C3", edgecolor="black", linewidth=0.5)
    axes[0].axhline(50, ls="--", lw=1, color="grey", label="50 genes truly changed")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(methods)
    axes[0].set_ylabel("genes")
    axes[0].set_title("What each method calls", fontsize=10)
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].bar(x, agg["fdp"], width=0.5,
                color=["C3" if m == "per_cell" else "C0" for m in methods],
                edgecolor="black", linewidth=0.5)
    axes[1].axhline(0.05, ls="--", lw=1, color="grey", label="nominal 5%")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(methods)
    axes[1].set_ylabel("observed false discovery proportion")
    axes[1].set_title("Both have ~100% power", fontsize=10)
    axes[1].legend(frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(FIG / "power_vs_fdp.png", dpi=150)
    plt.close(fig)


def fig_sample_sweep():
    """The diagnosis: per-cell is fine until sample-level variance exists."""
    df = _need("sim_sample_sweep.csv")

    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.plot(df["sample_sd"], df["per_cell"], "o-", color="C3",
            label="per-cell Wilcoxon")
    ax.plot(df["sample_sd"], df["pseudobulk"], "o-", color="C0",
            label="pseudobulk")
    ax.set(xlabel="sample-level variation (donor x condition, log SD)",
           ylabel="false positives (of 500 null genes)",
           title="Per-cell testing is valid only when samples differ by\n"
                 "sampling noise alone")
    ax.title.set_fontsize(10)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "sample_variation_sweep.png", dpi=150)
    plt.close(fig)


def fig_normalisation():
    """Composition bias gets worse with more donors, not better."""
    df = _need("sim_normalisation.csv")

    fig, ax = plt.subplots(figsize=(5.5, 4))
    for norm, g in df.groupby("normalisation"):
        g = g.sort_values("n_donors")
        ax.plot(g["n_donors"], g["fdp"], "o-",
                color="C3" if norm == "cpm" else "C0",
                label="total-count CPM" if norm == "cpm"
                      else "median-of-ratios")
    ax.axhline(0.05, ls="--", lw=1, color="grey", label="nominal 5%")
    ax.set(xlabel="number of donors",
           ylabel="observed false discovery proportion",
           title="Composition bias does not average out")
    ax.title.set_fontsize(10)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "normalisation.png", dpi=150)
    plt.close(fig)


def fig_real_null():
    """The same null experiment on the real dataset."""
    path = RESULTS / "false_positives.csv"
    if not path.exists():
        print("  (skipping real-data figure; run false_positive_test.py)")
        return
    df = pd.read_csv(path)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df["split"] + 1, df["per_cell"], "o-", color="C3",
            label="per-cell Wilcoxon")
    ax.plot(df["split"] + 1, df["pseudobulk"], "o-", color="C0",
            label="pseudobulk")
    n_genes = int(df["n_genes_sc"].iloc[0])
    ax.axhline(0.05 * n_genes, ls="--", lw=1, color="grey",
               label=f"5% of {n_genes} genes")
    ax.set(xlabel="random donor split (control cells only)",
           ylabel="genes called significant",
           title="Real data, Kang 2018 CD14+ monocytes:\n"
                 "donor labels shuffled, so every call is a false positive")
    ax.title.set_fontsize(10)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "real_null.png", dpi=150)
    plt.close(fig)


def fig_isg_ranking():
    """Ranking by q-value vs by effect size, on the real response."""
    path = RESULTS / "real_analysis_summary.csv"
    if not path.exists():
        print("  (skipping ISG figure; run real_analysis.py)")
        return
    df = pd.read_csv(path).sort_values("n_cells", ascending=False)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    x = np.arange(len(df))
    labels = [c.replace(" cells", "").replace(" Monocytes", " Mono")
              for c in df["cell_type"]]

    axes[0].bar(x - 0.2, df["median_isg_rank_by_q"], width=0.4,
                label="ranked by q-value", color="C3",
                edgecolor="black", linewidth=0.5)
    axes[0].bar(x + 0.2, df["median_isg_rank_by_effect"], width=0.4,
                label="ranked by effect size", color="C0",
                edgecolor="black", linewidth=0.5)
    axes[0].set_yscale("log")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    axes[0].set_ylabel("median rank of known interferon genes")
    axes[0].set_title("Where the real biology lands in the list", fontsize=10)
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].bar(x - 0.2, df["isgs_top12_by_q"], width=0.4,
                label="ranked by q-value", color="C3",
                edgecolor="black", linewidth=0.5)
    axes[1].bar(x + 0.2, df["isgs_top12_by_effect"], width=0.4,
                label="ranked by effect size", color="C0",
                edgecolor="black", linewidth=0.5)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    axes[1].set_ylabel("known interferon genes in the top 12")
    axes[1].set_title("Of 12 canonical ISGs", fontsize=10)
    axes[1].legend(frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(FIG / "isg_ranking.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    fig_null()
    fig_power()
    fig_sample_sweep()
    fig_normalisation()
    fig_real_null()
    fig_isg_ranking()
    print(f"figures written to {FIG}")
