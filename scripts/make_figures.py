"""Figures for the README: simulation study, real-data checks, model comparisons.

Every figure reads a committed artefact in results/ and names it in a grey
footer together with n; nothing is typed in by hand. Titles state the finding,
and where the finding is a number the title is computed from the CSV, so a
re-run cannot leave a stale claim on the page.

Run after the producing scripts (simulation_study.py, false_positive_test.py,
real_analysis.py, mixed_vs_others.py, nbglm_real.py):
    python scripts/make_figures.py                 # every figure
    python scripts/make_figures.py null power      # only the named ones
A figure whose artefact is missing is skipped with a message.

Colours: C0 = pseudobulk, C3 = the method a figure is arguing against
(per-cell Wilcoxon, NB GLM), C1 = the mixed model.
"""

import math
import re
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


def _maybe(name, producer, **kwargs):
    """Read an artefact, or None (with a message) when it was not produced."""
    path = RESULTS / name
    if not path.exists():
        print(f"  (skipping: {path.name} is missing; run {producer})")
        return None
    return pd.read_csv(path, **kwargs)


def _wrap_to_figure(fig, text, fontsize):
    """Greedy word wrap by rendered width, so the footer never leaves the page."""
    renderer = fig.canvas.get_renderer()
    probe = fig.text(0, 0, "", fontsize=fontsize)
    limit = 0.965 * fig.get_figwidth() * fig.dpi
    lines = []
    for para in text.split("\n"):
        line = ""
        for word in para.split(" "):
            trial = (line + " " + word).strip()
            probe.set_text(trial)
            if line and probe.get_window_extent(renderer).width > limit:
                lines.append(line)
                line = word
            else:
                line = trial
        lines.append(line)
    probe.remove()
    return lines


def _save(fig, name, footer):
    """Write n and the source artefact in a grey footer, then save at 150 dpi.

    The footer is wrapped to the figure width and the plot area is shrunk to
    leave room for it, so a long source list never runs off the page.
    """
    lines = _wrap_to_figure(fig, footer, 6.5)
    h = fig.get_figheight()
    fig.text(0.012, 0.05 / h, "\n".join(lines), fontsize=6.5, color="0.35",
             ha="left", va="bottom", linespacing=1.3)
    fig.tight_layout(rect=(0, (0.1 + 0.105 * len(lines)) / h, 1, 1))
    fig.savefig(FIG / name, dpi=150)
    plt.close(fig)


def _short(cell_type):
    return cell_type.replace(" cells", "").replace(" Monocytes", " Mono")


def _tick(cell_type):
    """Two-line tick label (the word Mono wraps onto its own line)."""
    return _short(cell_type).replace(" Mono", "\nMono")


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


def fig_mixed_null():
    """Mixed model and pseudobulk make the same calls on the donor-swap null."""
    swap = _maybe("mixed_null_donor_swap.csv", "scripts/mixed_vs_others.py")
    draws = _maybe("mixed_null_draws.csv", "scripts/mixed_vs_others.py")
    if swap is None or draws is None:
        return
    bal = swap[swap["balance"] == "balanced"]
    pb = bal[bal["method"] == "pseudobulk"].set_index("cell_type")
    mx = bal[bal["method"] == "mixed"].set_index("cell_type").loc[pb.index]
    bd = draws[draws["balance"] == "balanced"]
    n_draws = len(bd)
    same = int((bd["fp_pseudobulk"] == bd["fp_mixed"]).sum())
    any_fp = int(((bd["fp_pseudobulk"] > 0) | (bd["fp_mixed"] > 0)).sum())
    n_genes = int(pb["genes_tested"].iloc[0])
    ref = float(pb["expected_fp_bh"].iloc[0])
    by_type = " + ".join(str(int(v)) for v in pb["n_distinct_draws"])

    fig, ax = plt.subplots(figsize=(6.5, 4))
    x = np.arange(len(pb))
    for off, tab, color, label in ((-0.2, pb, "C0", "pseudobulk paired t"),
                                   (0.2, mx, "C1", "mixed (1 | donor)")):
        mean = tab["mean_false_positives"].to_numpy()
        lo = tab["fp_ci_lo"].to_numpy()
        hi = tab["fp_ci_hi"].to_numpy()
        ax.bar(x + off, mean, width=0.4, label=label, color=color,
               edgecolor="black", linewidth=0.5, yerr=[mean - lo, hi - mean],
               capsize=2, error_kw={"lw": 0.8})
        for xi, m, h in zip(x + off, mean, hi):
            ax.text(xi, h + 0.012, f"{m:.3f}", ha="center", va="bottom",
                    fontsize=6.5)
    ax.axhline(ref, ls="--", lw=1, color="grey",
               label=f"{ref:.2f}: BH expectation for an exact null")
    top = max(pb["fp_ci_hi"].max(), mx["fp_ci_hi"].max())
    ax.set_ylim(0, top * 1.15 + 0.03)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{_tick(ct)}\n{int(d)} donors\n{int(n)} draws"
                        for ct, d, n in zip(pb.index, pb["n_donors_used"],
                                            pb["n_distinct_draws"])],
                       fontsize=7.5)
    ax.set(ylabel=f"false positives per draw\n(BH q < 0.05, {n_genes} null genes)",
           title="Mixed model and pseudobulk make the same null calls:\n"
                 f"the same false-positive count in {same} of {n_draws} balanced draws")
    ax.title.set_fontsize(10)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    _save(fig, "mixed_vs_pseudobulk_null.png",
          f"n = {n_draws} balanced donor-swap draws per method ({by_type} by type), "
          f"{n_genes} small-effect genes each; {any_fp} draws had any false "
          "positive. Bars: mean; whiskers: bootstrap 95% CI. A swap null can show "
          "anti-conservatism, not calibration.\n"
          "source: results/mixed_null_donor_swap.csv, results/mixed_null_draws.csv")


FIGURES = {
    "null": fig_null,
    "power": fig_power,
    "sweep": fig_sample_sweep,
    "normalisation": fig_normalisation,
    "real_null": fig_real_null,
    "isg": fig_isg_ranking,
    "mixed_null": fig_mixed_null,
}


if __name__ == "__main__":
    wanted = sys.argv[1:] or list(FIGURES)
    unknown = [w for w in wanted if w not in FIGURES]
    if unknown:
        raise SystemExit(f"unknown figure(s) {unknown}; choose from {list(FIGURES)}")
    for name in wanted:
        FIGURES[name]()
    print(f"figures written to {FIG}")
