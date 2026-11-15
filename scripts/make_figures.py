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

# Replicate counts that scripts/simulation_study.py hard-codes (its CSVs store
# means only): experiment_sample_variation uses range(5), experiment_normalisation
# range(6).
SWEEP_REPS = 5
NORM_REPS = 6


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
    bh = _maybe("bh_null_expectation.csv", "scripts/bh_null_expectation.py")

    fig, ax = plt.subplots(figsize=(5.5, 4))
    x = np.arange(len(df))
    ax.bar(x - 0.2, df["per_cell"], width=0.4, label="per-cell Wilcoxon",
           color="C3", edgecolor="black", linewidth=0.5)
    ax.bar(x + 0.2, df["pseudobulk"], width=0.4, label="pseudobulk paired t",
           color="C0", edgecolor="black", linewidth=0.5)
    # Reference honesty (2026-11-15): this line is 5% of genes, which is what an
    # UNCORRECTED p < 0.05 threshold would reject. The calls plotted are BH
    # q < 0.05, whose expected count under an exact null is ~0.06 genes, not 25.
    # Drawing the uncorrected number on a corrected axis set the wrong bar and
    # made the methods look closer to calibrated than they are. The line is kept
    # (it is the threshold the README argues against) but relabelled as such, and
    # the BH expectation is drawn where it belongs - at zero.
    ax.axhline(0.05 * n_genes, ls=":", lw=1, color="0.55",
               label=f"{0.05 * n_genes:.0f} = 5% if p were UNCORRECTED (not this test)")
    if bh is not None:
        ax.axhline(float(bh["mean_rejections"].iloc[0]), ls="--", lw=1.2, color="black",
                   label=f"BH expectation under an exact null: "
                         f"{bh['mean_rejections'].iloc[0]:.2f} genes")
    ax.set(xlabel="simulation replicate",
           ylabel=f"genes called significant\n(of {n_genes}, BH q < 0.05)",
           title="Nothing changed, yet per-cell Wilcoxon calls "
                 f"{df['per_cell'].mean() / n_genes:.0%} of genes;\n"
                 f"pseudobulk calls {df['pseudobulk'].mean() / n_genes:.0%}")
    ax.title.set_fontsize(10)
    ax.set_xticks(x)
    ax.set_xticklabels(df["rep"] + 1)
    # Every pseudobulk bar is exactly 0, so the ten "0" labels are all at the same
    # height - and the dotted UNCORRECTED-5% line runs at 25, close enough that at
    # this scale the digits and the line touch. Place them above the line instead
    # of beside the bar: a mask cannot help when the line's height depends on
    # n_genes, and nudging by a few genes just moves the collision.
    label_y = 0.05 * n_genes * 1.25
    for xi, v in zip(x + 0.2, df["pseudobulk"]):
        ax.text(xi, label_y, f"{v:.0f}", ha="center", va="bottom", fontsize=7,
                color="0.25")
    ax.set_ylim(0, df["per_cell"].max() * 1.5)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    note = ""
    if bh is not None:
        note = (f" Under BH an exact null gives a mean of "
                f"{bh['mean_rejections'].iloc[0]:.3f} calls, not "
                f"{0.05 * n_genes:.0f}.")
    _save(fig, "null_false_positives.png",
          f"n = {len(df)} simulated replicates, {n_genes} genes, 0 truly changed; "
          f"mean per-cell {df['per_cell'].mean():.1f} ({df['per_cell'].mean() / n_genes:.1%}), "
          f"pseudobulk {df['pseudobulk'].mean():.1f}.{note}\n"
          "source: results/sim_null.csv, results/bh_null_expectation.csv")


def fig_power():
    """Power is not the problem; the false discovery proportion is."""
    df = _need("sim_power.csv")
    agg = df.groupby("method")[["tp", "fp", "power", "fdp"]].mean()

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))

    methods = list(agg.index)
    x = np.arange(len(methods))
    n_true = int((df["tp"] + df["fn"]).iloc[0])
    axes[0].bar(x - 0.2, agg["tp"], width=0.4, label="true positives",
                color="C0", edgecolor="black", linewidth=0.5)
    axes[0].bar(x + 0.2, agg["fp"], width=0.4, label="false positives",
                color="C3", edgecolor="black", linewidth=0.5)
    for off, col in ((-0.2, "tp"), (0.2, "fp")):
        for xi, v in zip(x + off, agg[col]):
            axes[0].text(xi, v + 5, f"{v:.1f}", ha="center", va="bottom",
                         fontsize=7.5)
    axes[0].axhline(n_true, ls="--", lw=1, color="grey",
                    label=f"{n_true} genes truly changed")
    names = {"per_cell": "per-cell Wilcoxon", "pseudobulk": "pseudobulk paired t"}
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([names.get(m, m) for m in methods])
    axes[0].set_ylim(0, agg[["tp", "fp"]].to_numpy().max() * 1.3)
    axes[0].set_ylabel("genes called (mean per replicate)")
    axes[0].set_title(
        f"Both find the true genes ({agg.at['per_cell', 'tp']:.1f} and "
        f"{agg.at['pseudobulk', 'tp']:.1f} of {n_true});\n"
        f"per-cell adds {agg.at['per_cell', 'fp']:.1f} false calls, pseudobulk "
        f"{agg.at['pseudobulk', 'fp']:.1f}", fontsize=10)
    axes[0].legend(frameon=False, fontsize=8, loc="upper right")

    axes[1].bar(x, agg["fdp"], width=0.5,
                color=["C3" if m == "per_cell" else "C0" for m in methods],
                edgecolor="black", linewidth=0.5)
    for xi, v in zip(x, agg["fdp"]):
        axes[1].text(xi, v + 0.02, f"{v:.2f}", ha="center", va="bottom",
                     fontsize=9, fontweight="bold")
    axes[1].axhline(0.05, ls="--", lw=1, color="grey", label="nominal 5%")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([names.get(m, m) for m in methods])
    axes[1].set_ylim(0, 1.05)
    axes[1].set_ylabel("observed false discovery proportion")
    axes[1].set_title(
        "Observed false discovery proportion:\n"
        f"per-cell {agg.at['per_cell', 'fdp']:.2f}, pseudobulk "
        f"{agg.at['pseudobulk', 'fdp']:.2f} (power {agg.at['per_cell', 'power']:.2f} "
        f"vs {agg.at['pseudobulk', 'power']:.2f})", fontsize=10)
    axes[1].legend(frameon=False, fontsize=8, loc="upper right")

    n_genes = int((df["tp"] + df["fp"] + df["fn"] + df["tn"]).iloc[0])
    _save(fig, "power_vs_fdp.png",
          f"n = {df['rep'].nunique()} simulated replicates per method, {n_true} of "
          f"{n_genes} genes given a real 2-fold effect (scripts/simulation_study.py), "
          "BH q < 0.05; bars are means.\n"
          "source: results/sim_power.csv")


def fig_sample_sweep():
    """The diagnosis: per-cell is fine until sample-level variance exists."""
    df = _need("sim_sample_sweep.csv")

    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.plot(df["sample_sd"], df["per_cell"], "o-", color="C3",
            label="per-cell Wilcoxon")
    ax.plot(df["sample_sd"], df["pseudobulk"], "o-", color="C0",
            label="pseudobulk paired t")
    for sd, v in zip(df["sample_sd"], df["per_cell"]):
        ax.text(sd, v + 12, f"{v:.1f}", ha="center", va="bottom", fontsize=7.5,
                bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.85))
    first, last = df.iloc[0], df.iloc[-1]
    ax.set(xlabel="sample-level variation (donor x condition, log SD)",
           ylabel="false positives (of 500 null genes, BH q < 0.05)",
           title=f"Per-cell false positives rise from {first['per_cell']:.0f} to "
                 f"{last['per_cell']:.0f} as sample-level\n"
                 f"variation grows; pseudobulk stays at {df['pseudobulk'].max():.0f}")
    ax.title.set_fontsize(10)
    ax.set_ylim(-15, df["per_cell"].max() * 1.15)
    ax.legend(frameon=False, fontsize=8, loc="center right")
    _save(fig, "sample_variation_sweep.png",
          f"n = {len(df)} sample-SD values, mean of {SWEEP_REPS} simulated replicates "
          "each (scripts/simulation_study.py); 500 null genes, nothing truly changed.\n"
          "source: results/sim_sample_sweep.csv")


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
        for n, v in zip(g["n_donors"], g["fdp"]):
            # The -0.008 offset put the TOP of the median-of-ratios text on the
            # dashed nominal-5% line (0.0592 - 0.008 = 0.0512 vs 0.050), so the
            # line grazed the digits. -0.016 clears it with room to spare.
            ax.text(n, v + (0.010 if norm == "cpm" else -0.016), f"{v:.3f}",
                    ha="center", va="bottom" if norm == "cpm" else "top",
                    fontsize=7.5)
    ax.axhline(0.05, ls="--", lw=1, color="grey", label="nominal 5%")
    cpm = df[df["normalisation"] == "cpm"].sort_values("n_donors")
    mor = df[df["normalisation"] == "median_ratio"]
    ax.set(xlabel="number of donors",
           ylabel="observed false discovery proportion",
           title="Total-count CPM bias grows with donors "
                 f"({cpm['fdp'].iloc[0]:.2f} to {cpm['fdp'].iloc[-1]:.2f});\n"
                 f"median-of-ratios stays at {mor['fdp'].min():.2f}-{mor['fdp'].max():.2f}")
    ax.title.set_fontsize(10)
    ax.set_xticks(sorted(df["n_donors"].unique()))
    ax.set_ylim(0, df["fdp"].max() * 1.2)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    _save(fig, "normalisation.png",
          f"n = {NORM_REPS} simulated replicates per point (mean), 50 of 500 genes truly "
          "changed (scripts/simulation_study.py), BH q < 0.05.\n"
          "source: results/sim_normalisation.csv")


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
            label="pseudobulk paired t")
    n_genes = int(df["n_genes_sc"].iloc[0])
    # Same correction as fig_null: the plotted calls are BH q < 0.05, so a 5%-of-
    # genes line is an uncorrected-p reference and must say so (2026-11-15).
    ax.axhline(0.05 * n_genes, ls=":", lw=1, color="0.55",
               label=f"{0.05 * n_genes:.0f} = 5% if p were UNCORRECTED (not this test)")
    ax.set(xlabel="random donor split (control cells only)",
           ylabel=f"genes called significant\n(of {n_genes:,}, BH q < 0.05)",
           title="Real data (Kang 2018, CD14+ monocytes), donors shuffled:\n"
                 f"per-cell median {df['per_cell'].median():.0f} false genes per "
                 f"split; pseudobulk median {df['pseudobulk'].median():.0f} "
                 f"(max {df['pseudobulk'].max():.0f})")
    ax.title.set_fontsize(10)
    ax.set_xticks(df["split"] + 1)
    # Upper left: the per-cell series climbs left-to-right, so the top-right
    # corner is exactly where its final point sits.
    ax.legend(frameon=False, fontsize=8, loc="upper center",
              bbox_to_anchor=(0.5, -0.18), ncol=2)
    share = df["per_cell"] / df["n_genes_sc"]
    _save(fig, "real_null.png",
          f"n = {len(df)} random donor splits of control cells only "
          f"({int(df['n_cells'].iloc[0]):,} cells). Per-cell: median "
          f"{df['per_cell'].median():.0f} ({share.median():.2%} of genes), mean "
          f"{df['per_cell'].mean():.1f} ({share.mean():.1%}), max {df['per_cell'].max()}; "
          f"pseudobulk max {df['pseudobulk'].max()}.\n"
          "source: results/false_positives.csv")


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
    for off, col in ((-0.2, "median_isg_rank_by_q"),
                     (0.2, "median_isg_rank_by_effect")):
        for xi, v in zip(x + off, df[col]):
            axes[0].text(xi, v * 1.12, f"{v:.0f}", ha="center", va="bottom",
                         fontsize=7)
    axes[0].set_yscale("log")
    axes[0].set_ylim(1, df["median_isg_rank_by_q"].max() * 2.5)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    axes[0].set_ylabel("median rank of known interferon genes\n"
                       "(among genes significant at BH q < 0.05)")
    axes[0].set_title(
        "Ranked by q-value, interferon genes sit at median rank "
        f"{df['median_isg_rank_by_q'].median():.0f};\n"
        f"by effect size {df['median_isg_rank_by_effect'].median():.0f} "
        f"(median over {len(df)} cell types)", fontsize=10)
    axes[0].legend(frameon=False, fontsize=8, loc="upper right")

    axes[1].bar(x - 0.2, df["isgs_top12_by_q"], width=0.4,
                label="ranked by q-value", color="C3",
                edgecolor="black", linewidth=0.5)
    axes[1].bar(x + 0.2, df["isgs_top12_by_effect"], width=0.4,
                label="ranked by effect size", color="C0",
                edgecolor="black", linewidth=0.5)
    for off, col in ((-0.2, "isgs_top12_by_q"), (0.2, "isgs_top12_by_effect")):
        for xi, v in zip(x + off, df[col]):
            axes[1].text(xi, v + 0.15, f"{v}", ha="center", va="bottom", fontsize=7.5)
    axes[1].set_ylim(0, 13)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    # One line, not two: the wrapped version pushed its last word against the
    # top edge of the canvas.
    axes[1].set_ylabel("known interferon genes in top 12")
    axes[1].set_title(
        "Interferon genes among the top 12: "
        f"{df['isgs_top12_by_q'].sum()} of {12 * len(df)} by q-value,\n"
        f"{df['isgs_top12_by_effect'].sum()} of {12 * len(df)} by effect size "
        f"(12 genes x {len(df)} cell types)", fontsize=10)
    axes[1].legend(frameon=False, fontsize=8, loc="upper right")

    _save(fig, "isg_ranking.png",
          f"n = {len(df)} cell types (5-8 donors each), 12 canonical interferon genes "
          "per type, paired pseudobulk test. Ranks are within the genes already "
          "significant at BH q < 0.05 (README 'Which rank').\n"
          "source: results/real_analysis_summary.csv")


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
           # Wording (2026-11-15): "make the same null calls" claimed the same
           # GENES are called; the artefacts store per-draw counts, not call
           # sets, so only the count is measured. Say what was measured.
           title="Mixed model and pseudobulk make the same number of null\n"
                 f"false calls: identical in {same} of {n_draws} balanced draws")
    ax.title.set_fontsize(10)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    _save(fig, "mixed_vs_pseudobulk_null.png",
          f"n = {n_draws} balanced donor-swap draws per method ({by_type} by type), "
          f"{n_genes} small-effect genes each; {any_fp} draws had any false "
          "positive. Bars: mean; whiskers: bootstrap 95% CI. A swap null can show "
          "anti-conservatism, not calibration.\n"
          "source: results/mixed_null_donor_swap.csv, results/mixed_null_draws.csv")


def fig_nbglm_null():
    """The NB GLM makes false calls on a balanced null in some types; pseudobulk does not."""
    df = _maybe("nbglm_real_null_balanced.csv", "scripts/nbglm_real.py")
    if df is None:
        return
    bh = _maybe("bh_null_expectation.csv", "scripts/bh_null_expectation.py")
    pb = df[df["method"] == "pseudobulk"].set_index("cell_type")
    nb = df[df["method"] == "nb_glm"].set_index("cell_type").loc[pb.index]
    ref = 0.05
    # highlighted: types whose Wilson lower bound clears the reference
    hot = [ct for ct in nb.index if nb.at[ct, "p_any_lo"] > ref]
    hot.sort(key=lambda ct: -nb.at[ct, "p_any"])

    def frac(tab, ct):
        return f"{int(tab.at[ct, 'n_draws_with_any_fp'])}/{int(tab.at[ct, 'draws'])}"

    fig, ax = plt.subplots(figsize=(6.5, 4))
    x = np.arange(len(pb))
    for ct in hot:
        i = list(pb.index).index(ct)
        ax.axvspan(i - 0.5, i + 0.5, color="C3", alpha=0.08, lw=0)
    for off, tab, color, label in ((-0.2, pb, "C0", "pseudobulk paired t"),
                                   (0.2, nb, "C3", "NB GLM")):
        p = tab["p_any"].to_numpy()
        lo = tab["p_any_lo"].to_numpy()
        hi = tab["p_any_hi"].to_numpy()
        ax.bar(x + off, p, width=0.4, label=label, color=color,
               edgecolor="black", linewidth=0.5, yerr=[p - lo, hi - p],
               capsize=2, error_kw={"lw": 0.8})
        for xi, ct, h in zip(x + off, tab.index, hi):
            bold = tab is nb and ct in hot
            ax.text(xi, h + 0.02, frac(tab, ct), ha="center", va="bottom",
                    fontsize=7.5 if bold else 6.5,
                    fontweight="bold" if bold else "normal")
    ax.axhline(ref, ls="--", lw=1, color="grey",
               label=f"{ref:.2f}: reference for an exact null")
    mean_lines = [f"{_short(ct)}: {nb.at[ct, 'mean_fp_small_effect']:.2f} "
                  f"({nb.at[ct, 'fp_ci_lo']:.2f}-{nb.at[ct, 'fp_ci_hi']:.2f})"
                  for ct in hot]
    if mean_lines:
        ax.text(0.01, 0.99, "NB GLM false calls per draw (95% CI)\n"
                + "\n".join(mean_lines) + "\npseudobulk: "
                + ", ".join(f"{pb.at[ct, 'mean_fp_small_effect']:.2f}" for ct in hot),
                transform=ax.transAxes, va="top", fontsize=7)
    ax.set_ylim(0, 1.62)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticks(x)
    ax.set_xticklabels([f"{_tick(ct)}\n{int(d)} donors\n{int(n)} draws"
                        for ct, d, n in zip(pb.index, pb["n_donors_used"],
                                            pb["draws"])], fontsize=7.5)
    parts = ", ".join(f"{frac(nb, ct)} {_short(ct)}" for ct in hot)
    if all(pb.at[ct, "n_draws_with_any_fp"] == 0 for ct in hot):
        pb_parts = "none"
    else:
        pb_parts = ", ".join(frac(pb, ct) for ct in hot)
    ax.set(ylabel="share of balanced draws with at least\none false positive (BH q < 0.05)",
           title="The NB GLM is anti-conservative on a balanced null:\n"
                 f"false calls in {parts} draws; pseudobulk: {pb_parts}")
    ax.title.set_fontsize(10)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    enumerated = bool((nb["n_distinct_draws"] == nb["n_plans_available"]).all())
    by_type = " + ".join(str(int(v)) for v in nb["draws"])
    ref_note = ""
    if bh is not None:
        ref_note = (f" The {ref:.2f} line is P(any rejection) for independent "
                    f"exact nulls (simulated {bh['p_any'].iloc[0]:.4f}); it is not "
                    "a bound for correlated genes.")
    _save(fig, "nbglm_null_anticonservative.png",
          f"n = {by_type} distinct balanced donor-swap draws "
          f"({'every plan enumerated' if enumerated else 'sampled plans'}), "
          "small-effect genes only (|log2FC| < 0.5, ISGs removed). Whiskers: "
          "Wilson 95% CI. Shaded: Wilson lower bound above the reference."
          + ref_note + "\n"
          "source: results/nbglm_real_null_balanced.csv, results/bh_null_expectation.csv")


def _model_slug(cell_type):
    """File stem used by scripts/mixed_vs_others.py ('CD14+ Monocytes' -> 'CD14_Monocytes')."""
    return (cell_type.replace("+", "").replace(" ", "_")
            .replace("-", "_").replace("__", "_"))


def fig_boundary():
    """Boundary fits carry most of the mixed model's extra calls."""
    summ = _maybe("mixed_vs_others_summary.csv", "scripts/mixed_vs_others.py")
    if summ is None:
        return
    parts, rows = [], []
    for _, s in summ.iterrows():
        path = RESULTS / "mixed_per_type" / f"{_model_slug(s['cell_type'])}_models.csv"
        if not path.exists():
            print(f"  (skipping: {path.name} is missing; run scripts/mixed_vs_others.py)")
            return
        m = pd.read_csv(path, index_col=0)
        extra = (m["mx_qval"] < 0.05) & ~(m["pb_qval"] < 0.05)
        row = {"cell_type": s["cell_type"], "genes": len(m),
               "boundary": int(m["mx_boundary"].sum()),
               "extra": int(extra.sum()),
               "extra_boundary": int((extra & m["mx_boundary"]).sum()),
               "pb_only": int(((m["pb_qval"] < 0.05) & ~(m["mx_qval"] < 0.05)).sum())}
        # the per-gene files and the summary must tell the same story
        if (row["extra"], row["extra_boundary"], row["pb_only"]) != (
                int(s["mixed_only_calls"]), int(s["mixed_only_calls_boundary"]),
                int(s["pseudobulk_only_calls"])):
            raise ValueError(f"{s['cell_type']}: models csv and summary disagree on "
                             "mixed-only calls")
        rows.append(row)
        m = m.assign(extra=extra)
        parts.append(m[["paired_t_stat", "mx_stat", "mx_boundary", "extra"]])
    tab = pd.DataFrame(rows).set_index("cell_type")
    allg = pd.concat(parts)
    tt = allg["paired_t_stat"].abs().clip(lower=1e-3)
    zz = allg["mx_stat"].abs().clip(lower=1e-3)
    bnd = allg["mx_boundary"].to_numpy()
    ratio = zz / tt
    r_bnd, r_int = ratio[bnd].max(), ratio[~bnd].max()
    n_extra, n_extra_b = int(tab["extra"].sum()), int(tab["extra_boundary"].sum())
    share_extra = n_extra_b / n_extra
    share_genes = tab["boundary"].sum() / tab["genes"].sum()

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3),
                             gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    ax.scatter(tt[~bnd], ratio[~bnd], s=3, color="C0", alpha=0.35, lw=0,
               rasterized=True, label=f"interior: donor variance estimated "
                                      f"(n = {int((~bnd).sum()):,})")
    ax.scatter(tt[bnd], ratio[bnd], s=3, color="C3", alpha=0.35, lw=0,
               rasterized=True, label=f"boundary: donor variance at 0 "
                                      f"(n = {int(bnd.sum()):,})")
    ex = allg["extra"].to_numpy()
    ax.scatter(tt[ex], ratio[ex], s=14, facecolors="none", edgecolors="black",
               linewidths=0.6, label=f"mixed-only call (n = {int(ex.sum())})")
    ax.axhline(1, ls="--", lw=1, color="grey", label="mixed = paired t")
    ax.set(xscale="log", xlim=(tt.min() * 0.8, tt.max() * 1.25),
           ylim=(0.97, r_bnd * 1.1),
           xlabel="pseudobulk paired t statistic (|t|, log scale)",
           ylabel="mixed Wald statistic / paired t",
           title="Boundary fits lift the mixed statistic up to "
                 f"{r_bnd:.2f}x above the paired t;\n"
                 f"interior fits agree within {100 * (r_int - 1):.2f}%")
    ax.title.set_fontsize(10)
    # Lower right: the point cloud is a diagonal band rising to the upper right,
    # so upper left is where the dense boundary points are, not empty space.
    # Outside the axes: the cloud is a diagonal band and every in-axes corner
    # holds points, so no rectangle can be guaranteed empty.
    ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.16),
              markerscale=2.5, frameon=False, ncol=2)

    ax = axes[1]
    x = np.arange(len(tab))
    interior = tab["extra"] - tab["extra_boundary"]
    ax.bar(x, tab["extra_boundary"], width=0.6, color="C3", edgecolor="black",
           linewidth=0.5, label="boundary fit")
    ax.bar(x, interior, width=0.6, bottom=tab["extra_boundary"], color="C0",
           edgecolor="black", linewidth=0.5, label="interior fit")
    for xi, b, e in zip(x, tab["extra_boundary"], tab["extra"]):
        ax.text(xi, e + 2, f"{b}/{e}", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(x)
    ax.set_xticklabels([_short(ct) for ct in tab.index], rotation=35,
                       ha="right", fontsize=8)
    ax.set_ylim(0, tab["extra"].max() * 1.18)
    ax.set(ylabel="genes called by mixed but not pseudobulk\n(BH q < 0.05)",
           title=f"{share_extra:.0%} of mixed-only calls ({n_extra_b} of "
                 f"{n_extra}) are boundary fits;\nboundary genes are "
                 f"{share_genes:.0%} of all genes tested")
    ax.title.set_fontsize(10)
    ax.text(0.98, 0.6, f"pseudobulk-only calls: {int(tab['pb_only'].sum())}",
            transform=ax.transAxes, ha="right", fontsize=8)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    _save(fig, "mixed_boundary_genes.png",
          f"n = {int(tab['genes'].sum()):,} gene fits over {len(tab)} cell types "
          f"(left); {n_extra} mixed-only calls (right: boundary/total per type). "
          "Both panels pool every tested gene. "
          "Boundary = fitted donor variance at 0. Mixed-only = BH q < 0.05 under "
          "the mixed model and not under the paired t.\n"
          "source: results/mixed_per_type/*_models.csv, "
          "results/mixed_vs_others_summary.csv")


def fig_plan_space():
    """Why '200 draws' cannot exist: the balanced plan space is tiny."""
    swap = _maybe("mixed_null_donor_swap.csv", "scripts/mixed_vs_others.py")
    if swap is None:
        return
    bal = swap[(swap["balance"] == "balanced") & (swap["method"] == "mixed")]
    cap = None
    driver = ROOT / "scripts" / "mixed_vs_others.py"
    if driver.exists():
        found = re.search(r"^NULL_PERMS\s*=\s*(\d+)", driver.read_text(), re.M)
        cap = int(found.group(1)) if found else None

    groups = []
    for used, g in sorted(bal.groupby("n_donors_used"), reverse=True):
        plans = g["n_plans_available"].unique()
        draws = g["n_distinct_draws"].unique()
        paired = g["n_paired_donors"].unique()
        if len(plans) != 1 or len(draws) != 1 or len(paired) != 1:
            raise ValueError(f"{used} donors used: types disagree on the plan space")
        used, plans, draws, paired = int(used), int(plans[0]), int(draws[0]), int(paired[0])
        holdout = math.comb(paired, used)          # which donors are held out
        if holdout * math.comb(used, used // 2) // 2 == plans:
            how = (f"{holdout} x " if holdout > 1 else "") + f"1/2 C({used},{used // 2})"
        else:
            how = ""                               # do not explain what we cannot reproduce
        groups.append({"used": used, "paired": paired, "plans": plans,
                       "draws": draws, "how": how, "types": list(g["cell_type"])})

    fig, ax = plt.subplots(figsize=(5.5, 4))
    x = np.arange(len(groups))
    plans = [g["plans"] for g in groups]
    ax.bar(x, plans, width=0.55, color="C0", edgecolor="black", linewidth=0.5,
           label="distinct balanced plans = draws run")
    for xi, g in zip(x, groups):
        ax.text(xi, g["plans"] + 4, f"= {g['how']}" if g["how"] else "",
                ha="center", va="bottom", fontsize=7.5)
        ax.text(xi, g["plans"] + 16, f"{g['plans']}", ha="center", va="bottom",
                fontsize=10, fontweight="bold")
    top = max(plans)
    if cap is not None:
        # Labelled at the line's right end instead of in a legend: any legend in
        # the upper area is drawn across this full-width line, and the line is the
        # point of the figure.
        ax.axhline(cap, ls="--", lw=1, color="grey")
        ax.annotate(f"{cap} draws requested (NULL_PERMS cap)", xy=(1.0, cap),
                    xycoords=("axes fraction", "data"), xytext=(-2, -3),
                    textcoords="offset points", ha="right", va="top",
                    fontsize=7.5, color="0.3",
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.5))
        top = max(top, cap)
    ax.set_ylim(0, top * 1.12)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g['used']} donors used\n(of {g['paired']} paired)\n"
                        + "\n".join(_short(t) for t in g["types"])
                        for g in groups], fontsize=7.5)
    ax.set(ylabel="distinct balanced donor-swap plans",
           title="The plan space is tiny: "
                 + ", ".join(str(p) for p in plans[:-1]) + f" and {plans[-1]} plans at "
                 + ", ".join(str(g["used"]) for g in groups[:-1])
                 + f" and {groups[-1]['used']} donors,\n"
                 + (f"so {cap} independent draws cannot exist" if cap else
                    "so many independent draws cannot exist"))
    ax.title.set_fontsize(10)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    enumerated = all(g["plans"] == g["draws"] for g in groups)
    _save(fig, "null_plan_space.png",
          f"n = {len(bal)} cell types; every plan was "
          f"{'enumerated, so draws run = plans that exist' if enumerated else 'not enumerated'}. "
          "A plan and its complement give the same |t|; with an odd donor count "
          "one donor is held out per draw.\n"
          "source: results/mixed_null_donor_swap.csv (n_plans_available, "
          "n_distinct_draws); cap: NULL_PERMS in scripts/mixed_vs_others.py")


def fig_log2fc_agreement():
    """One normalisation: the GLM and pseudobulk log2FC columns now coincide."""
    after = _maybe("nbglm_real_summary.csv", "scripts/nbglm_real.py")
    if after is None:
        return
    before = _maybe("nbglm_real_log2fc_before_fix.csv",
                    "(committed copy of commit 1c0c940)", comment="#")
    isg = _maybe("nbglm_real_isg.csv", "scripts/nbglm_real.py")
    types = list(after["cell_type"])
    x = np.arange(len(types))
    a_diff = after["mean_abs_diff"].to_numpy()
    floor = 1e-17

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3),
                             gridspec_kw={"width_ratios": [1.1, 1]})
    ax = axes[0]
    # dots, not bars: a bar on a log axis has no meaningful length
    shown = np.maximum(a_diff, 2 * floor)      # an exact 0 cannot sit on a log axis
    if before is not None:
        before = before.set_index("cell_type").loc[types]
        b_diff = before["mean_abs_diff"].to_numpy()
        for xi, b, a in zip(x, b_diff, shown):
            ax.plot([xi, xi], [a, b], color="grey", lw=0.8, zorder=0)
        ax.plot(x, b_diff, "o", color="C3", mec="black", mew=0.5,
                label="before: GLM on total-count CPM, pseudobulk on size factors")
        for xi, v in zip(x, b_diff):
            # Opaque (not alpha 0.85) and above the stems: the grey before/after
            # connectors are drawn at zorder 0 but a translucent mask still let
            # them show through the digits.
            ax.text(xi, v * 1.6, f"{v:.2f}", ha="center", va="bottom", fontsize=7,
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.9),
                    zorder=6)
    ax.plot(x, shown, "o", color="C0", mec="black", mew=0.5,
            label="now: both on median-of-ratios size factors")
    for xi, v, s in zip(x, a_diff, shown):
        ax.text(xi, s * 1.6, "0 (exact)" if v == 0 else f"{v:.0e}", ha="center",
                va="bottom", fontsize=7,
                bbox=dict(facecolor="white", edgecolor="none", pad=0.9), zorder=6)
    ax.set_yscale("log")
    # Bottom at floor/4, not floor: the exact-zero marker is drawn at 2*floor and
    # its label sat on the axis spine, clipped by the frame.
    # Top set from the data, not a constant: the highest red label sits at
    # max(b_diff)*1.6, so a fixed limit of 12 left "0.45" pressed against the frame
    # when the largest before-fix difference was small.
    ax.set_ylim(floor / 4, max(float(b_diff.max()) * 3.0, 5.0))
    ax.set_xlim(-0.6, len(types) - 0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([_tick(t) for t in types], fontsize=7.5)
    title = ("One shared normalisation: the two log2FC columns coincide\n"
             f"(mean |diff| at most {a_diff.max():.0e} log2)")
    if before is not None:
        title = ("One shared normalisation: the two log2FC columns coincide\n"
                 f"(mean |diff| {b_diff.min():.2f}-{b_diff.max():.2f} before, "
                 f"at most {a_diff.max():.0e} now)")
    ax.set(ylabel="mean |GLM log2FC - pseudobulk log2FC|\n(log2 units, log scale)",
           title=title)
    ax.title.set_fontsize(10)
    # Below the axes rather than in the middle: every x has a vertical stem, so
    # any in-axes position overlaps one.
    ax.legend(frameon=False, fontsize=7.5, loc="upper center",
              bbox_to_anchor=(0.5, -0.16), ncol=1)

    ax = axes[1]
    sep = after["glm_separated"].to_numpy()
    ax.bar(x, sep, width=0.6, color="C3", edgecolor="black", linewidth=0.5)
    for xi, v in zip(x, sep):
        ax.text(xi, v + 0.15, str(int(v)), ha="center", va="bottom", fontsize=8)
    n_sep, n_called = int(sep.sum()), int(after["glm_called_separated"].sum())
    n_failed = int(after["glm_failed"].sum())
    flagged = []
    if isg is not None:
        for _, row in after.iterrows():
            genes = str(row["glm_separated_genes"]).split(";")
            panel = isg[isg["cell_type"] == row["cell_type"]].set_index("gene")
            for g in genes:
                if g in panel.index:
                    q = panel.at[g, "q_pb"]
                    flagged.append((row["cell_type"], g, q, panel.at[g, "q_nb"]))
    ax.set_xticks(x)
    ax.set_xticklabels([f"{_tick(t)}\nr = {r:.2f}" for t, r in
                        zip(types, after["glm_coef_pearson_ok_fits"])],
                       fontsize=7.5)
    ax.set_ylim(0, max(sep.max(), 1) * 1.5)
    line2 = (f"{len(flagged)} is a panel interferon gene" if len(flagged) == 1
             else f"{len(flagged)} are panel interferon genes" if flagged
             else "none is a panel interferon gene")
    ax.set(ylabel="genes with a separated GLM fit",
           title=f"{n_sep} GLM fits are separated: {n_called or 'none'} called, "
                 f"{n_failed or 'none'} failed outright;\n{line2}")
    ax.title.set_fontsize(10)
    r_ok = after["glm_coef_pearson_ok_fits"]
    ax.text(0.02, 0.98, "r under each type: Pearson of the GLM's own\n"
            "coefficient vs the pseudobulk log2FC on the\n"
            f"fits that are estimates ({r_ok.min():.2f}-{r_ok.max():.2f})",
            transform=ax.transAxes, va="top", fontsize=7)
    for ct, g, q_pb, q_nb in flagged:
        xi = types.index(ct)
        ax.annotate(f"{g} ({_short(ct)}):\npseudobulk q {q_pb:.0e}, GLM q = {q_nb:.0f}",
                    xy=(xi, after.set_index("cell_type").at[ct, "glm_separated"]),
                    # Higher and masked: the leader line used to terminate on the
                    # bar's own "2" value label, which then read as struck through.
                    xytext=(xi + 0.4, max(sep) * 1.18), fontsize=7.5, ha="left",
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.8), zorder=6,
                    arrowprops={"arrowstyle": "-", "color": "grey", "lw": 0.8})
    pr = ("Pearson r between the two log2FC columns: "
          + (f"{before['pearson'].min():.3f}-{before['pearson'].max():.3f} before, "
             if before is not None else "")
          + f"{after['pearson'].min():.3f} now (identity under one normalisation, not "
          "independent agreement). Right: r = Pearson of the GLM's own coefficient vs "
          "the pseudobulk log2FC on the fits that are estimates.")
    _save(fig, "nbglm_log2fc_agreement.png",
          f"n = {len(types)} cell types ({int(after['donors'].min())}-"
          f"{int(after['donors'].max())} paired donors); {int(after['n_genes'].sum()):,} "
          f"gene fits. {pr}\n"
          "source: results/nbglm_real_summary.csv, "
          "results/nbglm_real_log2fc_before_fix.csv (verbatim columns of the same "
          "file at commit 1c0c940), results/nbglm_real_isg.csv")


def fig_megakaryocyte():
    """Megakaryocytes have no usable expression design: one donor supplies 84%."""
    counts = _maybe("megakaryocyte_donor_counts.csv", "scripts/megakaryocyte_check.py")
    sweep = _maybe("megakaryocyte_threshold_sweep.csv", "scripts/megakaryocyte_check.py")
    summ = _maybe("megakaryocyte_summary.csv", "scripts/megakaryocyte_check.py")
    if counts is None or sweep is None or summ is None:
        return
    s = summ.iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2),
                             gridspec_kw={"width_ratios": [1.1, 1]})
    ax = axes[0]
    d = counts.sort_values("total").reset_index(drop=True)
    y = np.arange(len(d))
    ax.barh(y, d["ctrl"], color="0.65", edgecolor="black", linewidth=0.5,
            label="control cells")
    ax.barh(y, d["stim"], left=d["ctrl"], color="C2", edgecolor="black",
            linewidth=0.5, label="stimulated cells")
    for yi, c, st, tot, share in zip(y, d["ctrl"], d["stim"], d["total"],
                                     d["share_of_type"]):
        ax.text(tot + 1.5, yi, f"{tot} = {c} + {st}"
                + (f" ({share:.0%})" if share > 0.5 else ""), va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(d["donor"].str.replace("_", " "), fontsize=8)
    ax.set_xlim(0, d["total"].max() * 1.5)
    ax.set(xlabel="megakaryocytes (cells)",
           title=f"One donor supplies {s['largest_donor_share']:.0%} of megakaryocytes\n"
                 f"({int(d['total'].max())} of {int(s['n_cells'])} cells, "
                 f"{len(d)} donors)")
    ax.title.set_fontsize(10)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_xlabel("megakaryocytes (cells; total = control + stimulated)")

    ax = axes[1]
    xs = np.arange(len(sweep))
    ok = sweep["test_possible"].to_numpy()
    n_don = sweep["donors_qualifying"].to_numpy()
    ax.plot(xs, n_don, "-", color="grey", lw=1, zorder=1)
    ax.scatter(xs[ok], n_don[ok], s=55, color="C0", edgecolor="black",
               linewidth=0.5, zorder=3, label="2 or more donors (test can run)")
    ax.scatter(xs[~ok], n_don[~ok], s=55, color="C3", edgecolor="black",
               linewidth=0.5, zorder=3, label="1 donor (no paired test)")
    ax.set_xticks(xs)
    ax.set_xticklabels(sweep["min_cells_per_condition"])
    ax.set_yticks(range(0, int(n_don.max()) + 1))
    ax.set_ylim(0, n_don.max() + 1.6)
    first_single = int(sweep.loc[~sweep["test_possible"], "min_cells_per_condition"].iloc[0])
    ax.set(xlabel="minimum cells per donor and condition",
           ylabel="donors that qualify",
           title=f"No threshold gives a usable design: {n_don[0]} donors at "
                 f">= {sweep['min_cells_per_condition'].iloc[0]} cell,\n"
                 f"{n_don[1]} at >= {sweep['min_cells_per_condition'].iloc[1]}, "
                 f"one donor from >= {first_single} cells up")
    ax.title.set_fontsize(10)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    # The design note lives in the footer, not inside the axes: at fontsize 7.5
    # it was wide enough to run through the grey connecting line and the markers.
    # The footer is wrapped to the figure width by _save(), so it cannot overflow.
    design_note = (
        f"loosest usable design: {int(s['donors_in_loosest_design'])} donors, "
        f"profiles of {int(s['cells_per_profile_min'])}-"
        f"{int(s['cells_per_profile_max'])} cells; "
        f"{int(s['genes_tested'])} genes tested, {int(s['genes_significant'])} called; "
        f"minimum detectable effect (80% power) "
        f"{s['mde_log2_exact_noncentral_t']:.2f} log2 "
        f"({s['mde_log2_exact_bonferroni']:.2f} with Bonferroni) vs largest observed "
        f"|log2FC| {s['largest_observed_abs_log2fc']:.2f}. ")
    _save(fig, "megakaryocyte_donors.png",
          design_note +
          f"n = {int(s['n_cells'])} megakaryocytes ({int(s['ctrl_cells'])} control, "
          f"{int(s['stim_cells'])} stimulated) across {len(d)} donors; the composition "
          "test (all 8 donors) is unaffected. The sweep counts donors with at least "
          "the given number of cells in both conditions.\n"
          "source: results/megakaryocyte_donor_counts.csv, "
          "results/megakaryocyte_threshold_sweep.csv, results/megakaryocyte_summary.csv")


def fig_nbglm_sim():
    """Simulated: the NB GLM's power advantage is confined to small donor counts."""
    s = _maybe("sim_nbglm_summary.csv", "scripts/nbglm_study.py")
    raw = _maybe("sim_nbglm.csv", "scripts/nbglm_study.py")
    if s is None or raw is None:
        return
    names = {"nb_glm": "NB GLM", "paired_t_logcpm": "paired t on log-CPM"}
    cols = {"nb_glm": "C3", "paired_t_logcpm": "C0"}
    donors = sorted(s["donors"].unique())
    x = np.arange(len(donors))
    reps = int(raw.groupby(["method", "donors"])["rep"].nunique().min())
    n_tested = int(raw["n_tested"].iloc[0])
    n_true = int((raw["tp"] + raw["fn"]).iloc[0])

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, col, ylabel in ((axes[0], "power", "power (share of changed genes found)"),
                            (axes[1], "fdp", "observed false discovery proportion")):
        for off, m in ((-0.06, "nb_glm"), (0.06, "paired_t_logcpm")):
            g = s[s["method"] == m].set_index("donors").loc[donors]
            # whiskers on power only: FDP is 0 whenever a replicate calls nothing,
            # so its SD at 3 donors says nothing useful
            ax.errorbar(x + off, g[f"{col}_mean"],
                        yerr=g[f"{col}_sd"] if col == "power" else None, fmt="o-",
                        color=cols[m], mec="black", mew=0.5, capsize=2, lw=1.2,
                        elinewidth=0.8, label=names[m])
        ax.set_xticks(x)
        ax.set_xticklabels(donors)
        ax.set_xlabel("donors (paired design)")
        ax.set_ylabel(ylabel)
    g6 = s[s["donors"] == 6].set_index("method")
    g12 = s[s["donors"] == 12].set_index("method")
    for off, m in ((-0.06, "nb_glm"), (0.06, "paired_t_logcpm")):
        axes[0].annotate(f"{g6.at[m, 'power_mean']:.0%}",
                         (list(donors).index(6) + off, g6.at[m, "power_mean"]),
                         textcoords="offset points",
                         xytext=(-30 if m == "nb_glm" else 10, 4 if m == "nb_glm" else -4),
                         fontsize=9, fontweight="bold")
    axes[0].set_ylim(-0.05, 1.12)
    axes[0].set_title(f"Simulated: at 6 donors the NB GLM finds {g6.at['nb_glm', 'power_mean']:.0%}\n"
                      f"of changed genes, the paired t {g6.at['paired_t_logcpm', 'power_mean']:.0%}; "
                      "equal from 12", fontsize=10)
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")
    axes[1].axhline(0.05, ls="--", lw=1, color="grey", label="nominal 5%")
    axes[1].set_ylim(-0.005, 0.14)
    axes[1].set_title("False discovery proportion: the GLM is lower at 6 donors\n"
                      f"({g6.at['nb_glm', 'fdp_mean']:.3f} vs "
                      f"{g6.at['paired_t_logcpm', 'fdp_mean']:.3f}) but exceeds 5% at 12 "
                      f"({g12.at['nb_glm', 'fdp_mean']:.3f})", fontsize=10)
    axes[1].legend(frameon=False, fontsize=8, loc="upper right")
    _save(fig, "nbglm_sim_power.png",
          f"n = {reps} simulated replicates per method and donor count (mean; power "
          f"whiskers +/-1 SD, FDP means only); {n_tested} genes, {n_true} truly changed at log2 effect 1.0, BH "
          "alpha 0.05 (scripts/nbglm_study.py). SIMULATED: the generator draws negative-"
          "binomial counts, the GLM's own model, and the dispersion-trend floor was "
          "tuned on it, so this is the GLM's best case and not out-of-sample.\n"
          "source: results/sim_nbglm_summary.csv, results/sim_nbglm.csv")


def _panel_letter(ax, letter):
    ax.text(-0.02, 1.16, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="bottom", ha="right")


def fig_at_a_glance():
    """Four panels that summarise the repository."""
    sim = _maybe("sim_null.csv", "scripts/simulation_study.py")
    real = _maybe("false_positives.csv", "scripts/false_positive_test.py")
    power = _maybe("sim_power.csv", "scripts/simulation_study.py")
    summ = _maybe("mixed_vs_others_summary.csv", "scripts/mixed_vs_others.py")
    draws = _maybe("mixed_null_draws.csv", "scripts/mixed_vs_others.py")
    nbnull = _maybe("nbglm_real_null_balanced.csv", "scripts/nbglm_real.py")
    nbsum = _maybe("nbglm_real_summary.csv", "scripts/nbglm_real.py")
    if any(t is None for t in (sim, real, power, summ, draws, nbnull, nbsum)):
        return
    fig, axes = plt.subplots(2, 2, figsize=(9, 7.4))

    # A: false positives on null data, simulated and real
    ax = axes[0, 0]
    sim_pc = 100 * sim["per_cell"] / sim["n_genes"]
    sim_pb = 100 * sim["pseudobulk"] / sim["n_genes"]
    real_pc = 100 * real["per_cell"] / real["n_genes_sc"]
    real_pb = 100 * real["pseudobulk"] / real["n_genes_pb"]
    groups = [(f"simulated\n{int(sim['n_genes'].iloc[0])} genes,\n{len(sim)} replicates",
               sim_pc, sim_pb),
              (f"real Kang CD14+ Mono\n{int(real['n_genes_sc'].iloc[0]):,} genes,\n"
               f"{len(real)} donor splits", real_pc, real_pb)]
    x = np.arange(len(groups))
    for off, k, color, label in ((-0.2, 1, "C3", "per-cell Wilcoxon"),
                                 (0.2, 2, "C0", "pseudobulk paired t")):
        vals = [g[k] for g in groups]
        mean = np.array([v.mean() for v in vals])
        lo = np.array([v.min() for v in vals])
        hi = np.array([v.max() for v in vals])
        ax.bar(x + off, mean, width=0.4, color=color, edgecolor="black",
               linewidth=0.5, label=label, yerr=[mean - lo, hi - mean],
               capsize=2, error_kw={"lw": 0.8})
        for xi, m, h in zip(x + off, mean, hi):
            ax.text(xi, h + 1.5, f"{m:.1f}%", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([g[0] for g in groups], fontsize=8)
    ax.set_ylim(0, 100)
    ax.set(ylabel="genes called significant\n(% of genes tested, BH q < 0.05)",
           title=f"Null data: per-cell testing calls {sim_pc.mean():.0f}% (simulated)\n"
                 f"and {real_pc.mean():.0f}% (real) of genes; pseudobulk "
                 f"{max(sim_pb.mean(), real_pb.mean()):.1f}%")
    ax.title.set_fontsize(9.5)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    _panel_letter(ax, "A")

    # B: power and false discovery proportion
    ax = axes[0, 1]
    agg = power.groupby("method")[["power", "fdp"]].mean()
    order = ["per_cell", "pseudobulk"]
    names = {"per_cell": "per-cell\nWilcoxon", "pseudobulk": "pseudobulk\npaired t"}
    cols = {"per_cell": "C3", "pseudobulk": "C0"}
    for xi, m in enumerate(order):
        ax.bar(xi, agg.at[m, "fdp"], width=0.5, color=cols[m], edgecolor="black",
               linewidth=0.5)
        ax.text(xi, agg.at[m, "fdp"] + 0.02, f"{agg.at[m, 'fdp']:.2f}", ha="center",
                va="bottom", fontsize=9, fontweight="bold")
    ax.axhline(0.05, ls="--", lw=1, color="grey", label="nominal 5%")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{names[m]}\npower {agg.at[m, 'power']:.2f}" for m in order],
                       fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set(ylabel="observed false discovery proportion",
           title=f"Same power ({agg.at['per_cell', 'power']:.2f} vs "
                 f"{agg.at['pseudobulk', 'power']:.2f}), but per-cell hit lists are\n"
                 f"{agg.at['per_cell', 'fdp']:.0%} false against "
                 f"{agg.at['pseudobulk', 'fdp']:.0%} for pseudobulk")
    ax.title.set_fontsize(9.5)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    _panel_letter(ax, "B")

    # C: mixed model vs pseudobulk
    ax = axes[1, 0]
    ratio = summ["mixed_sig"] / summ["pseudobulk_sig"]
    x = np.arange(len(summ))
    ax.bar(x, ratio, width=0.6, color="C1", edgecolor="black", linewidth=0.5,
           label="mixed (1 | donor)")
    for xi, v in zip(x, ratio):
        ax.text(xi, v + 0.02, f"{v:.2f}x", ha="center", va="bottom", fontsize=8)
    ax.axhline(1, ls="--", lw=1, color="grey", label="same calls as pseudobulk")
    ax.set_xticks(x)
    ax.set_xticklabels([_short(t) for t in summ["cell_type"]], rotation=35,
                       ha="right", fontsize=8)
    ax.set_ylim(0, 1.45)
    bd = draws[draws["balance"] == "balanced"]
    same = int((bd["fp_pseudobulk"] == bd["fp_mixed"]).sum())
    share = summ["mixed_only_calls_boundary"].sum() / summ["mixed_only_calls"].sum()
    ax.set(ylabel="genes called, mixed / pseudobulk\n(BH q < 0.05, real data)",
           title=f"Mixed model: {ratio.min():.2f}-{ratio.max():.2f}x pseudobulk's calls, "
                 f"{share:.0%} of the extras\non boundary fits; same null false-positive "
                 f"count in {same}/{len(bd)} draws")
    ax.title.set_fontsize(9.5)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    _panel_letter(ax, "C")

    # D: NB GLM calls more genes and is anti-conservative on two types
    ax = axes[1, 1]
    pb = nbnull[nbnull["method"] == "pseudobulk"].set_index("cell_type")
    nb = nbnull[nbnull["method"] == "nb_glm"].set_index("cell_type").loc[pb.index]
    x = np.arange(len(pb))
    for off, tab, color, label in ((-0.2, pb, "C0", "pseudobulk paired t"),
                                   (0.2, nb, "C3", "NB GLM")):
        p = tab["p_any"].to_numpy()
        lo = tab["p_any_lo"].to_numpy()
        hi = tab["p_any_hi"].to_numpy()
        ax.bar(x + off, p, width=0.4, color=color, edgecolor="black", linewidth=0.5,
               label=label, yerr=[p - lo, hi - p], capsize=2, error_kw={"lw": 0.8})
    hot = [ct for ct in nb.index if nb.at[ct, "p_any_lo"] > 0.05]
    hot.sort(key=lambda ct: -nb.at[ct, "p_any"])
    for ct in hot:
        i = list(pb.index).index(ct)
        ax.text(i + 0.2, nb.at[ct, "p_any_hi"] + 0.03,
                f"{int(nb.at[ct, 'n_draws_with_any_fp'])}/{int(nb.at[ct, 'draws'])}",
                ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax.axhline(0.05, ls="--", lw=1, color="grey", label="0.05 reference")
    ax.set_xticks(x)
    ax.set_xticklabels([_short(t) for t in pb.index], rotation=35, ha="right",
                       fontsize=8)
    ax.set_ylim(0, 1.3)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    more = 100 * (nbsum["n_nb"] / nbsum["n_pb"] - 1)
    # Title honesty (2026-11-15): this said "where pseudobulk makes none", which
    # is false for CD8 T cells - pseudobulk makes 1 false call in 15 draws
    # there, and CD8 is one of the highlighted types. State the measured worst
    # pseudobulk rate among the highlighted types instead of asserting zero.
    pb_hot = [(int(pb.at[ct, "n_draws_with_any_fp"]), int(pb.at[ct, "draws"]))
              for ct in hot]
    pb_worst = max(pb_hot, default=(0, 0), key=lambda nd: nd[0] / max(nd[1], 1))
    if pb_worst[0] == 0:
        pb_clause = "where pseudobulk makes none"
    else:
        pb_clause = (f"where pseudobulk makes at most {pb_worst[0]}/{pb_worst[1]}"
                     if pb_worst[0] <= 1 else
                     f", where pseudobulk makes up to {pb_worst[0]}/{pb_worst[1]}")
    ax.set(ylabel="share of balanced null draws\nwith >= 1 false positive",
           title=f"NB GLM: {more.min():.0f}-{more.max():.0f}% more genes called, but false calls\n"
                 "on a balanced null in " + " and ".join(_short(ct) for ct in hot)
                 + ",\n" + pb_clause)
    ax.title.set_fontsize(9.5)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    _panel_letter(ax, "D")

    _save(fig, "at_a_glance.png",
          f"A: n = {len(sim)} simulated replicates and {len(real)} donor splits; bars: "
          "mean, whiskers: min-max. B: n = "
          f"{power['rep'].nunique()} replicates per method. C: n = {len(summ)} cell "
          "types. D: balanced donor-swap draws (35/35/35/70/15/15/15 by type), whiskers: "
          "Wilson 95% CI; the 0.05 line is the exact-null reference, not a bound.\n"
          "source: results/sim_null.csv, false_positives.csv, sim_power.csv, "
          "mixed_vs_others_summary.csv, mixed_null_draws.csv, "
          "nbglm_real_null_balanced.csv, nbglm_real_summary.csv")


FIGURES = {
    "null": fig_null,
    "power": fig_power,
    "sweep": fig_sample_sweep,
    "normalisation": fig_normalisation,
    "real_null": fig_real_null,
    "isg": fig_isg_ranking,
    "mixed_null": fig_mixed_null,
    "nbglm_null": fig_nbglm_null,
    "boundary": fig_boundary,
    "plan_space": fig_plan_space,
    "log2fc": fig_log2fc_agreement,
    "glance": fig_at_a_glance,
    "mega": fig_megakaryocyte,
    "nbglm_sim": fig_nbglm_sim,
}


if __name__ == "__main__":
    wanted = sys.argv[1:] or list(FIGURES)
    unknown = [w for w in wanted if w not in FIGURES]
    if unknown:
        raise SystemExit(f"unknown figure(s) {unknown}; choose from {list(FIGURES)}")
    for name in wanted:
        FIGURES[name]()
    print(f"figures written to {FIG}")
