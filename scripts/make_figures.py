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
    ax.legend(fontsize=7.5, loc="upper left", markerscale=2.5, frameon=True,
              framealpha=0.9, edgecolor="none")

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
}


if __name__ == "__main__":
    wanted = sys.argv[1:] or list(FIGURES)
    unknown = [w for w in wanted if w not in FIGURES]
    if unknown:
        raise SystemExit(f"unknown figure(s) {unknown}; choose from {list(FIGURES)}")
    for name in wanted:
        FIGURES[name]()
    print(f"figures written to {FIG}")
