"""Item 3: mixed-effects model vs per-cell and pseudobulk on the same data.

Three models are fitted to the identical CD4 T cell subset of Kang 2018:

  per-cell     : Wilcoxon on cells, every cell a replicate (wrong, for scale)
  pseudobulk   : paired t-test on donor x condition sums (7-8 df)
  mixed        : expression ~ condition + (1|donor), fitted per gene

The mixed model is the one being evaluated, so it gets the harder questions:
does it converge, what does it cost, and does it agree with the method that
is already trusted on this data.

Cost is measured, not estimated. The per-cell fit is O(genes x iterative
optimisation) and is run on a deliberate subsample so the wall-clock is
recorded rather than extrapolated: the number reported is the measured time
for the genes and cells actually fitted, plus the measured per-gene rate
used to state what the full matrix would have cost.

Run:
    python scripts/mixed_vs_others.py
"""

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from scpipe import data, de, mixed, pseudobulk  # noqa: E402

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

ALPHA = 0.05
MIN_CELLS = 25
CELL_TYPE = "CD4 T cells"

# per-cell fits are the expensive arm: cap them deliberately
PERCELL_GENES = 400
PERCELL_MAX_CELLS = 3200
NULL_GENES = 500
NULL_PERMS = 10
SEED = 0


def _subsample(adata, n_max, seed):
    """Take cells evenly across donor x condition, never at random globally.

    A uniform random subsample would let a donor with many cells dominate
    the fit; stratifying preserves the paired design, which is the only
    thing the random effect is estimated from.
    """
    rng = np.random.default_rng(seed)
    obs = adata.obs
    keys = obs["donor"].astype(str) + "||" + obs["condition"].astype(str)
    per = max(1, n_max // keys.nunique())
    idx = []
    for _, group in keys.groupby(keys).groups.items():
        g = np.asarray(group)
        idx.extend(rng.choice(g, size=min(per, len(g)), replace=False).tolist())
    return np.sort(np.asarray(idx))


def main():
    t_all = time.time()
    adata = data.load_kang()
    sub = adata[adata.obs["cell_type"] == CELL_TYPE].copy()
    print(f"{CELL_TYPE}: {sub.n_obs:,} cells x {sub.n_vars:,} genes, "
          f"{sub.obs['donor'].nunique()} donors")
    print(sub.obs.groupby(["condition"], observed=True).size().to_string())

    # ---- pseudobulk profiles, which two of the three arms share ---------
    mat, meta = pseudobulk.aggregate(sub, "donor", "condition",
                                     min_cells=MIN_CELLS)
    mat = pseudobulk.filter_genes(mat, min_count=10, min_samples=3)
    logcpm = pseudobulk.normalise(mat, log=True)
    print(f"\npseudobulk profiles: {logcpm.shape[0]} samples x "
          f"{logcpm.shape[1]:,} genes tested")

    # ---- arm 1: per-cell Wilcoxon --------------------------------------
    t0 = time.time()
    per_cell = de.single_cell_de(sub, "condition", "stim", "ctrl")
    t_percell = time.time() - t0
    per_cell = per_cell.loc[logcpm.columns]        # same genes for all arms
    print(f"per-cell Wilcoxon: {de.n_significant(per_cell):,} significant, "
          f"{t_percell:.1f}s")

    # ---- arm 2: pseudobulk paired t ------------------------------------
    t0 = time.time()
    pb = de.pseudobulk_de(sub, "condition", "donor", "stim", "ctrl",
                          paired=True, min_cells=MIN_CELLS)
    pb = pb.loc[logcpm.columns]
    t_pb = time.time() - t0
    print(f"pseudobulk paired t: {de.n_significant(pb):,} significant, "
          f"{t_pb:.1f}s")

    # ---- arm 3: mixed model per gene -----------------------------------
    t0 = time.time()
    mx = mixed.mixed_de(logcpm, "donor", "condition", "stim", "ctrl")
    t_mx = time.time() - t0
    n_conv = int(mx["converged"].sum())
    print(f"mixedlm per gene: {de.n_significant(mx):,} significant, "
          f"{t_mx:.1f}s, converged {n_conv:,}/{len(mx):,} "
          f"({100 * n_conv / len(mx):.1f}%)")
    print(f"  measured per-gene cost: "
          f"{1000 * t_mx / len(mx):.2f} ms")
    print(f"  donor variance component: median "
          f"{mx['donor_var'].median():.4f}, "
          f"{int((mx['donor_var'] < 1e-6).sum()):,} genes at the boundary")

    # ---- agreement -----------------------------------------------------
    def called(frame, col="qval"):
        return set(frame.index[frame[col] < ALPHA])

    s_pc, s_pb, s_mx = called(per_cell), called(pb), called(mx)
    both = pd.concat([pb["log2fc"].rename("pb"), mx["log2fc"].rename("mx")],
                     axis=1).dropna()
    r_lfc = stats.pearsonr(both["pb"], both["mx"])
    rho = stats.spearmanr(both["pb"], both["mx"])

    def jaccard(a, b):
        return len(a & b) / len(a | b) if (a | b) else np.nan

    print("\n--- calls at FDR 5% ---")
    print(f"  per-cell : {len(s_pc):,}")
    print(f"  pseudobulk: {len(s_pb):,}")
    print(f"  mixed    : {len(s_mx):,}  (of {n_conv:,} converged)")
    print(f"  Jaccard pseudobulk vs mixed : {jaccard(s_pb, s_mx):.3f}")
    print(f"  Jaccard per-cell vs mixed   : {jaccard(s_pc, s_mx):.3f}")
    print(f"  log2FC pearson r (pb, mx)   : {r_lfc[0]:.3f} "
          f"(p={r_lfc[1]:.2e})")
    print(f"  log2FC spearman rho         : {rho.statistic:.3f}")

    isg = [g for g in data.ISG_MARKERS if g in logcpm.columns]
    print(f"\nISG recovery of {len(isg)} markers present:")
    for name, frame in (("per-cell", per_cell), ("pseudobulk", pb),
                        ("mixed", mx)):
        sig = [g for g in isg if frame.loc[g, "qval"] < ALPHA]
        ranks = [int(np.where(frame.sort_values("qval").index == g)[0][0]) + 1
                 for g in isg]
        print(f"  {name:>10}: {len(sig)}/{len(isg)} significant, "
              f"median rank {int(np.median(ranks))}")

    # ---- null: swap the two aliquots of each donor ---------------------
    # Truth is "no difference", because the swap is applied identically to
    # every gene and leaves the donor structure intact. Any call is a false
    # positive. This is the test that decides between methods.
    rng = np.random.default_rng(SEED)
    null_genes = list(logcpm.columns[
        rng.choice(len(logcpm.columns), size=min(NULL_GENES,
                                                 len(logcpm.columns)),
                   replace=False)])
    donors = meta["donor"].to_numpy()
    cond = meta["condition"].to_numpy()
    # the per-cell arm is deliberately absent: swapping whole samples does
    # not expose its failure mode, which is treating cells as replicates.
    # That failure is measured by the simulated null in scripts/simulation_study.py.
    fps = {"pseudobulk": [], "mixed": []}
    t0 = time.time()
    for _ in range(NULL_PERMS):
        swap = rng.random(len(set(donors))) < 0.5
        swap_map = dict(zip(sorted(set(donors)), swap))
        flipped = np.array([("ctrl" if c == "stim" else "stim") if
                            swap_map[d] else c
                            for d, c in zip(donors, cond)])
        is_stim = flipped == "stim"
        A = logcpm.loc[is_stim, null_genes]
        B = logcpm.loc[~is_stim, null_genes]
        donors_a = meta.loc[is_stim, "donor"].to_numpy()
        donors_b = meta.loc[~is_stim, "donor"].to_numpy()
        shared = sorted(set(donors_a) & set(donors_b))
        ai = [np.where(donors_a == d)[0][0] for d in shared]
        bi = [np.where(donors_b == d)[0][0] for d in shared]
        Aa, Bb = A.iloc[ai], B.iloc[bi]
        _, p_t = stats.ttest_rel(Aa.to_numpy(), Bb.to_numpy(), axis=0)
        fps["pseudobulk"].append(int((de._bh(p_t) < ALPHA).sum()))

        # mixed on the same permuted labels, through the module so both
        # arms use the same solver and the same failure handling
        permuted = logcpm.copy()
        permuted.index = pd.Index(
            [f"{d}||{('stim' if s else 'ctrl')}"
             for d, s in zip(donors, is_stim)])
        pm = mixed.mixed_de(permuted[null_genes], "donor", "condition",
                            "stim", "ctrl")["pval"].to_numpy()
        fps["mixed"].append(int((de._bh(pm) < ALPHA).sum()))
    t_null = time.time() - t0

    null_df = pd.DataFrame({
        "method": list(fps),
        "mean_false_positives": [float(np.mean(v)) for v in fps.values()],
        "genes_tested": NULL_GENES,
        "perms": NULL_PERMS,
    })
    null_df["mean_fdp"] = null_df["mean_false_positives"] / NULL_GENES
    print(f"\n--- donor-swap null ({NULL_PERMS} permutations, "
          f"{len(null_genes)} genes, {t_null:.1f}s) ---")
    print(null_df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    null_df.to_csv(RESULTS / "mixed_null_donor_swap.csv", index=False)

    # ---- per-cell mixed model: the expensive arm ----------------------
    # Measured, not extrapolated: fit a fixed number of genes on a
    # donor-stratified cell subsample and report the actual seconds.
    idx = _subsample(sub, PERCELL_MAX_CELLS, SEED)
    cell_sub = sub[idx]
    Xc = cell_sub.X
    Xc = Xc.toarray() if hasattr(Xc, "toarray") else np.asarray(Xc)
    genes_pc = list(logcpm.columns[
        np.random.default_rng(SEED).choice(
            len(logcpm.columns), size=min(PERCELL_GENES, len(logcpm.columns)),
            replace=False)])
    col_idx = [cell_sub.var_names.get_loc(g) for g in genes_pc]
    vals = Xc[:, col_idx].astype(float)
    lib = Xc.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    vals = np.log2(vals / lib * 1e6 + 1.0)

    t0 = time.time()
    mx_cell = mixed.mixed_de_percell(
        vals, genes_pc, cell_sub.obs["donor"].to_numpy(),
        cell_sub.obs["condition"].to_numpy(), "stim", "ctrl")
    t_mx_cell = time.time() - t0
    per_gene = t_mx_cell / len(genes_pc)
    projected = per_gene * len(logcpm.columns)
    print(f"\n--- per-cell mixed model (down-sampled) ---")
    print(f"  cells fitted: {len(idx):,} of {sub.n_obs:,} "
          f"({100 * len(idx) / sub.n_obs:.0f}%)")
    print(f"  genes fitted: {len(genes_pc)} of {len(logcpm.columns):,}")
    print(f"  measured: {t_mx_cell:.1f}s = {1000 * per_gene:.1f} ms/gene")
    print(f"  projected full gene set at that cost: {projected:.0f}s "
          f"({projected / 60:.1f} min)")
    print(f"  converged {int(mx_cell['converged'].sum())}/{len(mx_cell)}")
    print(f"  significant {de.n_significant(mx_cell)}/{len(mx_cell)} genes")

    # agreement with the pseudobulk mixed model on the shared genes
    shared_genes = [g for g in genes_pc if g in mx.index]
    cmp_small = pd.DataFrame({
        "mx_pseudobulk": mx.loc[shared_genes, "log2fc"],
        "mx_percell": mx_cell.loc[shared_genes, "log2fc"],
    }).dropna()
    r_cell = stats.pearsonr(cmp_small["mx_pseudobulk"],
                            cmp_small["mx_percell"])
    agree = int(((mx.loc[shared_genes, "qval"] < ALPHA)
                 == (mx_cell.loc[shared_genes, "qval"] < ALPHA)).sum())
    print(f"  log2FC pearson r vs pseudobulk mixed: {r_cell[0]:.3f} "
          f"on {len(cmp_small)} genes")
    print(f"  same significance call: {agree}/{len(shared_genes)} "
          f"({100 * agree / len(shared_genes):.0f}%)")

    # ---- write ------------------------------------------------------
    mx.sort_values("qval").to_csv(RESULTS / "mixed_cd4.csv")
    mx_cell.sort_values("qval").to_csv(RESULTS / "mixed_percell_cd4.csv")

    summary = pd.DataFrame([{
        "cell_type": CELL_TYPE,
        "n_cells": int(sub.n_obs),
        "n_genes_tested": int(len(logcpm)),
        "percell_sig": len(s_pc),
        "pseudobulk_sig": len(s_pb),
        "mixed_sig": len(s_mx),
        "mixed_converged": n_conv,
        "mixed_pct_converged": round(100 * n_conv / len(mx), 2),
        "jaccard_pb_mx": round(jaccard(s_pb, s_mx), 4),
        "jaccard_percell_mx": round(jaccard(s_pc, s_mx), 4),
        "lfc_pearson_pb_mx": round(float(r_lfc[0]), 4),
        "lfc_spearman_pb_mx": round(float(rho.statistic), 4),
        "mixed_ms_per_gene_pseudobulk": round(1000 * t_mx / len(mx), 2),
        "percell_mixed_cells_fitted": int(len(idx)),
        "percell_mixed_genes_fitted": int(len(genes_pc)),
        "percell_mixed_s": round(t_mx_cell, 1),
        "percell_mixed_projected_full_s": round(projected, 1),
        "percell_mixed_pct_converged": round(
            100 * mx_cell["converged"].mean(), 2),
        "null_fp_pseudobulk": round(float(np.mean(fps["pseudobulk"])), 2),
        "null_fp_mixed": round(float(np.mean(fps["mixed"])), 2),
        "null_genes": NULL_GENES,
        "timestamp": pd.Timestamp.now().isoformat(timespec="seconds"),
    }])
    summary.to_csv(RESULTS / "mixed_vs_others_summary.csv", index=False)
    print("\n--- summary ---")
    print(summary.T.to_string())
    print(f"\ntotal {time.time() - t_all:.1f}s; wrote {RESULTS}")


if __name__ == "__main__":
    main()
