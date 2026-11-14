"""Item 2: aggregate within cell type before testing.

The question this answers is whether the expression changes reported for the
interferon response are expression changes at all. In this dataset the
conditions also differ in which cells were recovered, so a pooled test can
report the composition difference as if it were transcription. Stratifying
by cell type removes that by construction, and the difference between the
two answers is the size of the confound.

Run:
    python scripts/celltype_analysis.py
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from scpipe import data, de, pseudobulk  # noqa: E402

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

ALPHA = 0.05
MIN_CELLS = 25
MIN_DONORS = 4


def skipped_detail(adata, skipped, group_a="stim", group_b="ctrl"):
    """Print the per-donor counts behind every skipped cell type.

    A skip is a claim about the data, so it has to be shown rather than
    asserted: "1 donor meets the threshold" is only convincing next to the
    counts that made it true, and it is the only way a reader can tell a
    threshold problem from a data problem.
    """
    if not skipped:
        return
    obs = adata.obs
    print("\n--- why each skipped type was skipped ---")
    for ct, why in skipped.items():
        sub = obs[obs["cell_type"].astype(str) == ct]
        tab = pd.crosstab(sub["donor"], sub["condition"])
        tab["total"] = tab.sum(axis=1)
        tab = tab.sort_values("total", ascending=False)
        n_loose = int(((tab[group_a] >= 1) & (tab[group_b] >= 1)).sum())
        print(f"{ct}: {why}")
        print(f"  {len(sub)} cells, {len(tab)} donors, "
              f"{n_loose} donors with >=1 cell in both conditions")
        print("  " + tab.to_string().replace("\n", "\n  "))


def main(min_cells=MIN_CELLS, min_donors=MIN_DONORS):
    t0 = time.time()
    adata = data.load_kang()
    data.describe(adata)

    # ---- 1. does the cell composition move between conditions? ----------
    props, shift = pseudobulk.cell_type_composition(adata)
    shift = shift.sort_values("log2_ratio")
    print("\n--- cell type composition per donor sample ---")
    print(shift.to_string(float_format=lambda v: f"{v:.3f}"))
    shift.to_csv(RESULTS / "celltype_composition_shift.csv")
    props.to_csv(RESULTS / "celltype_composition_proportions.csv")

    # is the composition difference bigger than donor-to-donor noise in it?
    cond = props.index.get_level_values("condition")
    comp_t = {}
    for ct in props.columns:
        a = props.loc[cond == "stim", ct]
        b = props.loc[cond == "ctrl", ct]
        from scipy import stats
        t, p = stats.ttest_rel(a.to_numpy(), b.to_numpy())
        comp_t[ct] = {"t": t, "pval": p,
                      "mean_stim": a.mean(), "mean_ctrl": b.mean()}
    comp_df = pd.DataFrame(comp_t).T.sort_values("pval")
    print("\npaired t-test on composition (stim vs ctrl, per donor):")
    print(comp_df.to_string(float_format=lambda v: f"{v:.4g}"))
    comp_df.to_csv(RESULTS / "celltype_composition_paired.csv")
    max_shift = float(shift["log2_ratio"].abs().max())
    print(f"\nlargest composition shift: {max_shift:.2f} log2 "
          f"({2 ** max_shift:.1f}x)")

    # ---- 2. pooled test, all cells of a donor together ------------------
    pooled = de.pseudobulk_de(adata, "condition", "donor", "stim", "ctrl",
                              paired=True, min_cells=min_cells)
    pooled_sig = set(pooled.index[pooled["qval"] < ALPHA])
    print(f"\npooled (all cells): {len(pooled_sig):,} of {len(pooled):,} "
          f"genes significant")
    pooled.sort_values("qval").to_csv(
        RESULTS / "de_pooled_all_cells.csv")

    # ---- 3. stratified test --------------------------------------------
    t1 = time.time()
    results, _, skipped = pseudobulk.pseudobulk_by_celltype(
        adata, group_a="stim", group_b="ctrl", min_cells=min_cells,
        min_donors=min_donors)
    elapsed = time.time() - t1
    print(f"\nstratified test (min_cells={min_cells}, "
          f"min_donors={min_donors}) on "
          f"{results.index.get_level_values('cell_type').nunique()}"
          f" cell types in {elapsed:.1f}s")
    skipped_detail(adata, skipped)
    for ct, why in skipped.items():
        print(f"  skipped {ct}: {why}")

    counts_local = pseudobulk.significant_by_celltype(results, ALPHA, "qval")
    counts_global = pseudobulk.significant_by_celltype(
        results, ALPHA, "qval_global")
    union = set(results.index[results["qval"] < ALPHA].get_level_values("gene"))
    union_global = set(
        results.index[results["qval_global"] < ALPHA].get_level_values("gene"))

    per_type = pd.DataFrame({
        "genes_tested": results.groupby(level="cell_type").size(),
        "n_donors": results.groupby(level="cell_type")["n_donors"].first(),
        "sig_within_type_bh": counts_local,
        "sig_global_bh": counts_global,
    })
    print("\n--- per cell type ---")
    print(per_type.to_string())
    per_type.to_csv(RESULTS / "celltype_stratified_counts.csv")

    results.to_csv(RESULTS / "de_stratified_by_celltype.csv")

    # ---- 4. what the pooled test adds, and what it is ------------------
    novel = pooled_sig - union
    lost = union - pooled_sig
    print(f"\ndistinct genes significant, stratified (within-type BH): "
          f"{len(union):,}")
    print(f"pooled: {len(pooled_sig):,}")
    print(f"pooled-only genes (not significant in any cell type): "
          f"{len(novel):,}")
    print(f"significant in >=1 cell type but not pooled: {len(lost):,}")

    # A composition artefact should still be significant at the pooled level
    # and should not be significant within any cell type. The gene's pooled
    # fold change should also be predictable from the composition shift:
    # when a cell type becomes more common, its markers rise in the pooled
    # comparison without any change in their own expression.
    top_ct = shift["log2_ratio"].idxmax()
    if novel:
        nov = pooled.loc[sorted(novel)].sort_values("log2fc", ascending=False)
        print(f"\ntop pooled-only genes by log2fc:")
        print(nov.head(15).to_string(
            columns=["log2fc", "qval"], float_format=lambda v: f"{v:.3g}"))

        # are the pooled-only hits the markers of the type that expanded?
        var_names = adata.var_names
        marker_scores = {}
        for ct in results.index.get_level_values("cell_type").unique():
            sub = adata[(adata.obs["cell_type"] == ct)
                        & (adata.obs["condition"] == "ctrl")]
            X = sub.X
            mean = np.asarray(X.mean(axis=0)).ravel() if not hasattr(X, "toarray") \
                else np.asarray(X.toarray().mean(axis=0)).ravel()
            others = adata[(adata.obs["cell_type"] != ct)
                           & (adata.obs["condition"] == "ctrl")]
            Xo = others.X
            mean_o = np.asarray(Xo.mean(axis=0)).ravel() if not hasattr(Xo, "toarray") \
                else np.asarray(Xo.toarray().mean(axis=0)).ravel()
            score = pd.Series(mean - mean_o, index=var_names)
            marker_scores[ct] = score
        markers = pd.DataFrame(marker_scores)
        for ct in per_type.index:
            top = markers[ct].sort_values(ascending=False).head(30).index
            overlap = len(set(top) & novel)
            print(f"  pooled-only genes among top 30 markers of {ct}: "
                  f"{overlap}/30")

    # ---- 5. the known interferon genes --------------------------------
    summary_rows = []
    for ct in per_type.index:
        sub = results.loc[ct]
        sig = sub[sub["qval"] < ALPHA]
        found = [g for g in data.ISG_MARKERS if g in sig.index]
        ranks = [int(np.where(sig.sort_values("qval").index == g)[0][0]) + 1
                 for g in data.ISG_MARKERS if g in sig.index]
        summary_rows.append({
            "cell_type": ct,
            "n_isg_found": len(found),
            "median_isg_rank": int(np.median(ranks)) if ranks else -1,
            "sig": int(len(sig)),
        })
    isg = pd.DataFrame(summary_rows).set_index("cell_type")
    print("\n--- known interferon-stimulated genes ---")
    print(isg.to_string())
    isg.to_csv(RESULTS / "celltype_isg_recovery.csv")

    pooled_found = [g for g in data.ISG_MARKERS if g in pooled_sig]
    print(f"\npooled test finds {len(pooled_found)}/{len(data.ISG_MARKERS)} ISGs")

    summary = pd.DataFrame([{
        "n_cells": int(adata.n_obs),
        "n_cell_types_tested": len(per_type),
        "max_composition_shift_log2": max_shift,
        "pooled_sig": len(pooled_sig),
        "stratified_union_sig": len(union),
        "stratified_union_sig_global_bh": len(union_global),
        "pooled_only_sig": len(novel),
        "stratified_only_sig": len(lost),
        "pooled_sig_marker_share": np.nan,
        "pooled_finds_isgs": len(pooled_found),
        "stratified_finds_isgs_anywhere": int(
            (isg["n_isg_found"] > 0).sum()),
        "elapsed_s": round(time.time() - t0, 1),
    }])
    summary.to_csv(RESULTS / "celltype_analysis_summary.csv", index=False)
    print("\n--- summary ---")
    print(summary.to_string(index=False))
    print(f"\nwrote {RESULTS}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-cells", type=int, default=MIN_CELLS,
                    help="cells required per donor per condition (default 25)")
    ap.add_argument("--min-donors", type=int, default=MIN_DONORS,
                    help="donors required per cell type (default 4)")
    args = ap.parse_args()
    main(min_cells=args.min_cells, min_donors=args.min_donors)
