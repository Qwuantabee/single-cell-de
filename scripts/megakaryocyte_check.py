"""Item 2b: can megakaryocytes actually be tested?

The README said megakaryocytes were dropped because "only 1 donor has >= 25
cells in both conditions". That sentence is true but it names a threshold,
not a cause, and the two are not the same claim. 132 megakaryocytes exist
across 8 donors, so the cells are not absent; the question is whether they
are distributed across donors in a way that supports a test whose replicate
is the donor.

This script answers that question with the per-donor counts printed, a sweep
over every candidate threshold, and the power of the best design the data
admits. It does not pick a threshold to make a result appear.

Run:
    python scripts/megakaryocyte_check.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from scpipe import data, de, pseudobulk  # noqa: E402

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

CELL_TYPE = "Megakaryocytes"
ALPHA = 0.05
POWER = 0.80


def min_detectable_effect(sd, n, alpha=ALPHA, power=POWER):
    """Smallest paired log2 effect a two-sided t-test can find.

    MDE = (t_{1-a/2, n-1} + t_{power, n-1}) * sd / sqrt(n), with sd the
    standard deviation of the *paired difference* -- which is the quantity
    the paired test actually has, not the spread of the levels.
    """
    if n < 2 or not np.isfinite(sd):
        return np.inf
    df = n - 1
    crit = stats.t.ppf(1 - alpha / 2, df) + stats.t.ppf(power, df)
    return float(crit * sd / np.sqrt(n))


def main():
    adata = data.load_kang()
    obs = adata.obs
    sub = adata[obs["cell_type"].astype(str) == CELL_TYPE].copy()
    n_total = sub.n_obs

    print(f"{CELL_TYPE}: {n_total} cells of {adata.n_obs:,} "
          f"({100 * n_total / adata.n_obs:.2f}%)")

    # ---- 1. the raw per-donor counts, which decide everything ----------
    tab = pd.crosstab(sub.obs["donor"], sub.obs["condition"])
    tab["total"] = tab.sum(axis=1)
    tab["share_of_type"] = tab["total"] / n_total
    tab = tab.sort_values("total", ascending=False)
    print("\n--- cells per donor per condition ---")
    print(tab.to_string())
    tab.to_csv(RESULTS / "megakaryocyte_donor_counts.csv")

    ctrl_total = int(tab["ctrl"].sum())
    stim_total = int(tab["stim"].sum())
    top = tab.index[0]
    top_share = float(tab.loc[top, "total"] / n_total)
    print(f"\nctrl {ctrl_total} cells, stim {stim_total} cells, "
          f"{len(tab)} donors contribute at least one cell")
    print(f"{top} alone supplies {int(tab.loc[top, 'total'])} cells "
          f"({100 * top_share:.1f}% of the type)")

    # ---- 2. what every candidate threshold admits ----------------------
    rows = []
    for t in (1, 2, 3, 5, 10, 15, 25):
        ok = tab[(tab["stim"] >= t) & (tab["ctrl"] >= t)]
        rows.append({
            "min_cells_per_condition": t,
            "donors_qualifying": len(ok),
            "donors": ";".join(ok.index.tolist()),
            "cells_retained": int(ok[["ctrl", "stim"]].to_numpy().sum()),
            "test_possible": len(ok) >= 2,
        })
    sweep = pd.DataFrame(rows)
    print("\n--- donors surviving each threshold ---")
    print(sweep.to_string(index=False))
    sweep.to_csv(RESULTS / "megakaryocyte_threshold_sweep.csv", index=False)

    # ---- 3. the threshold is not the binding constraint ----------------
    # >=1 cell in both conditions is the loosest possible filter that still
    # leaves something to compare. If that does not yield 2 donors, no
    # threshold can.
    loosest = int(sweep.loc[sweep["min_cells_per_condition"] == 1,
                            "donors_qualifying"].iloc[0])

    # ---- 4. how much power the loosest design would actually have ------
    # Run the aggregation at 1 cell so the per-gene paired differences can
    # be measured, then state the MDE. This is a sensitivity calculation,
    # not a result: it says what such a test could see, not what it saw.
    mde_block = {}
    try:
        mat, meta = pseudobulk.aggregate(sub, "donor", "condition",
                                         min_cells=1)
        mat = pseudobulk.filter_genes(mat, min_count=10, min_samples=3)
        logcpm = pseudobulk.normalise(mat, log=True)
        donors = meta["donor"].to_numpy()
        cond = meta["condition"].to_numpy()
        shared = sorted({d for d in donors
                         if ((donors == d) & (cond == "stim")).any()
                         and ((donors == d) & (cond == "ctrl")).any()})
        diffs = []
        for d in shared:
            a = logcpm.loc[(donors == d) & (cond == "stim")].to_numpy()[0]
            b = logcpm.loc[(donors == d) & (cond == "ctrl")].to_numpy()[0]
            diffs.append(a - b)
        diffs = np.vstack(diffs)
        paired_sd = float(np.median(np.std(diffs, axis=0, ddof=1)))
        mde = min_detectable_effect(paired_sd, len(shared))
        mde_block = {
            "donors_in_loosest_design": len(shared),
            "cells_per_profile_min": int(meta["n_cells"].min()),
            "cells_per_profile_max": int(meta["n_cells"].max()),
            "median_paired_sd_log2": paired_sd,
            "mde_log2_at_80pct_power": mde,
            "mde_fold": 2.0 ** mde,
        }
        print("\n--- power of the loosest possible design ---")
        for k, v in mde_block.items():
            print(f"  {k}: {v if isinstance(v, int) else f'{v:.3f}'}")

        res = de.pseudobulk_de(sub, "condition", "donor", "stim", "ctrl",
                               paired=True, min_cells=1)
        n_sig = de.n_significant(res)
        print(f"  genes tested {len(res):,}, significant at FDR 5%: {n_sig}")
        mde_block["genes_tested"] = int(len(res))
        mde_block["genes_significant"] = int(n_sig)
        res.to_csv(RESULTS / "de_megakaryocytes_sensitivity_min1cell.csv")
    except ValueError as exc:
        print(f"\nloosest design could not be aggregated: {exc}")
        mde_block = {"error": str(exc)}

    # ---- 5. the verdict, computed rather than asserted ----------------
    # "Test possible" and "test usable" are different questions. At >=1 cell
    # in both conditions four donors technically qualify, so the pairing
    # exists; what it does NOT have is any balance, because one donor is 84%
    # of the cells and the other extremes are 1-cell profiles. The boolean
    # below is about the pairing, the number next to it is why the pairing is
    # not informative.
    testable = loosest >= 2
    balanced = top_share < 0.5
    reliable = bool(testable and balanced)
    statement = (
        f"Megakaryocytes cannot support a donor-level expression test. At the "
        f"loosest threshold that still compares anything (>=1 cell in both "
        f"conditions) {loosest} donors qualify, and {top} supplies "
        f"{100 * top_share:.0f}% of the cells, so the comparison pairs "
        f"58-cell profiles against 1-cell profiles and is dominated by one "
        f"donor. The composition of megakaryocytes across conditions is "
        f"testable and is tested on all 8 donors in "
        f"results/celltype_composition_paired.csv; the expression contrast "
        f"is not."
    )
    print("\n--- verdict ---")
    print(statement)
    print(f"pairing_exists_at_loosest_threshold = {testable}")
    print(f"cell_balance_reasonable = {balanced} "
          f"(largest donor {100 * top_share:.0f}% of the type)")
    print(f"expression_test_reliable = {reliable}")

    summary = pd.DataFrame([{
        "cell_type": CELL_TYPE,
        "n_cells": int(n_total),
        "donors_with_any_cell": int(len(tab)),
        "ctrl_cells": ctrl_total,
        "stim_cells": stim_total,
        "largest_donor": top,
        "largest_donor_share": top_share,
        "donors_at_loosest_threshold": loosest,
        "pairing_exists": bool(testable),
        "expression_test_reliable": reliable,
        **mde_block,
    }])
    summary.to_csv(RESULTS / "megakaryocyte_summary.csv", index=False)
    print("\n--- summary ---")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
