"""The positive control: does pseudobulk find the interferon response?

A method that calls nothing on null data is easy to build -- return no
results and the false positive rate is zero. The other half of the claim is
that it still finds what is really there.

Interferon-beta stimulation produces one of the best characterised
transcriptional responses in immunology. ISG15, IFIT1/3, MX1, OAS1 and
friends should be at the top of any correct analysis, and if they are not,
the method is broken regardless of how well behaved it looked on the null.

Run:
    python scripts/real_analysis.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from scpipe import data, de  # noqa: E402

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

ALPHA = 0.05
MIN_CELLS = 25


def main():
    adata = data.load_kang()
    data.describe(adata)

    cell_types = (adata.obs["cell_type"].value_counts()
                  .loc[lambda s: s > 500].index.tolist())
    print(f"\nanalysing {len(cell_types)} cell types with >500 cells")

    all_rows = []
    for ct in cell_types:
        sub = adata[adata.obs["cell_type"] == ct].copy()

        # every donor must contribute to both conditions for a paired test
        tab = pd.crosstab(sub.obs["donor"], sub.obs["condition"])
        usable = tab[(tab >= MIN_CELLS).all(axis=1)].index.tolist()
        if len(usable) < 4:
            print(f"\n{ct}: only {len(usable)} usable donors, skipping")
            continue
        sub = sub[sub.obs["donor"].isin(usable)].copy()

        pb = de.pseudobulk_de(sub, "condition", "donor", "stim", "ctrl",
                              paired=True, min_cells=MIN_CELLS)
        sc = de.single_cell_de(sub, "condition", "stim", "ctrl")

        n_pb, n_sc = de.n_significant(pb, ALPHA), de.n_significant(sc, ALPHA)
        print(f"\n{ct}  ({sub.n_obs:,} cells, {len(usable)} donors)")
        print(f"  pseudobulk {n_pb:>5} / {len(pb):>5} genes significant")
        print(f"  per-cell   {n_sc:>5} / {len(sc):>5} genes significant")

        # Two rankings of the same significant genes.
        #
        # By q-value is the reflex, and with 8 donors it is close to
        # useless for prioritising: a paired t-test on 7 degrees of freedom
        # has a floor on how small p can get, so hundreds of genes pile up
        # against it and their order is decided by noise. A gene with a
        # modest but very consistent change outranks one with a tenfold
        # change that wobbles slightly between donors.
        #
        # By effect size among genes that already passed the significance
        # filter answers the question actually being asked: of the things
        # that changed, which changed most.
        sig = pb[pb["qval"] < ALPHA]
        by_q = sig.sort_values("qval")
        by_effect = sig.reindex(
            sig["log2fc"].abs().sort_values(ascending=False).index)

        found_q = [g for g in data.ISG_MARKERS if g in by_q.index[:12]]
        found_e = [g for g in data.ISG_MARKERS if g in by_effect.index[:12]]
        print(f"  smallest q reached: {pb['qval'].min():.2e}")
        print(f"  top 12 by q-value:    {', '.join(by_q.index[:12])}")
        print(f"      known ISGs: {len(found_q)}/12")
        print(f"  top 12 by effect size: {', '.join(by_effect.index[:12])}")
        print(f"      known ISGs: {len(found_e)}/12 -> {found_e}")

        def ranks(frame):
            return [int(np.where(frame.index == g)[0][0]) + 1
                    for g in data.ISG_MARKERS if g in frame.index]

        rq, re_ = ranks(by_q), ranks(by_effect)
        all_rows.append({
            "cell_type": ct, "n_cells": int(sub.n_obs),
            "n_donors": len(usable),
            "pseudobulk_sig": n_pb, "per_cell_sig": n_sc,
            "genes_tested_pb": len(pb), "genes_tested_sc": len(sc),
            "min_qval": float(pb["qval"].min()),
            "isgs_top12_by_q": len(found_q),
            "isgs_top12_by_effect": len(found_e),
            "median_isg_rank_by_q": int(np.median(rq)) if rq else -1,
            "median_isg_rank_by_effect": int(np.median(re_)) if re_ else -1,
        })
        pb.sort_values("qval").to_csv(
            RESULTS / f"de_{ct.replace(' ', '_').replace('+', '')}.csv")

    summary = pd.DataFrame(all_rows)
    summary.to_csv(RESULTS / "real_analysis_summary.csv", index=False)

    print("\n--- summary ---")
    print(summary.to_string(index=False))
    print(f"\nTotal significant calls: pseudobulk "
          f"{summary['pseudobulk_sig'].sum():,}, "
          f"per-cell {summary['per_cell_sig'].sum():,}")
    print(f"\nMedian ISG rank, by q-value:    "
          f"{summary['median_isg_rank_by_q'].median():.0f}")
    print(f"Median ISG rank, by effect size: "
          f"{summary['median_isg_rank_by_effect'].median():.0f}")
    print(f"\nwrote {RESULTS / 'real_analysis_summary.csv'}")


if __name__ == "__main__":
    main()
