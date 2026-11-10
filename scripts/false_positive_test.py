"""Measure the false positive rate of per-cell vs pseudobulk differential
expression, on data where the true answer is known.

The design is a permutation null. Take control cells only -- no stimulation,
so there is no real biological difference anywhere in the data. Randomly
split the 8 donors into two groups of 4 and call them "condition A" and
"condition B". Any gene that comes out significant is, by construction, a
false positive.

Repeat over many random splits and count. A method controlling FDR at 5%
should call almost nothing, because there is nothing to call.

Run:
    python scripts/false_positive_test.py
    python scripts/false_positive_test.py --n-splits 50 --cell-type "CD14+ Monocytes"
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from scpipe import data, de  # noqa: E402

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)


def run(n_splits=20, cell_type=None, alpha=0.05, seed=0, max_genes=None):
    adata = data.load_kang()
    data.describe(adata)

    # control cells only: there is no real difference between any two
    # groups of these donors
    ctrl_label = _pick_control_label(adata)
    sub = adata[adata.obs["condition"] == ctrl_label].copy()
    print(f"\nusing {sub.n_obs:,} {ctrl_label!r} cells")

    if cell_type and "cell_type" in sub.obs.columns:
        sub = sub[sub.obs["cell_type"].astype(str) == cell_type].copy()
        print(f"restricted to {cell_type}: {sub.n_obs:,} cells")

    if max_genes:
        # keep the most expressed genes, to keep the per-cell test tractable
        expr = np.asarray(sub.X.sum(axis=0)).ravel()
        keep = np.argsort(-expr)[:max_genes]
        sub = sub[:, np.sort(keep)].copy()
        print(f"restricted to the top {max_genes} expressed genes")

    donors = sorted(sub.obs["donor"].unique())
    print(f"{len(donors)} donors: {donors}")
    if len(donors) < 6:
        raise SystemExit("need at least 6 donors to split into two groups")

    rng = np.random.default_rng(seed)
    rows = []

    for i in range(n_splits):
        shuffled = list(rng.permutation(donors))
        half = len(shuffled) // 2
        group_a, group_b = set(shuffled[:half]), set(shuffled[half:2 * half])

        fake = sub.obs["donor"].map(
            lambda d: "A" if d in group_a else ("B" if d in group_b else None))
        keep = fake.notna().to_numpy()
        trial = sub[keep].copy()
        trial.obs["fake_condition"] = fake[keep].astype(str).to_numpy()

        sc_res = de.single_cell_de(trial, "fake_condition", "A", "B")
        pb_res = de.pseudobulk_de(trial, "fake_condition", "donor", "A", "B")

        n_sc = de.n_significant(sc_res, alpha)
        n_pb = de.n_significant(pb_res, alpha)
        rows.append({"split": i, "per_cell": n_sc, "pseudobulk": n_pb,
                     "n_genes_sc": len(sc_res), "n_genes_pb": len(pb_res),
                     "n_cells": int(trial.n_obs)})
        print(f"  split {i + 1:>2}/{n_splits}: per-cell {n_sc:>5} "
              f"false positives, pseudobulk {n_pb:>3}", flush=True)

    df = pd.DataFrame(rows)
    out = RESULTS / "false_positives.csv"
    df.to_csv(out, index=False)

    print(f"\n--- false positives at FDR {alpha} over {n_splits} null splits ---")
    for method in ("per_cell", "pseudobulk"):
        v = df[method]
        n_tested = df["n_genes_sc" if method == "per_cell" else "n_genes_pb"]
        rate = (v / n_tested).mean()
        print(f"  {method:<12} median {v.median():>6.0f}  mean {v.mean():>7.1f}  "
              f"max {v.max():>6.0f}   ({100 * rate:.2f}% of genes tested)")
    print(f"\nwrote {out}")
    return df


def _pick_control_label(adata):
    labels = sorted(adata.obs["condition"].unique())
    for cand in ("ctrl", "control", "CTRL", "unstim", "NT"):
        if cand in labels:
            return cand
    print(f"could not identify a control label among {labels}, using {labels[0]!r}")
    return labels[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-splits", type=int, default=20)
    ap.add_argument("--cell-type", default=None)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--max-genes", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    run(n_splits=a.n_splits, cell_type=a.cell_type, alpha=a.alpha,
        seed=a.seed, max_genes=a.max_genes)
