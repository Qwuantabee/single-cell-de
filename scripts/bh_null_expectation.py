"""What Benjamini-Hochberg rejects on an all-null family, by simulation.

`scripts/mixed_vs_others.py` scores false positives against a reference of
about 0.05 rejected genes per family at q = 0.05. That figure used to be
described as "measured by simulation in scripts/false_positive_test.py",
which is a control-cell split experiment and contains no such simulation.
This script is the simulation: m independent uniform p-values (an exact null
for independent tests), BH at q, repeated, and the mean number of rejections
and P(at least one rejection) are recorded with Monte-Carlo standard errors.

Theory it should reproduce: for independent null p-values BH controls the
family-wise error rate at exactly q, so P(any rejection) = q = 0.05, and the
mean number of rejections is slightly above that because a draw that rejects
at all usually rejects more than one gene.

Caveat that matters for how the number is used: genes in one dataset are
not independent, and the swap null in mixed_vs_others.py uses a fixed gene
set across draws, so 0.05 is the reference for independent exact nulls, not
a derived bound for that experiment.

Run:
    python scripts/bh_null_expectation.py            # m=500, 5000 draws
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"


def bh_rejections(p, q=0.05):
    """Number of BH rejections per row of a (draws x m) p-value matrix."""
    p = np.sort(np.asarray(p, dtype=float), axis=1)
    m = p.shape[1]
    ok = p <= (np.arange(1, m + 1) / m) * q
    # largest k with p_(k) <= k q / m; zero when none qualifies
    ks = np.where(ok.any(axis=1), m - np.argmax(ok[:, ::-1], axis=1), 0)
    return ks


def simulate(m=500, draws=5000, q=0.05, seed=0):
    """Mean rejections and P(any) over iid uniform null families."""
    rng = np.random.default_rng(seed)
    k = bh_rejections(rng.random((draws, m)), q)
    any_ = (k > 0).astype(float)
    return {
        "m": m, "draws": draws, "q": q, "seed": seed,
        "mean_rejections": float(k.mean()),
        "mean_rejections_se": float(k.std(ddof=1) / np.sqrt(draws)),
        "p_any": float(any_.mean()),
        "p_any_se": float(np.sqrt(any_.mean() * (1 - any_.mean()) / draws)),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--m", type=int, default=500)
    ap.add_argument("--draws", type=int, default=5000)
    ap.add_argument("--q", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)
    row = simulate(args.m, args.draws, args.q, args.seed)
    RESULTS.mkdir(exist_ok=True)
    pd.DataFrame([row]).to_csv(RESULTS / "bh_null_expectation.csv", index=False)
    print(pd.DataFrame([row]).to_string(index=False))
    return row


if __name__ == "__main__":
    main()
