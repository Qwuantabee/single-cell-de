"""Does a proper count model beat a t-test on log counts at small donor counts?

The README's first "still to do" item asked this and it had never been
measured. The paired t-test on log2 CPM (`de.pseudobulk_de`) assumes the
per-donor values are Gaussian with constant variance. A negative-binomial
GLM (`nbglm.nbglm_de`) instead models the counts, estimates a per-gene
dispersion, and runs a Wald test on the condition coefficient -- the
DESeq2/edgeR approach.

Design. Same simulator as the rest of the study (`simulate.simulate`), which
draws from a negative binomial with a per-donor and per-sample random
effect. The generator is therefore already the model the GLM assumes, which
is the friendliest possible case for it; if it does not win here it will not
win on real data either.

  4 donor counts x 10 replicates x 2 methods, 300 genes, 30 truly changed
  with a log2 effect of 1.0 (2x). Both methods get the same pseudobulk
  matrix, the same gene filter and the same BH correction.

Reported per cell of the grid: mean observed FDP (false discovery proportion
among calls, i.e. what the FDR is supposed to bound at 0.05), mean power,
and the number of genes called. FDP is reported as measured rather than as
"1 - precision" so the numbers can be read directly against alpha.

Run:
    python scripts/nbglm_study.py
"""

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from scpipe import de, nbglm, simulate  # noqa: E402

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

ALPHA = 0.05
N_REPS = 10
DONOR_GRID = [3, 6, 12, 24]
N_GENES = 300
N_DE = 30


def one_rep(n_donors, seed):
    adata, truth = simulate.simulate(n_genes=N_GENES, n_donors=n_donors,
                                     n_de=N_DE, effect=1.0, seed=seed)
    rows = []
    pb = de.pseudobulk_de(adata, "condition", "donor", "stim", "ctrl",
                          paired=True, min_count=5, min_samples=3)
    t = simulate.confusion(pb, truth, ALPHA)
    rows.append({"method": "paired_t_logcpm", **t, "n_tested": len(pb)})

    nb = nbglm.nbglm_de(adata, "condition", "donor", "stim", "ctrl",
                        paired=True, min_count=5, min_samples=3)
    g = simulate.confusion(nb, truth, ALPHA)
    rows.append({"method": "nb_glm", **g, "n_tested": len(nb)})
    return rows


def main():
    out = []
    t0 = time.time()
    for n_donors in DONOR_GRID:
        acc = []
        for rep in range(N_REPS):
            rows = one_rep(n_donors, seed=rep)
            for r in rows:
                r.update(donors=n_donors, rep=rep)
            acc.extend(rows)
            print(f"  donors {n_donors:>2} rep {rep + 1:>2}: "
                  + ", ".join(f"{r['method']} fdp={r['fdp']:.2f} "
                              f"pow={r['power']:.2f}" for r in rows),
                  flush=True)
        out.extend(acc)

    df = pd.DataFrame(out)
    df.to_csv(RESULTS / "sim_nbglm.csv", index=False)

    summary = (df.groupby(["donors", "method"])
                 .agg(fdp_mean=("fdp", "mean"),
                      fdp_sd=("fdp", "std"),
                      power_mean=("power", "mean"),
                      power_sd=("power", "std"),
                      n_called_mean=("n_called", "mean"),
                      tp_mean=("tp", "mean"),
                      fp_mean=("fp", "mean"),
                      n_tested=("n_tested", "mean"))
                 .reset_index())
    summary.to_csv(RESULTS / "sim_nbglm_summary.csv", index=False)

    print(f"\n{time.time() - t0:.0f}s. summary:")
    print(summary.round(3).to_string(index=False))
    return df, summary


if __name__ == "__main__":
    main()
