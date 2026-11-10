"""How badly does the per-cell test fail, and why?

Three experiments on simulated data where the truth is known:

1. Null: no genes change. Count false positives from each method.
2. Power: 50 of 500 genes change. Measure true positives and the observed
   false discovery proportion.
3. Diagnosis: sweep the sample-level variation from zero upwards. If the
   per-cell test is failing because it ignores sample-level structure, it
   should be fine at sample_sd=0 and degrade as sample_sd rises. If it fails
   everywhere, the explanation is wrong.

4. Normalisation: the same pseudobulk test with total-count CPM instead of
   median-of-ratios size factors, to isolate how much damage composition
   bias does on its own.

Run:
    python scripts/simulation_study.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from scpipe import de, simulate  # noqa: E402

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

ALPHA = 0.05
N_REPS = 10


def both_methods(adata):
    sc = de.single_cell_de(adata, "condition", "stim", "ctrl")
    pb = de.pseudobulk_de(adata, "condition", "donor", "stim", "ctrl",
                          paired=True, min_count=5, min_samples=3)
    return sc, pb


def experiment_null():
    print(f"\n--- 1. null: nothing changes, {N_REPS} replicate simulations ---")
    rows = []
    for rep in range(N_REPS):
        adata, _ = simulate.simulate(n_de=0, seed=rep)
        sc, pb = both_methods(adata)
        rows.append({"rep": rep,
                     "per_cell": de.n_significant(sc, ALPHA),
                     "pseudobulk": de.n_significant(pb, ALPHA),
                     "n_genes": len(sc)})
        print(f"  rep {rep + 1:>2}: per-cell {rows[-1]['per_cell']:>4} "
              f"false positives, pseudobulk {rows[-1]['pseudobulk']:>3}",
              flush=True)

    df = pd.DataFrame(rows)
    print(f"\n  per-cell   mean {df['per_cell'].mean():>7.1f} of "
          f"{df['n_genes'].iloc[0]} genes "
          f"({100 * df['per_cell'].mean() / df['n_genes'].iloc[0]:.1f}%)")
    print(f"  pseudobulk mean {df['pseudobulk'].mean():>7.1f} "
          f"({100 * df['pseudobulk'].mean() / df['n_genes'].iloc[0]:.1f}%)")
    print(f"  expected under correct FDR control: about "
          f"{ALPHA * df['n_genes'].iloc[0]:.0f} or fewer")
    df.to_csv(RESULTS / "sim_null.csv", index=False)
    return df


def experiment_power():
    print(f"\n--- 2. power: 50 of 500 genes truly change ---")
    rows = []
    for rep in range(N_REPS):
        adata, truth = simulate.simulate(n_de=50, effect=1.0, seed=100 + rep)
        sc, pb = both_methods(adata)
        for name, res in (("per_cell", sc), ("pseudobulk", pb)):
            m = simulate.confusion(res, truth, ALPHA)
            rows.append({"rep": rep, "method": name, **m})

    df = pd.DataFrame(rows)
    print(f"\n{'method':<12} {'called':>7} {'true pos':>9} {'false pos':>10} "
          f"{'power':>7} {'FDP':>7}")
    for method, g in df.groupby("method"):
        print(f"{method:<12} {g['n_called'].mean():>7.1f} {g['tp'].mean():>9.1f} "
              f"{g['fp'].mean():>10.1f} {g['power'].mean():>7.2f} "
              f"{g['fdp'].mean():>7.2f}")
    print(f"\n  FDP should sit near {ALPHA}. A method with high power and a "
          f"high FDP is not finding more, it is calling more.")
    df.to_csv(RESULTS / "sim_power.csv", index=False)
    return df


def experiment_sample_variation():
    print("\n--- 3. diagnosis: sweep sample-level (donor x condition) variation ---")
    print("  if the problem is unmodelled sample-level structure, per-cell")
    print("  should be fine at 0 and degrade as it rises\n")
    print(f"{'sample_sd':>10} {'per-cell FP':>12} {'pseudobulk FP':>14}")

    rows = []
    for sample_sd in (0.0, 0.1, 0.2, 0.3, 0.5):
        sc_fp, pb_fp = [], []
        for rep in range(5):
            adata, _ = simulate.simulate(n_de=0, sample_sd=sample_sd,
                                         seed=200 + rep)
            sc, pb = both_methods(adata)
            sc_fp.append(de.n_significant(sc, ALPHA))
            pb_fp.append(de.n_significant(pb, ALPHA))
        rows.append({"sample_sd": sample_sd,
                     "per_cell": np.mean(sc_fp),
                     "pseudobulk": np.mean(pb_fp)})
        print(f"{sample_sd:>10.2f} {np.mean(sc_fp):>12.1f} "
              f"{np.mean(pb_fp):>14.1f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "sim_sample_sweep.csv", index=False)
    return df


def experiment_normalisation():
    """Composition bias, and why it gets worse as the study gets bigger.

    Normalising by total count lets a handful of strongly changed genes
    shift every other gene's apparent level. That shift is a bias: it is the
    same in every sample, so averaging over more donors does not reduce it.
    Sampling noise does shrink with more donors. So the ratio of bias to
    noise grows, and a normalisation choice that looks harmless on 8 donors
    becomes the dominant error on 40.
    """
    print("\n--- 4. normalisation: composition bias vs study size ---")
    print("  50 of 500 genes change. Bias does not average out; noise does.\n")
    print(f"{'donors':>7} {'normalisation':<15} {'called':>7} {'true pos':>9} "
          f"{'false pos':>10} {'FDP':>7}")

    rows = []
    for n_donors in (8, 20, 40):
        for norm in ("median_ratio", "cpm"):
            agg = []
            for rep in range(6):
                adata, truth = simulate.simulate(n_de=50, effect=1.0,
                                                 n_donors=n_donors,
                                                 seed=300 + rep)
                res = de.pseudobulk_de(adata, "condition", "donor", "stim",
                                       "ctrl", paired=True, min_count=5,
                                       min_samples=3, normalisation=norm)
                agg.append(simulate.confusion(res, truth, ALPHA))
            d = pd.DataFrame(agg)
            rows.append({"n_donors": n_donors, "normalisation": norm,
                         **d.mean().to_dict()})
            print(f"{n_donors:>7} {norm:<15} {d['n_called'].mean():>7.1f} "
                  f"{d['tp'].mean():>9.1f} {d['fp'].mean():>10.1f} "
                  f"{d['fdp'].mean():>7.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "sim_normalisation.csv", index=False)
    return df


if __name__ == "__main__":
    experiment_null()
    experiment_power()
    experiment_sample_variation()
    experiment_normalisation()
    print(f"\nresults written to {RESULTS}")
