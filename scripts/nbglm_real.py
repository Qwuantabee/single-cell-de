"""The NB GLM against the paired t-test, on the real Kang 2018 counts.

`scripts/nbglm_study.py` measured the GLM's advantage on simulated data. The
simulator draws from a negative binomial with a per-donor and per-sample
random effect, which is the GLM's own model, so the 97%-vs-68% power gap at
6 donors is a best case by construction. The README says so and then adds a
sentence that had never been tested: that on real data the gap is unlikely to
be larger than this, and could be smaller.

This script tests that sentence. Both methods get the same pseudobulk matrix
from the same real cells, the same gene filter, the same BH correction and the
same 8 paired donors (ctrl vs 6 h interferon-beta). Cell-type stratification
matters here for the reason the rest of the repo gives: a composition shift
between conditions becomes an expression change when cells are pooled across
types, so CD4 T cells is the primary comparison and every other type with
enough cells is reported beside it.

Four questions, in the order they should be read:

  1. How many genes does each method call at BH q < 0.05, and how much do the
     two lists overlap? Jaccard is reported because "both find a lot" and
     "both find the same lot" are different statements.
  2. Do the effect sizes agree? Pearson and Spearman on log2 fold change,
     pairing only genes both methods tested.
  3. Are the 12 canonical interferon-stimulated genes found, and at what
     rank? A method that misses ISG15 on this dataset is broken regardless of
     its call count.
  4. What does each method do when there is provably nothing to find?

The null in (4) swaps condition labels *within donors* via
`scpipe.permutation.swap_cells`, so every cell keeps its donor and the
donor-to-donor spread is untouched. Two families of draws are run, because
they mean different things:

  fair coin   -- each donor flipped with probability 0.5. This is the
                 module's native behaviour. Only k = n/2 flips leave the real
                 condition effect absent from *every* donor; at k = 6 of 8 the
                 two unflipped donors still carry the genuine interferon
                 response, so a rejection there may be sensitivity rather than
                 miscalibration.
  balanced    -- rejection-sampled seeds until exactly n/2 donors are flipped.
                 Still `swap_cells` doing the swapping; only the seed is
                 redrawn. This is the strictly-null family.

A calibrated procedure returns a mean of about 0.05 rejections per family on
all-null data: BH bounds the false-discovery *proportion*, so with no true
positives any rejection at all is a failure (README section 1). A call count
alone cannot distinguish a sensitive test from an anti-conservative one; the
null can.

Run:
    python scripts/nbglm_real.py
    python scripts/nbglm_real.py --null-on "CD14" --n-perm 20
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from scpipe import data, de, nbglm, permutation, pseudobulk  # noqa: E402

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

ALPHA = 0.05
MIN_CELLS = 25
N_PERM = 10
PRIMARY = "CD4 T cells"
ISG = data.ISG_MARKERS

# Only genes with no measured response can be treated as null by a permutation
# on this dataset. The interferon programme genuinely moves, so a permuted call
# on an ISG is not a false positive.
NULL_FC_CUTOFF = 0.5


# ---------------------------------------------------------------- helpers


def called_genes(result, alpha=ALPHA, column="qval"):
    """Genes called significant, as a set. NaN q-values are not calls."""
    if result is None or len(result) == 0:
        return set()
    q = result[column]
    return set(result.index[q.notna() & (q < alpha)])


def overlap_stats(pb, nb, alpha=ALPHA):
    """Call counts, overlap and Jaccard for two result frames."""
    a = called_genes(pb, alpha)
    b = called_genes(nb, alpha)
    both = a & b
    union = a | b
    return {
        "n_tested_pb": int(len(pb)),
        "n_tested_nb": int(len(nb)),
        "n_pb": len(a),
        "n_nb": len(b),
        "n_both": len(both),
        "n_pb_only": len(a - b),
        "n_nb_only": len(b - a),
        "jaccard": (len(both) / len(union)) if union else float("nan"),
        "pb_only_genes": sorted(a - b),
        "nb_only_genes": sorted(b - a),
    }


def log2fc_agreement(pb, nb):
    """Pearson and Spearman on log2FC, over genes both methods tested.

    `log2fc` is the log2 CPM difference in both frames (nbglm also carries the
    model coefficient as `log2fc_glm`; the shared column is compared so this
    is a like-for-like check of the summary statistic).
    """
    m = pb[["log2fc"]].join(nb[["log2fc"]], how="inner", lsuffix="_pb",
                            rsuffix="_nb")
    m = m.dropna()
    x = m["log2fc_pb"].to_numpy(dtype=float)
    y = m["log2fc_nb"].to_numpy(dtype=float)
    if len(m) < 3:
        return {"n_genes": int(len(m)), "pearson": float("nan"),
                "spearman": float("nan"), "mean_abs_diff": float("nan")}
    from scipy import stats as st
    flat = np.std(x) == 0 or np.std(y) == 0
    pr = float("nan") if flat else float(st.pearsonr(x, y)[0])
    sp = float("nan") if flat else float(st.spearmanr(x, y)[0])
    return {"n_genes": int(len(m)), "pearson": pr, "spearman": sp,
            "mean_abs_diff": float(np.mean(np.abs(x - y)))}


def rank_by_q(result, column="qval"):
    """Gene -> rank by q-value, 1 = most significant. NaN q-values last."""
    q = result[column]
    order = q.sort_values(kind="mergesort").index
    return pd.Series(np.arange(1, len(order) + 1), index=order)


def isg_table(pb, nb, markers=ISG, alpha=ALPHA):
    """Per-marker q-value, rank and call, for both methods."""
    rp = rank_by_q(pb)
    rn = rank_by_q(nb)
    cp = called_genes(pb, alpha)
    cn = called_genes(nb, alpha)
    rows = []
    for g in markers:
        rows.append({
            "gene": g,
            "q_pb": (float(pb["qval"].get(g, np.nan))
                     if g in pb.index else float("nan")),
            "rank_pb": int(rp[g]) if g in rp.index else -1,
            "called_pb": g in cp,
            "q_nb": (float(nb["qval"].get(g, np.nan))
                     if g in nb.index else float("nan")),
            "rank_nb": int(rn[g]) if g in rn.index else -1,
            "called_nb": g in cn,
        })
    return pd.DataFrame(rows).set_index("gene")


def isg_median_rank(tab):
    """Median rank and call count for the markers present in both rankings."""
    both = tab[(tab["rank_pb"] > 0) & (tab["rank_nb"] > 0)]
    return {
        "n_isg": int(len(both)),
        "median_rank_pb": float(both["rank_pb"].median()) if len(both)
        else float("nan"),
        "median_rank_nb": float(both["rank_nb"].median()) if len(both)
        else float("nan"),
        "n_called_pb": int(tab["called_pb"].sum()),
        "n_called_nb": int(tab["called_nb"].sum()),
    }


def null_gene_set(pb, nb, markers=ISG, cutoff=NULL_FC_CUTOFF):
    """Genes a permutation can treat as null: tested by both, small effect.

    Excluding the interferon programme makes the null stricter, not looser:
    any false positive left sits on a gene with no measured response.
    """
    common = pb.index.intersection(nb.index)
    small = common[pb.loc[common, "log2fc"].abs() < cutoff]
    return [g for g in small if g not in set(markers)]


def _n_donors_flipped(before, after):
    """How many donors had their condition labels exchanged."""
    db = np.asarray(before.obs["donor"].astype(str))
    cb = np.asarray(before.obs["condition"].astype(str))
    ca = np.asarray(after.obs["condition"].astype(str))
    flipped = 0
    for unit in sorted(set(db)):
        m = db == unit
        if not np.array_equal(cb[m], ca[m]):
            flipped += 1
    return flipped


def check_cell_swap(donors_before, donors_after, conditions_before,
                    conditions_after):
    """Two properties a valid cell-level swap must have, checked directly.

    `permutation.donor_structure_preserved` asserts each unit's multiset of
    labels is unchanged. That is the right invariant for pseudobulk rows (one
    row per donor per condition) but not for cells: a donor with 353 ctrl and
    467 stim cells legitimately becomes 467 ctrl and 353 stim after a swap, so
    the multiset changes while the design does not. The cell-level invariant
    is: no cell changes donor, and each donor's per-condition counts are
    exchanged exactly.
    """
    db = np.asarray(donors_before).astype(str)
    da = np.asarray(donors_after).astype(str)
    cb = np.asarray(conditions_before).astype(str)
    ca = np.asarray(conditions_after).astype(str)

    problems = []
    if not np.array_equal(db, da):
        moved = int((db != da).sum())
        problems.append(f"{moved} cells changed donor")
    for unit in sorted(set(db)):
        m = db == unit
        before = {lab: int((cb[m] == lab).sum()) for lab in set(cb[m])}
        after = {lab: int((ca[m] == lab).sum()) for lab in set(ca[m])}
        if set(before) != set(after):
            problems.append(f"{unit}: levels {sorted(before)} -> "
                            f"{sorted(after)}")
            break
        if before == after:
            continue                    # this donor was not swapped
        if len(before) == 2:
            (k1, v1), (k2, v2) = list(before.items())
            if {k1: v2, k2: v1} != after:
                problems.append(f"{unit}: counts {before} -> {after} "
                                f"is not an exchange")
                break
        else:
            problems.append(f"{unit}: counts {before} -> {after}")
            break
    return (not problems), problems


def _sample_level_rows(ad):
    """Per-(donor, condition) rows, for the sample-level invariant."""
    donors = np.asarray(ad.obs["donor"].astype(str))
    conds = np.asarray(ad.obs["condition"].astype(str))
    keys = sorted({(d, c) for d, c in zip(donors, conds)})
    return (np.array([k[0] for k in keys], dtype=object),
            np.array([k[1] for k in keys], dtype=object))


def draw_balanced_plan(obs, rng, n_donors, group_a="stim", group_b="ctrl",
                       max_draws=10000):
    """A swap plan with exactly half the donors flipped.

    Independent fair coins give k ~ Binomial(n, 0.5); only k = n/2 leaves the
    real effect absent from every donor. The plan still comes from
    `scpipe.permutation.swap_cells` -- no new shuffler -- and the seed is
    redrawn until that module's own plan has the right k.
    """
    for _ in range(max_draws):
        labels, plan = permutation.swap_cells(
            obs["donor"].to_numpy(), obs["condition"].to_numpy(), rng,
            group_a=group_a, group_b=group_b)
        if sum(bool(v) for v in plan.values()) == n_donors // 2:
            return labels, plan
    raise RuntimeError("could not draw a balanced swap plan")


def permuted_labels(obs, rng, group_a="stim", group_b="ctrl"):
    """Condition labels swapped within each donor, cells never moved."""
    labels, plan = permutation.swap_cells(obs["donor"].to_numpy(),
                                          obs["condition"].to_numpy(), rng,
                                          group_a=group_a, group_b=group_b)
    return labels, plan


def permute_adata(adata, rng, balanced=False):
    """Copy of adata with only the condition column permuted.

    Only labels move; the expression matrix is shared, not rewritten.
    """
    import anndata as ad

    new = ad.AnnData(X=adata.X, obs=adata.obs.copy(), var=adata.var.copy())
    if balanced:
        labels, _ = draw_balanced_plan(new.obs, rng,
                                       new.obs["donor"].nunique())
    else:
        labels, _ = permuted_labels(new.obs, rng)
    new.obs["condition"] = labels
    return new


def null_summary(records):
    """False positives per permutation, overall and on the small-effect set."""
    rows = []
    for method in ("pseudobulk", "nb_glm"):
        vals = [r[f"fp_{method}"] for r in records]
        sub = [r[f"fp_{method}_null"] for r in records]
        bal = [r for r in records if r["balanced"]]
        rows.append({
            "method": method,
            "n_perm": len(records),
            "mean_fp": float(np.mean(vals)),
            "median_fp": float(np.median(vals)),
            "max_fp": int(np.max(vals)),
            "min_fp": int(np.min(vals)),
            "mean_fp_small_effect": float(np.mean(sub)),
            "max_fp_small_effect": int(np.max(sub)),
            "n_balanced": len(bal),
            "mean_fp_balanced": (float(np.mean([r[f"fp_{method}"]
                                                for r in bal]))
                                 if bal else float("nan")),
            "max_fp_balanced": (int(np.max([r[f"fp_{method}"] for r in bal]))
                                if bal else -1),
            "n_tested": int(records[0][f"n_{method}"]),
            # BH on an all-null family calls about 0.05 genes per family, not
            # 5% of them: it bounds the false-discovery PROPORTION, and with
            # no true positives any proportion above zero is a failure.
            "n_expected_calibrated": 0.05,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- drivers


def run_one(sub, label):
    """Both methods on one cell-type subset. Returns frames and timings."""
    t0 = time.perf_counter()
    pb = de.pseudobulk_de(sub, "condition", "donor", "stim", "ctrl",
                          paired=True, min_cells=MIN_CELLS)
    t_pb = time.perf_counter() - t0

    t0 = time.perf_counter()
    nb = nbglm.nbglm_de(sub, "condition", "donor", "stim", "ctrl",
                        paired=True, min_cells=MIN_CELLS)
    t_nb = time.perf_counter() - t0

    # A gene that carries counts but still has p = 1 hit the module's guard,
    # i.e. the Wald fit failed. Those are convergence failures, not nulls.
    mat, _ = pseudobulk.aggregate(sub, sample_key="donor",
                                  group_key="condition", min_cells=MIN_CELLS)
    has_counts = set(mat.columns[(mat.sum(axis=0) > 0).to_numpy()])
    stuck = nb[(nb["pval"] >= 1.0) & nb.index.isin(has_counts)].index.tolist()

    return {"label": label, "pb": pb, "nb": nb, "t_pb": t_pb, "t_nb": t_nb,
            "failed": stuck}


def usable_subset(adata, ct):
    """Cells of one type, restricted to donors with both conditions covered."""
    sub = adata[adata.obs["cell_type"] == ct].copy()
    tab = pd.crosstab(sub.obs["donor"], sub.obs["condition"])
    ok = tab[(tab >= MIN_CELLS).all(axis=1)].index.tolist()
    if len(ok) < 4:
        return None, len(ok)
    return sub[sub.obs["donor"].isin(ok)].copy(), len(ok)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--null-on", default="CD4", type=str,
                    help="cell type substring for the permutation null")
    ap.add_argument("--n-perm", default=N_PERM, type=int,
                    help="draws per family (fair-coin and balanced each)")
    args = ap.parse_args(argv)

    adata = data.load_kang()
    data.describe(adata)
    cell_types = (adata.obs["cell_type"].value_counts()
                  .loc[lambda s: s > 500].index.tolist())
    print(f"\ncomparing on {len(cell_types)} cell types: {cell_types}")

    frames, summary, pairs, isg_rows = {}, [], {}, {}
    for ct in cell_types:
        sub, n_usable = usable_subset(adata, ct)
        if sub is None:
            print(f"\n{ct}: only {n_usable} usable donors, skipping")
            continue
        r = run_one(sub, ct)
        frames[ct] = r
        ov = overlap_stats(r["pb"], r["nb"])
        agr = log2fc_agreement(r["pb"], r["nb"])
        it = isg_median_rank(isg_table(r["pb"], r["nb"]))
        n_donors = sub.obs["donor"].nunique()

        print(f"\n{ct}  ({sub.n_obs:,} cells, {n_donors} donors, "
              f"{n_usable} usable)")
        print(f"  pseudobulk  {ov['n_pb']:>5} / {ov['n_tested_pb']:>5} "
              f"significant at q<{ALPHA}   [{r['t_pb']:.1f}s]")
        print(f"  nb_glm      {ov['n_nb']:>5} / {ov['n_tested_nb']:>5} "
              f"significant at q<{ALPHA}   [{r['t_nb']:.1f}s]")
        print(f"  overlap {ov['n_both']}, Jaccard {ov['jaccard']:.3f}, "
              f"pb-only {ov['n_pb_only']}, glm-only {ov['n_nb_only']}")
        print(f"  log2FC: Pearson {agr['pearson']:.3f}, "
              f"Spearman {agr['spearman']:.3f} on {agr['n_genes']} genes")
        print(f"  ISGs called {it['n_called_pb']}/12 vs "
              f"{it['n_called_nb']}/12, median rank "
              f"{it['median_rank_pb']:.0f} vs {it['median_rank_nb']:.0f}")
        if r["failed"]:
            print(f"  GLM did not fit {len(r['failed'])} genes that carry "
                  f"counts (first: {r['failed'][:5]})")

        summary.append({
            "cell_type": ct, "cells": int(sub.n_obs),
            "donors": int(n_donors),
            **{k: v for k, v in ov.items() if not k.endswith("_genes")},
            **agr, **it,
            "t_pb_s": round(r["t_pb"], 2), "t_nb_s": round(r["t_nb"], 2),
            "glm_failures": len(r["failed"]),
            "glm_failed_genes": ";".join(r["failed"]),
        })
        pairs[ct] = r["pb"][["log2fc", "qval"]].join(
            r["nb"][["log2fc", "qval", "log2fc_glm"]], how="outer",
            lsuffix="_pb", rsuffix="_nb")
        isg_rows[ct] = isg_table(r["pb"], r["nb"])

    if not frames:
        raise SystemExit("no cell type had enough donors")

    # ---- permutation null on the primary cell type ----------------------
    target = [ct for ct in frames if args.null_on.lower() in ct.lower()]
    if not target:
        raise SystemExit(f"no cell type matches {args.null_on!r}")
    primary = target[0]
    sub, _ = usable_subset(adata, primary)

    small_set = set(null_gene_set(frames[primary]["pb"], frames[primary]["nb"]))
    n_donor_total = sub.obs["donor"].nunique()
    print(f"\npermutation null on {primary}: {args.n_perm} fair-coin draws + "
          f"{args.n_perm} balanced draws (k={n_donor_total // 2}), "
          f"{len(small_set)} small-effect genes tracked")

    records = []
    for balanced in (False, True):
        for i in range(args.n_perm):
            rng = np.random.default_rng(1000 + i)
            perm = permute_adata(sub, rng, balanced=balanced)
            flipped = _n_donors_flipped(sub, perm)
            ok, problems = check_cell_swap(
                sub.obs["donor"].to_numpy(), perm.obs["donor"].to_numpy(),
                sub.obs["condition"].to_numpy(),
                perm.obs["condition"].to_numpy())
            if not ok:
                raise SystemExit(f"invalid cell swap: {problems}")
            if balanced and flipped != n_donor_total // 2:
                raise SystemExit(f"balanced draw flipped {flipped} of "
                                 f"{n_donor_total} donors")
            # the sample-level invariant the module was written for
            b_d, b_c = _sample_level_rows(sub)
            a_d, a_c = _sample_level_rows(perm)
            ok_s, problems_s = permutation.donor_structure_preserved(
                b_d, a_d, b_c, a_c)
            if not ok_s:
                raise SystemExit(f"donor structure broken: {problems_s}")

            p_pb = de.pseudobulk_de(perm, "condition", "donor", "stim",
                                    "ctrl", paired=True, min_cells=MIN_CELLS)
            p_nb = nbglm.nbglm_de(perm, "condition", "donor", "stim",
                                  "ctrl", paired=True, min_cells=MIN_CELLS)
            c_pb = called_genes(p_pb)
            c_nb = called_genes(p_nb)
            rec = {
                "perm": i,
                "balanced": balanced,
                "n_donors_flipped": flipped,
                "fp_pseudobulk": len(c_pb),
                "fp_nb_glm": len(c_nb),
                "fp_pseudobulk_null": len(c_pb & small_set),
                "fp_nb_glm_null": len(c_nb & small_set),
                "n_pseudobulk": int(len(p_pb)),
                "n_nb_glm": int(len(p_nb)),
                "n_small_set": len(small_set),
            }
            records.append(rec)
            tag = "balanced" if balanced else "coin    "
            print(f"  {tag} {i + 1:>2} (k={flipped}): pseudobulk "
                  f"{rec['fp_pseudobulk']} ({rec['fp_pseudobulk_null']} "
                  f"small-effect), nb_glm {rec['fp_nb_glm']} "
                  f"({rec['fp_nb_glm_null']})", flush=True)

    null = null_summary(records)
    print("\nnull false positives per permutation:")
    print(null.to_string(index=False))
    bal_only = null[["method", "n_balanced", "mean_fp_balanced",
                     "max_fp_balanced"]]
    print("\nbalanced draws only (k = n/2, the strictly-null plans):")
    print(bal_only.to_string(index=False))

    # ---- persist --------------------------------------------------------
    pd.DataFrame(summary).to_csv(RESULTS / "nbglm_real_summary.csv",
                                 index=False)
    pd.concat(pairs.values(), keys=pairs.keys(),
              names=["cell_type"]).to_csv(RESULTS / "nbglm_real_pairs.csv")
    pd.concat(isg_rows.values(), keys=isg_rows.keys(),
              names=["cell_type"]).to_csv(RESULTS / "nbglm_real_isg.csv")
    pd.DataFrame(records).to_csv(RESULTS / "nbglm_real_null.csv", index=False)

    prim = frames[primary]
    p_ov = overlap_stats(prim["pb"], prim["nb"])
    verdict = {
        "primary_cell_type": primary,
        "primary_n_pb": p_ov["n_pb"],
        "primary_n_nb": p_ov["n_nb"],
        "primary_jaccard": p_ov["jaccard"],
        "primary_log2fc_pearson": log2fc_agreement(prim["pb"],
                                                   prim["nb"])["pearson"],
        "null_mean_fp_pseudobulk": float(
            null.loc[null.method == "pseudobulk", "mean_fp"].iloc[0]),
        "null_mean_fp_nb_glm": float(
            null.loc[null.method == "nb_glm", "mean_fp"].iloc[0]),
        "null_mean_fp_balanced_nb_glm": float(
            null.loc[null.method == "nb_glm", "mean_fp_balanced"].iloc[0]),
        "null_mean_fp_balanced_pseudobulk": float(
            null.loc[null.method == "pseudobulk",
                     "mean_fp_balanced"].iloc[0]),
        "simulated_claim": "97% vs 68% power at 6 donors, from "
                           "simulate.simulate which draws NB data -- the "
                           "GLM's own model",
    }
    with open(RESULTS / "nbglm_real_verdict.json", "w") as fh:
        json.dump(verdict, fh, indent=2)
    print("\n" + json.dumps(verdict, indent=2))
    return verdict


if __name__ == "__main__":
    main()
