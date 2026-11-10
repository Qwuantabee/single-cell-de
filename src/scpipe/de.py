"""Differential expression, done two ways so they can be compared.

`single_cell_de` is what most tutorials do: treat every cell as an
independent sample and run a rank test per gene.

`pseudobulk_de` aggregates to one profile per donor per condition first and
tests across donors.

Both return the same shape of result so they can be put side by side. The
null experiment in scripts/false_positive_test.py is what decides between
them, by running both on data where the answer is known to be "nothing".
"""

import numpy as np
import pandas as pd
from scipy import stats

from .pseudobulk import aggregate, cpm, filter_genes, normalise


def _bh(pvals):
    """Benjamini-Hochberg FDR. Returns q-values in the original order."""
    p = np.asarray(pvals, dtype=float)
    ok = np.isfinite(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    q[ok] = out
    return q


def single_cell_de(adata, condition_key, group_a, group_b, layer=None):
    """Wilcoxon rank-sum per gene, every cell treated as a replicate.

    This is the approach whose false positive rate the null experiment
    measures. It is included to be compared against, not to be used.
    """
    import scipy.sparse as sp

    X = adata.layers[layer] if layer else adata.X
    labels = adata.obs[condition_key].astype(str).to_numpy()
    ia = np.flatnonzero(labels == group_a)
    ib = np.flatnonzero(labels == group_b)
    if len(ia) == 0 or len(ib) == 0:
        raise ValueError(f"one of {group_a!r}/{group_b!r} has no cells")

    Xa = X[ia]
    Xb = X[ib]
    if sp.issparse(Xa):
        Xa, Xb = Xa.toarray(), Xb.toarray()

    stat, pvals = stats.mannwhitneyu(Xa, Xb, axis=0, alternative="two-sided")
    mean_a = Xa.mean(axis=0)
    mean_b = Xb.mean(axis=0)

    return pd.DataFrame({
        "gene": adata.var_names,
        "mean_a": mean_a,
        "mean_b": mean_b,
        "log2fc": np.log2((mean_a + 1e-9) / (mean_b + 1e-9)),
        "pval": pvals,
        "qval": _bh(pvals),
        "n_a": len(ia),
        "n_b": len(ib),
    }).set_index("gene")


def pseudobulk_de(adata, condition_key, sample_key, group_a, group_b,
                  layer=None, min_cells=10, paired=False,
                  min_count=10, min_samples=3, normalisation="median_ratio"):
    """Aggregate per donor, then test across donors.

    paired=True uses a paired t-test, which is right when every donor
    contributes to both conditions -- as in a stimulate/don't-stimulate
    design on cells from the same person. Pairing removes between-donor
    variation, which is usually the largest source of variance in the whole
    experiment, so it is worth a lot of statistical power.

    normalisation: "median_ratio" (default) or "cpm". Use median_ratio.
    Total-count CPM is vulnerable to composition bias and the simulation
    study measures exactly how much damage that does.
    """
    mat, meta = aggregate(adata, sample_key=sample_key,
                          group_key=condition_key, layer=layer,
                          min_cells=min_cells)
    mat = filter_genes(mat, min_count=min_count, min_samples=min_samples)
    if normalisation == "median_ratio":
        logcpm = normalise(mat, log=True)
    elif normalisation == "cpm":
        logcpm = cpm(mat, log=True)
    else:
        raise ValueError(f"unknown normalisation {normalisation!r}")

    meta = meta.loc[logcpm.index]
    is_a = meta[condition_key] == group_a
    is_b = meta[condition_key] == group_b

    A = logcpm.loc[is_a]
    B = logcpm.loc[is_b]
    if len(A) < 2 or len(B) < 2:
        raise ValueError(
            f"need at least 2 donors per condition, got {len(A)} and {len(B)}")

    if paired:
        sa = meta.loc[is_a, sample_key].to_numpy()
        sb = meta.loc[is_b, sample_key].to_numpy()
        shared = [s for s in sa if s in set(sb)]
        if len(shared) < 2:
            raise ValueError("paired=True but donors do not overlap")
        A = A.loc[meta.loc[is_a, sample_key].isin(shared).to_numpy()]
        B = B.loc[meta.loc[is_b, sample_key].isin(shared).to_numpy()]
        A = A.iloc[np.argsort(meta.loc[A.index, sample_key].to_numpy())]
        B = B.iloc[np.argsort(meta.loc[B.index, sample_key].to_numpy())]
        stat, pvals = stats.ttest_rel(A.to_numpy(), B.to_numpy(), axis=0)
    else:
        stat, pvals = stats.ttest_ind(A.to_numpy(), B.to_numpy(), axis=0,
                                      equal_var=False)

    mean_a = A.mean(axis=0).to_numpy()
    mean_b = B.mean(axis=0).to_numpy()

    return pd.DataFrame({
        "gene": logcpm.columns,
        "mean_a": mean_a,
        "mean_b": mean_b,
        "log2fc": mean_a - mean_b,       # already log2 CPM
        "pval": pvals,
        "qval": _bh(pvals),
        "n_a": len(A),
        "n_b": len(B),
    }).set_index("gene")


def n_significant(result, alpha=0.05):
    return int((result["qval"] < alpha).sum())
