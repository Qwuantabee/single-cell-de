"""Pseudobulk aggregation for differential expression.

This module exists because of a specific, well documented mistake.

The obvious way to find genes that change between two conditions in
single-cell data is to take all the cells of one type, split them by
condition, and run a Wilcoxon test per gene. Scanpy makes that a one-liner
and it is what most tutorials show. It is also wrong in a way that matters:
cells from the same individual are not independent observations. Treating
10,000 cells from 8 donors as 10,000 independent samples inflates the
effective sample size by three orders of magnitude, and the p-values that
come out are not interpretable.

Squair et al. (2021, Nat Commun 12:5692) showed that this produces large
numbers of false positives, and that the genes it falsely calls are
systematically the highly expressed ones. Aggregating counts per
donor-per-condition first ("pseudobulk") and then testing across donors
gives the right number of degrees of freedom, because the replicate is the
donor, which is what was actually randomised.

So the comparison this module supports is not just methodological trivia:
it changes which genes you would take to a bench scientist.
"""

import numpy as np
import pandas as pd
import pandas as pd


def aggregate(adata, sample_key, group_key=None, layer=None, min_cells=10):
    """Sum raw counts within each (sample, group) combination.

    Summing rather than averaging is deliberate: the downstream count models
    (and the library-size normalisation that follows) expect counts, and
    summing preserves the relationship between depth and variance that makes
    those models work.

    sample_key: the biological replicate. Almost always donor or animal, not
        the 10x lane, and getting this wrong is the whole problem.
    group_key: optional, usually cell type, so each cell type is aggregated
        separately.
    min_cells: combinations with fewer cells than this are dropped. A
        pseudobulk profile built from 3 cells is noise wearing the costume
        of a replicate.
    """
    import scipy.sparse as sp

    X = adata.layers[layer] if layer else adata.X
    obs = adata.obs

    keys = [sample_key] + ([group_key] if group_key else [])
    combos = obs[keys].astype(str).agg("||".join, axis=1)

    rows, index = [], []
    counts = []
    for combo, idx in combos.groupby(combos).groups.items():
        mask = obs.index.get_indexer(idx)
        if len(mask) < min_cells:
            continue
        sub = X[mask]
        total = np.asarray(sub.sum(axis=0)).ravel() if sp.issparse(sub) \
            else np.asarray(sub.sum(axis=0)).ravel()
        rows.append(total)
        index.append(combo)
        counts.append(len(mask))

    if not rows:
        raise ValueError(
            f"no {sample_key} x {group_key} combination had >= {min_cells} cells")

    mat = pd.DataFrame(np.vstack(rows), index=index, columns=adata.var_names)
    meta = pd.DataFrame(
        [c.split("||") for c in index], index=index, columns=keys)
    meta["n_cells"] = counts
    return mat, meta


def cpm(matrix, log=True, prior_count=1.0):
    """Counts per million by total count, optionally log2.

    Kept because it is the obvious thing to do, and because the simulation
    study uses it to demonstrate why it is not good enough. Prefer
    size_factors() + normalise() for anything being tested.

    prior_count damps the variance of genes with near-zero counts, which
    otherwise dominate any fold-change ranking.
    """
    lib = matrix.sum(axis=1).to_numpy()[:, None]
    lib[lib == 0] = 1
    out = matrix / lib * 1e6
    if log:
        out = np.log2(out + prior_count)
    return out


def size_factors(matrix):
    """Median-of-ratios size factors, as used by DESeq2.

    Normalising by total count assumes the total amount of RNA is the same
    in every sample. It is not, and more importantly it does not have to be:
    if a handful of genes go up strongly, they take up more of the library
    and every other gene's share goes down even though nothing happened to
    it. That is composition bias, and on simulated data where only 50 of 500
    genes truly change it shifted every null gene by -0.07 log2 CPM --
    consistently enough across donors that a paired test called 163 of them
    significant.

    The median of ratios is robust to it: a minority of genes can move as
    much as they like without shifting the median.
    """
    counts = matrix.to_numpy(dtype=float)
    with np.errstate(divide="ignore"):
        log_counts = np.log(counts)
    # reference is the geometric mean across samples, over genes present
    # in every sample
    usable = np.all(counts > 0, axis=0)
    if usable.sum() < 10:
        # Too sparse for a geometric-mean reference: no gene is detected in
        # every sample, so there is nothing robust to anchor on. Fall back to
        # library size and accept the composition bias, because the
        # alternative is a reference computed from a handful of genes, which
        # is worse. An empty sample gets a factor of 1 rather than a
        # division by zero.
        lib = counts.sum(axis=1)
        lib = np.where(lib > 0, lib, 1.0)
        return pd.Series(lib / np.exp(np.mean(np.log(lib))),
                         index=matrix.index)

    log_ref = log_counts[:, usable].mean(axis=0)
    ratios = log_counts[:, usable] - log_ref
    factors = np.exp(np.median(ratios, axis=1))
    return pd.Series(factors, index=matrix.index)


def normalise(matrix, log=True, prior_count=1.0):
    """Divide by median-of-ratios size factors, then optionally log2."""
    sf = size_factors(matrix).to_numpy()[:, None]
    sf[sf <= 0] = 1.0
    out = matrix / sf
    if log:
        out = np.log2(out + prior_count)
    return out


def filter_genes(matrix, min_count=10, min_samples=3):
    """Drop genes that are not detected often enough to be testable.

    Filtering before testing rather than after is important: it reduces the
    multiple-testing burden honestly, whereas filtering on the result would
    be selecting on the outcome.
    """
    detected = (matrix >= min_count).sum(axis=0)
    keep = detected >= min_samples
    return matrix.loc[:, keep]
