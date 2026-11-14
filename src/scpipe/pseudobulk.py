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
    # .to_numpy() can hand back a read-only view (pandas 3.0 copy-on-write,
    # and any DataFrame built on a read-only array), so copy before editing.
    lib = np.array(matrix.sum(axis=1).to_numpy(), dtype=float)[:, None]
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
    sf = np.array(size_factors(matrix).to_numpy(), dtype=float)[:, None]
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


def cell_type_composition(adata, cell_type_key="cell_type",
                          sample_key="donor", condition_key="condition"):
    """Proportion of each cell type within each (sample, condition) sample.

    This exists to be looked at *before* any per-gene test, because the two
    questions are different and only one of them is about expression.

    If interferon-beta stimulation, or a disease, changes the mix of cells
    that are recovered -- more monocytes, fewer B cells, a new activated
    subset -- then a test that aggregates all cells of a donor together
    compares a monocyte-heavy sample to a B-cell-heavy one. Every monocyte
    marker comes out "up" and every B-cell marker "down", with p-values that
    are perfectly valid for the question "do these two cell mixtures differ"
    and completely wrong for the question "did any gene's expression
    change". The false positives are not noise; they are the composition
    difference, measured accurately and reported as biology.

    Returns (proportions, shift):
        proportions: samples x cell types, each row summing to 1
        shift: one row per cell type with the mean proportion in each
            condition and the log2 ratio between them. A shift of +1 means
            that type is twice as common in the first condition.
    """
    obs = adata.obs
    tab = pd.crosstab([obs[sample_key], obs[condition_key]],
                      obs[cell_type_key])
    proportions = tab.div(tab.sum(axis=1), axis=0)

    conds = proportions.index.get_level_values(1).unique().tolist()
    if len(conds) != 2:
        raise ValueError(
            f"expected 2 conditions for a shift estimate, found {conds}")
    a, b = conds[0], conds[1]
    mask = proportions.index.get_level_values(1)
    mean_a = proportions[mask == a].mean()
    mean_b = proportions[mask == b].mean()

    eps = 1e-12
    shift = pd.DataFrame({
        "mean_prop_a": mean_a,
        "mean_prop_b": mean_b,
        "log2_ratio": np.log2((mean_a + eps) / (mean_b + eps)),
        "n_a": tab[mask == a].sum(),
        "n_b": tab[mask == b].sum(),
    })
    shift.index.name = cell_type_key
    return proportions, shift


def pseudobulk_by_celltype(adata, cell_type_key="cell_type",
                           sample_key="donor", condition_key="condition",
                           group_a="stim", group_b="ctrl", min_cells=25,
                           min_donors=4, min_count=10, min_samples=3,
                           normalisation="median_ratio", paired=True):
    """Aggregate within each cell type, then test within each cell type.

    The single all-cells aggregate cannot answer an expression question on
    data whose cell composition differs between conditions -- see
    cell_type_composition for why. Splitting by cell type first removes the
    composition difference from the comparison by construction: within a
    cell type, condition is the only thing that can differ, so a change in
    the *mix* of cells has nowhere to show up except in the composition
    table, where it belongs.

    Cell types with fewer than `min_donors` donors contributing
    `min_cells` cells to *both* conditions are skipped rather than tested on
    a handful of donors, and the reason is returned instead of dropped
    silently.

    Multiple testing: qval is Benjamini-Hochberg within each cell type
    (the usual per-contrast correction). qval_global is BH across every
    cell-type x gene test together, which is the honest correction if the
    claim being made is "this gene changed somewhere".

    Returns (results, composition, skipped):
        results: long frame indexed by (cell_type, gene)
        composition: the per-sample proportions that motivated the split
        skipped: {cell_type: reason}
    """
    from .de import _bh, pseudobulk_de

    composition, _ = cell_type_composition(
        adata, cell_type_key=cell_type_key, sample_key=sample_key,
        condition_key=condition_key)

    order = adata.obs[cell_type_key].value_counts().index.tolist()
    frames, skipped = [], {}

    for ct in order:
        sub = adata[adata.obs[cell_type_key] == ct].copy()
        tab = pd.crosstab(sub.obs[sample_key], sub.obs[condition_key])
        if group_a not in tab.columns or group_b not in tab.columns:
            skipped[ct] = "missing a condition"
            continue
        usable = tab[(tab[group_a] >= min_cells)
                     & (tab[group_b] >= min_cells)].index.tolist()
        if len(usable) < min_donors:
            skipped[ct] = (f"{len(usable)} donors with >= {min_cells} cells "
                           f"in both conditions (< {min_donors})")
            continue
        sub = sub[sub.obs[sample_key].isin(usable)].copy()
        try:
            res = pseudobulk_de(
                sub, condition_key, sample_key, group_a, group_b,
                min_cells=min_cells, paired=paired, min_count=min_count,
                min_samples=min_samples, normalisation=normalisation)
        except ValueError as exc:
            skipped[ct] = str(exc)
            continue
        res.insert(0, "cell_type", ct)
        res.insert(1, "n_donors", len(usable))
        frames.append(res)

    if not frames:
        raise ValueError(
            "no cell type had enough donors to test; skipped: "
            f"{skipped}")

    results = pd.concat(frames)
    results.index.name = "gene"
    results = results.set_index("cell_type", append=True)
    results = results.reorder_levels(["cell_type", "gene"]).sort_index()
    results["qval_global"] = _bh(results["pval"].to_numpy())
    return results, composition, skipped


def significant_by_celltype(results, alpha=0.05, column="qval"):
    """How many genes each cell type calls, using one q-value column."""
    called = results[results[column] < alpha]
    per_type = called.groupby(level="cell_type").size()
    order = results.index.get_level_values("cell_type").unique()
    return per_type.reindex(order).fillna(0).astype(int)
