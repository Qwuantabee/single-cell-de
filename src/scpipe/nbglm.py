"""Negative-binomial GLM differential expression, the way DESeq2/edgeR do it.

`de.pseudobulk_de` tests log-transformed, size-factor-normalised pseudobulk
counts with a t-test. That is a normal approximation: it assumes the
per-sample values are Gaussian with constant variance and it never looks at
the mean-variance relationship of counts. On a handful of donors the
approximation is the weakest link -- log(x + 1) of a gene with a mean of 3
counts is a lumpy, skewed quantity, and the t-test pretends otherwise.

This module instead models the counts directly:

    y_ij ~ NB(mu_ij, alpha_g)
    log(mu_ij) = log(sf_i) + beta0 + beta_cond * x_i  (+ donor effects)
    var(y_ij) = mu_ij + alpha_g * mu_ij^2

with median-of-ratios size factors as the offset (the same ones
`pseudobulk.size_factors` uses), a per-gene dispersion alpha_g estimated by
method of moments and then shrunk toward a mean-dispersion trend, and a Wald
z-test on the condition coefficient.

Three pieces, each of which exists because skipping it is a known way to get
the wrong answer:

1. Size factors as an *offset*, not a pre-division. Dividing counts by a
   non-integer factor destroys the Poisson/NB variance structure the model
   relies on.
2. Dispersion shrinkage. Per-gene method-of-moments dispersion from 6
   pseudobulk samples is extremely noisy, and when it lands low the Wald
   statistic is inflated. Genes are therefore shrunk toward a trend fitted
   over log mean expression, weighted by how much information the gene
   actually carries. DESeq2 calls the same step "dispersion estimation";
   the weighting here is a simpler version of their empirical Bayes.
3. A paired design keeps the donor term. When every donor contributes to
   both conditions, `~ donor + condition` removes between-donor variance
   from the condition estimate exactly as the paired t-test does.

What comes out is a p-value per gene and BH q-values, so the result can be
dropped straight next to `de.pseudobulk_de` in a comparison.
"""

import warnings

import numpy as np
import pandas as pd

from scipy import stats

from .de import _bh
from .pseudobulk import aggregate, filter_genes, size_factors

DISP_FLOOR = 1e-4

# A fitted condition coefficient larger than this many log2 units (a
# million-fold) is not an estimate: it is a (quasi-)separated gene, where one
# condition has zero counts in every donor and the likelihood has no finite
# maximum. The Wald test is not valid there. The threshold is a flag, not an
# inferential choice: separated genes are reported as their own category and
# nothing is dropped.
SEPARATION_LOG2 = 20.0


def _log2norm_means(counts, groups_a, groups_b, sf):
    """Mean log2 of size-factor-normalised counts per group (fold change).

    Uses the SAME normalisation as `de.pseudobulk_de` (counts / median-of-
    ratios size factor, +1, log2), so `log2fc` of the two methods is a
    like-for-like quantity. An earlier version used total-count CPM here, which
    differs from the pseudobulk arm by a per-sample composition offset: the two
    log2FC columns then disagreed by a mean of 0.45 log2 on CD14+ monocytes
    while the README said they agree to three decimals.
    """
    sf = np.asarray(sf, dtype=float)
    l2 = np.log2(counts / sf[:, None] + 1.0)
    return l2[groups_a].mean(axis=0), l2[groups_b].mean(axis=0)


def _nb_irls(y, X, offset, alpha, max_iter=100, tol=1e-8):
    """Fit an NB GLM with fixed dispersion by IRLS. Log link.

    With alpha fixed this is a plain GLM: score X'(y-mu)/V and working
    weights mu^2/V, V = mu + alpha*mu^2. Returns (beta, cov).
    """
    n, p = X.shape
    y = np.asarray(y, dtype=float)
    # start from a log-scale least squares fit; +0.5 keeps zeros defined
    with np.errstate(divide="ignore", invalid="ignore"):
        y0 = np.log(np.maximum(y, 0.0) + 0.5) - offset
    beta, *_ = np.linalg.lstsq(X, y0, rcond=None)
    if not np.all(np.isfinite(beta)):
        beta = np.zeros(p)

    for _ in range(max_iter):
        eta = np.clip(X @ beta + offset, -100.0, 100.0)
        mu = np.exp(eta)
        w = mu / (1.0 + alpha * mu)            # = mu^2 / V
        z = eta - offset + (y - mu) / mu
        xtwx = X.T @ (X * w[:, None])
        xtwz = X.T @ (w * z)
        try:
            new = np.linalg.solve(xtwx, xtwz)
        except np.linalg.LinAlgError:
            new, *_ = np.linalg.lstsq(xtwx, xtwz, rcond=None)
        if not np.all(np.isfinite(new)):
            break
        # damped update: IRLS on a log link can overshoot at large mu
        delta = np.clip(new - beta, -5.0, 5.0)
        beta = np.clip(beta + delta, -200.0, 200.0)
        if np.max(np.abs(delta)) < tol:
            break

    eta = np.clip(X @ beta + offset, -100.0, 100.0)
    mu = np.exp(eta)
    w = mu / (1.0 + alpha * mu)
    xtwx = X.T @ (X * w[:, None])
    try:
        cov = np.linalg.inv(xtwx)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(xtwx)
    return beta, cov, mu


def _dispersion_trend(base_mean, disp):
    """Smooth log(disp) against log(mean): a straight line, robustly fitted.

    Replaces a lowess curve, which at these sample sizes is mostly noise.
    Trims the extreme 10% at each end before fitting so a few very noisy
    genes cannot tilt the line.
    """
    ok = (base_mean > 0) & np.isfinite(disp) & np.isfinite(base_mean)
    lx = np.log(base_mean[ok])
    ly = np.log(np.clip(disp[ok], DISP_FLOOR, 10.0))
    if len(lx) < 10:
        const = np.median(ly) if len(ly) else np.log(0.1)
        return np.full(len(base_mean), np.exp(const))
    lo, hi = np.quantile(lx, [0.05, 0.95])
    keep = (lx >= lo) & (lx <= hi)
    if keep.sum() < 10:
        keep = np.ones_like(lx, dtype=bool)
    try:
        slope, intercept = np.polyfit(lx[keep], ly[keep], 1)
        slope = float(np.clip(slope, -1.0, 0.0))    # dispersion only falls
    except Exception:
        slope, intercept = 0.0, float(np.median(ly[keep]))
    pred = intercept + slope * np.log(np.maximum(base_mean, 1e-8))
    return np.exp(np.clip(pred, np.log(DISP_FLOOR), np.log(10.0)))


def _shrink(disp_mom, disp_trend, n_samples, n_params):
    """Weighted geometric-ish blend, trend-heavy when information is thin.

    n_eff = samples - parameters. A gene fitted from 6 pseudobulk samples
    with 2 parameters gets a weight of ~0.3 on its own estimate; at 48
    samples it gets ~0.8. This is the empirical-Bayes step in miniature:
    its job is to stop a gene whose MoM dispersion came out too low from
    producing a spuriously small standard error.
    """
    n_eff = max(n_samples - n_params, 1)
    w = n_eff / (n_eff + 10.0)
    blended = np.exp(w * np.log(np.clip(disp_mom, DISP_FLOOR, 10.0))
                     + (1 - w) * np.log(np.clip(disp_trend, DISP_FLOOR, 10.0)))
    # Never let shrinkage push a gene BELOW the trend. Shrinking downward is
    # the direction that manufactures false positives. CAVEAT on provenance:
    # this floor was added while tuning on the repo's own simulator at 6
    # donors (an earlier iteration measured an FDR of 0.066 with the floor and
    # 0.13 without it; those are not re-measured figures), so the simulated
    # power and FDP in scripts/nbglm_study.py are not out-of-sample for this
    # choice. The simulator also draws negative-binomial counts, the GLM's own
    # model. Real-data behaviour is reported separately (scripts/nbglm_real.py).
    return np.clip(np.maximum(blended, disp_trend), DISP_FLOOR, 10.0)


def nbglm_de(adata, condition_key, sample_key, group_a, group_b, layer=None,
             min_cells=10, paired=True, min_count=10, min_samples=3,
             alpha=None):
    """Negative-binomial GLM test of group_a vs group_b across donors.

    Returns the same columns as `de.pseudobulk_de` so the two can be put in
    the same table: gene, mean_a, mean_b, log2fc, pval, qval, n_a, n_b,
    plus dispersion, log2fc_glm and fit_status.

    log2fc is on the same normalisation as de.pseudobulk_de (median-of-ratios
    size factors), so the two are directly comparable. log2fc_glm is the model
    coefficient and is NaN for failed fits; for separated genes
    (fit_status == "separated", |log2fc_glm| > SEPARATION_LOG2) it is
    unbounded and must not be read as an effect size.
    """
    mat, meta = aggregate(adata, sample_key=sample_key,
                          group_key=condition_key, layer=layer,
                          min_cells=min_cells)
    mat = filter_genes(mat, min_count=min_count, min_samples=min_samples)

    keep = meta[condition_key].isin([group_a, group_b])
    mat = mat.loc[keep]
    meta = meta.loc[keep]
    if len(mat) < 4:
        raise ValueError(f"need at least 4 pseudobulk samples, got {len(mat)}")

    sf = np.array(size_factors(mat).to_numpy(), dtype=float)
    bad_sf = int((~(sf > 0)).sum())
    if bad_sf:
        # a zero / negative / NaN size factor means a sample with no usable
        # counts; replacing it by 1 keeps the fit running but the sample is
        # then not depth-corrected, so say so rather than mask it
        warnings.warn(f"{bad_sf} sample(s) with a non-positive size factor "
                      "were set to 1.0; their depth is NOT corrected",
                      RuntimeWarning, stacklevel=2)
    sf[~(sf > 0)] = 1.0
    offset = np.log(sf)
    counts = mat.to_numpy(dtype=float)

    cond = (meta[condition_key] == group_b).to_numpy(dtype=float)
    donors = meta[sample_key].astype(str).to_numpy()
    if paired:
        # A donor with only one condition contributes no information about
        # the contrast and would blow up its own dummy level, so drop it
        # before the design matrix is built.
        keep_rows = np.ones(len(counts), dtype=bool)
        for d in np.unique(donors):
            m = donors == d
            if len(np.unique(cond[m])) < 2:
                keep_rows &= ~m
        if keep_rows.sum() < 4:
            raise ValueError("too few complete donor pairs for a paired GLM")
        counts = counts[keep_rows]
        offset = offset[keep_rows]
        cond = cond[keep_rows]
        donors = donors[keep_rows]

    cols = [np.ones((len(counts), 1)), cond[:, None]]
    if paired:
        for lev in sorted(set(donors))[1:]:
            cols.append((donors == lev).astype(float)[:, None])
    X = np.hstack(cols)

    n_samples, n_params = X.shape
    genes = mat.columns.to_numpy()

    # --- dispersion -------------------------------------------------------
    if alpha is not None:
        disp = np.full(len(genes), float(alpha))
    else:
        beta_p, mus = [], []
        for g in range(len(genes)):     # Poisson pass, for the MoM baseline
            b, _, m = _nb_irls(counts[:, g], X, offset, 1e-8, max_iter=40)
            beta_p.append(b)
            mus.append(m)
        mu_all = np.vstack(mus).T
        df = max(n_samples - n_params, 1)
        # E[RSS] = sum_i V_i (1 - h_ii), and h_ii averages p/n, so
        # RSS * n / (n - p) is an unbiased estimate of sum_i V_i. Skipping
        # the n/(n-p) factor is what makes a naive MoM dispersion too small
        # and the Wald test anti-conservative at 6 donors.
        scale = n_samples / df
        numer = np.maximum(((counts - mu_all) ** 2).sum(axis=0) * scale
                           - mu_all.sum(axis=0), 0.0)
        denom = (mu_all ** 2).sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            disp_mom = numer / np.maximum(denom, 1e-12)
        disp_mom = np.clip(np.nan_to_num(disp_mom, nan=0.1), DISP_FLOOR, 10.0)

        base_mean = mu_all.mean(axis=0)
        trend = _dispersion_trend(base_mean, disp_mom)
        disp = _shrink(disp_mom, trend, n_samples, n_params)

    # --- Wald test on the condition coefficient ---------------------------
    pvals = np.full(len(genes), 1.0)
    coefs = np.full(len(genes), np.nan)
    ses = np.full(len(genes), np.nan)
    status = np.full(len(genes), "ok", dtype=object)
    for g in range(len(genes)):
        y = counts[:, g]
        if y.max() <= 0:
            pvals[g] = 1.0                  # nothing to test
            status[g] = "no_counts"
            continue
        b, cov, _ = _nb_irls(y, X, offset, disp[g])
        var = cov[1, 1]
        if not np.isfinite(var) or var <= 0:
            # p is left at 1 (so the BH family keeps its size and the gene is
            # never called) but the coefficient is NOT reported as 0: that
            # would invent an estimate. These are failed fits, not nulls.
            pvals[g] = 1.0
            status[g] = "failed"
            continue
        coefs[g] = b[1]
        if abs(b[1]) / np.log(2.0) > SEPARATION_LOG2:
            status[g] = "separated"
        ses[g] = np.sqrt(var)
        z = b[1] / ses[g]
        # Student-t, not normal. The dispersion is estimated from few
        # samples, so the Wald statistic is only t-distributed; using a
        # normal reference here is what makes NB GLMs anti-conservative at
        # 3-6 donors (measured: 18% of null genes at p<0.05 with 6 donors).
        pvals[g] = 2.0 * stats.t.sf(abs(z), df=max(n_samples - n_params, 1))

    is_a = cond == 0.0
    is_b = cond == 1.0
    mean_a, mean_b = _log2norm_means(counts, is_a, is_b, np.exp(offset))

    return pd.DataFrame({
        "gene": genes,
        "mean_a": mean_a,
        "mean_b": mean_b,
        "log2fc": mean_a - mean_b,    # log2 of size-factor-normalised counts
        "log2fc_glm": -coefs / np.log(2.0),        # model coefficient
        "fit_status": status,    # ok / separated / failed / no_counts
        "pval": pvals,
        "qval": _bh(pvals),
        "n_a": int(is_a.sum()),
        "n_b": int(is_b.sum()),
        "dispersion": disp,
    }, index=range(len(genes))).set_index("gene")


def _norm_sf(z):
    """Standard normal survival function, vectorised over scalars."""
    from scipy.special import erfc
    return 0.5 * erfc(np.asarray(z) / np.sqrt(2.0))
