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

import numpy as np
import pandas as pd

from .de import _bh
from .pseudobulk import aggregate, filter_genes, size_factors

DISP_FLOOR = 1e-4


def _log2cpm_means(counts, groups_a, groups_b):
    """Mean log2 CPM per group, for reporting fold change only."""
    lib = counts.sum(axis=1)
    lib = np.where(lib > 0, lib, 1.0)
    l2 = np.log2(counts / lib[:, None] * 1e6 + 1.0)
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
    # never let shrinkage take a gene below a fraction of the trend: that is
    # the direction that manufactures false positives
    return np.clip(np.maximum(blended, 0.25 * disp_trend), DISP_FLOOR, 10.0)


def nbglm_de(adata, condition_key, sample_key, group_a, group_b, layer=None,
             min_cells=10, paired=True, min_count=10, min_samples=3,
             alpha=None):
    """Negative-binomial GLM test of group_a vs group_b across donors.

    Returns the same columns as `de.pseudobulk_de` so the two can be put in
    the same table: gene, mean_a, mean_b, log2fc, pval, qval, n_a, n_b,
    plus dispersion and log2fc_glm.
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
    sf[sf <= 0] = 1.0
    offset = np.log(sf)
    counts = mat.to_numpy(dtype=float)

    cond = (meta[condition_key] == group_b).to_numpy(dtype=float)[:, None]
    cols = [np.ones((len(mat), 1)), cond]
    if paired:
        donors = meta[sample_key].astype(str).to_numpy()
        # keep donors contributing to both conditions; drop the rest
        for d in np.unique(donors):
            m = donors == d
            if len(np.unique(cond[m, 0])) == 2:
                continue
            counts = counts[~m]
            offset = offset[~m]
            cond = cond[~m]
            donors = donors[~m]
        if len(counts) < 4:
            raise ValueError("too few complete donor pairs for a paired GLM")
        levels = sorted(set(donors))
        for lev in levels[1:]:
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
    coefs = np.zeros(len(genes))
    ses = np.full(len(genes), np.nan)
    for g in range(len(genes)):
        y = counts[:, g]
        if y.max() <= 0:
            pvals[g] = 1.0                  # nothing to test
            continue
        b, cov, _ = _nb_irls(y, X, offset, disp[g])
        var = cov[1, 1]
        if not np.isfinite(var) or var <= 0:
            pvals[g] = 1.0
            continue
        coefs[g] = b[1]
        ses[g] = np.sqrt(var)
        z = b[1] / ses[g]
        pvals[g] = 2.0 * _norm_sf(abs(z))

    is_a = (meta[condition_key] == group_a).to_numpy()
    is_b = (meta[condition_key] == group_b).to_numpy()
    mean_a, mean_b = _log2cpm_means(counts, is_a, is_b)

    return pd.DataFrame({
        "gene": genes,
        "mean_a": mean_a,
        "mean_b": mean_b,
        "log2fc": mean_a - mean_b,                 # log2 CPM, as in de.py
        "log2fc_glm": -coefs / np.log(2.0),        # model coefficient
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
