"""Mixed-effects differential expression: one random intercept per donor.

Pseudobulk works by throwing cells away. It collapses a donor's sample into
a single profile so that the replicate count equals the number of donors,
which is the right number of degrees of freedom but also discards all of the
within-sample information. A mixed-effects model keeps the cells and still
lets the donor carry its own variance component:

    expression ~ condition + (1 | donor)

The fixed effect of condition is then estimated from the within-donor
contrast, which is what a paired design actually randomises. Donor-to-donor
spread is not treated as replication.

Two things here were measured rather than assumed, and both are why this
module is not a three-line wrapper around statsmodels.

1. The optimiser matters. On this design -- two samples per donor, so the
   random intercept is estimated from one pair per donor -- the default L-BFGS
   solver in statsmodels raises `LinAlgError: Singular matrix` on about two
   thirds of genes, because the variance component hits the boundary and the
   information matrix becomes numerically singular. Powell converges on the
   same genes without the failure. The default here is therefore Powell, and
   the alternative is available for anyone who wants to check the claim.

2. The variance component is often genuinely near zero. With one pair per
   donor the data cannot tell a small donor variance apart from none, and the
   reported value should not be read as an estimate of donor spread. It is
   reported anyway, because a caller that wants to know whether the random
   effect is doing anything needs to see it.

Cost: fitting one model per gene is orders of magnitude slower than the
vectorised t-test, so per-cell fits are expected to be run on a deliberate
subsample. scripts/mixed_vs_others.py measures both walls.
"""

import warnings

import numpy as np
import pandas as pd
from scipy import stats

from .de import _bh

# a variance component below this is a boundary fit, not an estimate
BOUNDARY_VAR = 1e-6


def _fit(y, cond, donor, method="powell", maxiter=1000, df=np.inf):
    """One gene. Returns (coef, pvalue, converged, donor_variance, boundary).

    df is the reference degrees of freedom for the fixed-effect test. The
    Wald p-value statsmodels returns uses a NORMAL reference (df=inf). On this
    design the coefficient IS the paired mean difference, so the right
    reference is Student t on n_pairs-1: re-referencing the same CD4 estimates
    from normal to t(7) cut BH-significant calls from 1872 to 1033 (paired t
    gives 986). The "anti-conservative inflation" was a distributional
    default, not a property of mixed models.
    """
    import statsmodels.formula.api as smf

    frame = pd.DataFrame({"y": y, "cond": cond, "donor": donor})
    with warnings.catch_warnings():
        # a boundary variance component is expected on this design, not news
        warnings.simplefilter("ignore")
        try:
            fit = smf.mixedlm("y ~ cond", frame, groups=frame["donor"]).fit(
                method=method, maxiter=maxiter)
        except (np.linalg.LinAlgError, ValueError, np.linalg.Error):
            return np.nan, np.nan, False, np.nan, False

    coef = fit.params.get("cond", np.nan)
    try:
        donor_var = float(fit.cov_re.iloc[0, 0])
    except (AttributeError, IndexError):
        donor_var = np.nan

    # the standard error from the reported Wald bse, then re-reference on t
    try:
        se = float(fit.bse.get("cond", np.nan))
    except (AttributeError, TypeError):
        se = np.nan
    ok = (bool(fit.converged) and np.isfinite(coef) and np.isfinite(se)
          and se > 0)
    if not ok:
        # a p-value from a fit that did not converge is not a p-value
        coef_out = float(coef) if np.isfinite(coef) else np.nan
        return coef_out, np.nan, False, donor_var, False

    t_stat = coef / se
    if np.isfinite(df):
        pval = float(2 * stats.t.sf(abs(t_stat), df))
    else:
        pval = float(2 * stats.norm.sf(abs(t_stat)))

    # boundary fit: the random intercept collapsed, so the model is really
    # pooled. Reported separately from convergence, not folded into it.
    boundary = bool(not np.isfinite(donor_var) or donor_var < BOUNDARY_VAR)
    return float(coef), pval, True, donor_var, boundary


def _split_index(index, sample_key, condition_key):
    """Recover sample/condition labels from 'sample||condition' keys."""
    # "||" is a regex alternation, not a literal, so the default regex=True
    # would tear "d0||ctrl" into one character per column.
    parts = pd.Series(index).str.split("||", expand=True, regex=False)
    if parts.shape[1] < 2:
        raise ValueError(
            "index must look like 'sample||condition'; got "
            f"{list(index[:2])}")
    return (parts[0].to_numpy().astype(str),
            parts[1].to_numpy().astype(str))


def mixed_de(logcpm, sample_key, condition_key, group_a, group_b,
             method="powell", maxiter=1000, donors=None):
    """Fit `expression ~ condition + (1 | sample)` once per gene.

    logcpm: samples x genes, already normalised (log2 CPM from
        pseudobulk.normalise). Passing raw counts fits a Gaussian model to
        counts, which is a different and worse model.
    sample_key: the donor. Use the 10x lane here and this becomes the
        per-cell test wearing a hat.
    donors: optional whitelist of donors. Use this to fit the mixed model on
        exactly the profiles a paired test uses: without it the model also
        absorbs single-condition profiles from donors that cannot be paired,
        which changes both the estimate and the size factors, and then the
        two arms are no longer comparable.

    Returns one row per gene with the condition coefficient, its p-value,
    whether the fit converged, and the donor variance component.
    Non-converged fits get pval=NaN rather than a discarded number: dropping
    them silently would bias the result towards the easy genes.
    """
    donor_labs, conditions = _split_index(logcpm.index, sample_key,
                                          condition_key)
    is_a = conditions == group_a
    is_b = conditions == group_b
    if is_a.sum() < 2 or is_b.sum() < 2:
        raise ValueError(
            f"need >=2 samples per condition, got {is_a.sum()}/{is_b.sum()}")

    keep = is_a | is_b
    if donors is not None:
        keep &= np.isin(donor_labs, np.asarray(donors).astype(str))
    y_all = logcpm.to_numpy(dtype=float)[keep]
    cond = is_a[keep].astype(float)
    donor = donor_labs[keep]

    # degrees of freedom: the number of DONORS contributing both conditions,
    # i.e. the number of pairs actually fitted. The fixed effect is estimated
    # from the within-donor contrast, so that count -- not the row count -- is
    # the replication, and the reference distribution has to be t on it.
    n_pairs = 0
    for d in np.unique(donor):
        own = cond[donor == d]
        if (own == 1.0).any() and (own == 0.0).any():
            n_pairs += 1
    df = n_pairs - 1 if n_pairs >= 2 else np.inf

    rows = []
    for j, gene in enumerate(logcpm.columns):
        coef, pval, ok, dvar, boundary = _fit(
            y_all[:, j], cond, donor, method, maxiter, df=df)
        rows.append((gene, coef, pval, ok, dvar, boundary))

    out = pd.DataFrame(rows, columns=[
        "gene", "coef", "pval", "converged", "donor_var", "boundary"]
    ).set_index("gene")
    out["qval"] = _bh(out["pval"].to_numpy())
    out["log2fc"] = out["coef"]
    out["df"] = df
    out["n_pairs"] = n_pairs
    out["n_a"] = int((is_a & keep).sum())
    out["n_b"] = int((is_b & keep).sum())
    return out


def mixed_de_percell(values, genes, donor, condition, group_a, group_b,
                     method="powell", maxiter=1000):
    """Mixed model on individual cells rather than pseudobulk profiles.

    values: cells x genes matrix, dense, already normalised. The caller
        densifies, which is one reason per-cell fits are run on a subsample.
    donor / condition: per-cell labels.

    The contrast is estimated within donor, so cells are not treated as
    independent replicates -- that is the difference from single_cell_de.
    """
    values = np.asarray(values, dtype=float)
    donor = np.asarray(donor).astype(str)
    condition = np.asarray(condition).astype(str)

    is_a = condition == group_a
    is_b = condition == group_b
    keep = is_a | is_b
    values, donor, cond = values[keep], donor[keep], is_a[keep].astype(float)

    rows = []
    pairs = [d for d in np.unique(donor)
             if ((cond == 1.0) & (donor == d)).any()
             and ((cond == 0.0) & (donor == d)).any()]
    df = len(pairs) - 1 if len(pairs) >= 2 else np.inf
    for j, gene in enumerate(genes):
        coef, pval, ok, _dvar, _bd = _fit(values[:, j], cond, donor, method,
                                          maxiter, df=df)
        rows.append((gene, coef, pval, ok))

    out = pd.DataFrame(rows, columns=["gene", "log2fc", "pval", "converged"]
                       ).set_index("gene")
    out["qval"] = _bh(out["pval"].to_numpy())
    out["n_a"] = int(is_a.sum())
    out["n_b"] = int(is_b.sum())
    return out
