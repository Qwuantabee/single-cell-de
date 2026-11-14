"""Tests for the mixed-effects DE path.

A mixed model that "works" is easy to write and easy to get wrong: it can
return the pooled estimate while looking like it accounts for donors, or it
can silently return numbers for genes where the optimiser never converged.
These tests pin both, because the interesting failure is a plausible-looking
wrong answer, not an exception.
"""

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from scpipe import mixed


def _profiles(n_donors=8, n_genes=12, effect=1.5, donor_sd=1.0, seed=0,
              noise=0.15):
    """Pseudobulk-like frame indexed donor||condition."""
    rng = np.random.default_rng(seed)
    donors = [f"d{i}" for i in range(n_donors)]
    row_index = [f"{d}||{c}" for d in donors for c in ("ctrl", "stim")]
    genes = [f"g{i:02d}" for i in range(n_genes)]

    donor_shift = rng.normal(0, donor_sd, size=(n_donors, n_genes))
    vals = {}
    for i, d in enumerate(donors):
        for c in ("ctrl", "stim"):
            base = donor_shift[i] + 5.0
            shift = effect if c == "stim" else 0.0
            vals[f"{d}||{c}"] = base + shift + rng.normal(0, noise, n_genes)

    frame = pd.DataFrame(vals, index=genes).T
    return frame


def test_recovers_the_planted_effect():
    """A constant 1.5 shift in every gene must come back as ~1.5.

    The donor offsets are large (sd 1.0) and the within-pair noise small, so
    a model that ignored donors would still find the effect but with a much
    larger standard error. The coefficient itself is what is checked here.
    """
    frame = _profiles(effect=1.5, donor_sd=1.0, noise=0.15)
    out = mixed.mixed_de(frame, "donor", "condition", "stim", "ctrl")
    assert out["converged"].mean() > 0.9
    assert out["coef"].mean() == pytest.approx(1.5, abs=0.2)
    assert (out["qval"].dropna() < 0.05).mean() > 0.9


def test_no_effect_gives_no_calls():
    """Null calibration, across seeds rather than one draw.

    One seed cannot decide this: a single FDR-controlling test is allowed to
    call something. The claim being tested is that the *rate* of calls under
    the null stays at or below alpha, so the rate is what is asserted, over
    enough genes for the rate to mean something.
    """
    calls = total = 0
    for seed in range(6):
        frame = _profiles(n_donors=10, n_genes=20, effect=0.0,
                          donor_sd=0.3, noise=0.5, seed=seed)
        out = mixed.mixed_de(frame, "donor", "condition", "stim", "ctrl")
        calls += int((out["qval"] < 0.05).sum())
        total += len(out)
    assert calls / total <= 0.05, f"{calls}/{total} false calls under null"


def test_smaller_effect_is_smaller_coefficient():
    """Monotonic in the truth, so the number means something."""
    small = mixed.mixed_de(_profiles(effect=0.5, seed=1), "donor",
                           "condition", "stim", "ctrl")
    large = mixed.mixed_de(_profiles(effect=2.0, seed=1), "donor",
                           "condition", "stim", "ctrl")
    assert large["coef"].mean() > small["coef"].mean() + 1.0


def test_sign_and_orientation():
    """stim minus ctrl, not the other way round, and not |effect|."""
    frame = _profiles(effect=1.5)
    out = mixed.mixed_de(frame, "donor", "condition", "stim", "ctrl")
    assert (out["coef"] > 0).all()
    flipped = mixed.mixed_de(frame, "donor", "condition", "ctrl", "stim")
    assert (flipped["coef"] < 0).all()


def test_index_and_columns():
    frame = _profiles(n_genes=6)
    out = mixed.mixed_de(frame, "donor", "condition", "stim", "ctrl")
    assert list(out.index) == list(frame.columns)
    for col in ("log2fc", "pval", "qval", "converged", "donor_var", "n_a",
                "n_b"):
        assert col in out.columns
    assert (out["n_a"] == 8).all() and (out["n_b"] == 8).all()
    # log2fc is the fitted coefficient, not a recomputed mean difference
    assert np.allclose(out["log2fc"], out["coef"])


def test_qvalues_bounded_by_one():
    out = mixed.mixed_de(_profiles(n_genes=40, seed=2), "donor", "condition",
                         "stim", "ctrl")
    q = out["qval"].dropna()
    assert (q >= 0).all() and (q <= 1).all()


def test_donor_variance_is_reported():
    """The random effect has to actually be there; a model with a zero
    variance component is a pooled model wearing the name."""
    frame = _profiles(effect=1.5, donor_sd=1.0, noise=0.1)
    out = mixed.mixed_de(frame, "donor", "condition", "stim", "ctrl")
    assert (out["donor_var"].dropna() > 0).mean() > 0.5


def test_rejects_bad_index():
    frame = pd.DataFrame({"g0": [1.0, 2.0, 3.0, 4.0]},
                         index=["a", "b", "c", "d"])
    with pytest.raises(ValueError, match="must look like"):
        mixed.mixed_de(frame, "donor", "condition", "stim", "ctrl")


def test_rejects_single_sample_per_condition():
    frame = _profiles(n_donors=1, n_genes=4)
    with pytest.raises(ValueError, match="need >=2 samples per condition"):
        mixed.mixed_de(frame, "donor", "condition", "stim", "ctrl")


def test_mixed_pvalue_equals_paired_t_on_balanced_fixture():
    """The mixed p-value must be the PAIRED t p-value, not a normal one.

    On a balanced two-sample-per-donor design the fitted condition
    coefficient IS the paired mean difference, so the only defensible
    reference distribution is Student t on n_pairs-1. statsmodels returns a
    Wald p on a NORMAL reference, which is the bug this test pins: with 8
    donors, normal vs t(7) is the difference between 1872 and 1033 CD4 calls.

    donor_sd is positive so the random effect is genuinely estimated and the
    equality is not a degenerate boundary case.
    """
    frame = _profiles(n_donors=8, n_genes=15, effect=1.2, donor_sd=1.0,
                      noise=0.4, seed=3)
    out = mixed.mixed_de(frame, "donor", "condition", "stim", "ctrl")
    assert out["df"].iloc[0] == 7
    assert (out["donor_var"].dropna() > 0).mean() > 0.5

    donors = [f"d{i}" for i in range(8)]
    stim = frame.loc[[f"{d}||stim" for d in donors]].to_numpy(dtype=float)
    ctrl = frame.loc[[f"{d}||ctrl" for d in donors]].to_numpy(dtype=float)
    _, p_paired = stats.ttest_rel(stim, ctrl, axis=0)

    assert out["pval"].to_numpy() == pytest.approx(p_paired, rel=1e-6,
                                                   abs=1e-9)
    # and the coefficient is the mean difference, which is why it holds
    assert out["coef"].to_numpy() == pytest.approx(
        (stim - ctrl).mean(axis=0), abs=1e-9)


def test_percell_agrees_with_profiles_when_cells_are_identical():
    """The per-cell fit is a different data path to the same estimand.

    Cells are built so that each donor's sample values are exactly the
    pseudobulk profile repeated, which makes the two models estimands
    identical and only the variance structure different.
    """
    profiles = _profiles(n_donors=6, n_genes=5, effect=1.2, noise=0.05)
    rows, donor, cond = [], [], []
    for key in profiles.index:
        d, c = key.split("||")
        for _ in range(6):
            rows.append(profiles.loc[key].to_numpy(dtype=float))
            donor.append(d)
            cond.append(c)
    out = mixed.mixed_de_percell(np.vstack(rows), profiles.columns, donor,
                                 cond, "stim", "ctrl")
    expected = mixed.mixed_de(profiles, "donor", "condition", "stim", "ctrl")
    assert out["log2fc"].to_numpy() == pytest.approx(
        expected.loc[out.index, "log2fc"].to_numpy(), abs=0.05)
    assert out["converged"].all()
