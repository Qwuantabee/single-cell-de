"""Differential expression tests.

The important ones are the calibration tests, which assert that each method
behaves correctly on data where the right answer is known. They are the
reason this package exists.
"""

import numpy as np
import pytest

from scpipe import de, simulate


def test_bh_matches_known_values():
    """Benjamini-Hochberg on a hand-checkable example."""
    p = np.array([0.01, 0.02, 0.03, 0.04, 0.05])
    q = de._bh(p)
    # q_i = min over j>=i of (p_j * n / j)
    assert q[0] == pytest.approx(0.05)
    assert np.all(np.diff(q) >= -1e-12)       # monotone non-decreasing
    assert np.all(q <= 1.0)


def test_bh_handles_nan():
    p = np.array([0.01, np.nan, 0.5])
    q = de._bh(p)
    assert np.isnan(q[1])
    assert np.isfinite(q[0]) and np.isfinite(q[2])


def test_pseudobulk_controls_false_positives_under_the_null():
    """Nothing changes, so almost nothing should be called.

    This is the headline claim of the project and the test that would catch
    a regression in the aggregation or the normalisation.
    """
    called = []
    for seed in range(5):
        adata, _ = simulate.simulate(n_genes=300, n_de=0, seed=seed)
        res = de.pseudobulk_de(adata, "condition", "donor", "stim", "ctrl",
                               paired=True, min_count=5, min_samples=3)
        called.append(de.n_significant(res, 0.05))
    # expect roughly 5% of 300 = 15 at worst; in practice far fewer
    assert np.mean(called) < 15, f"called {np.mean(called):.1f} on null data"


def test_per_cell_test_is_anticonservative_under_the_null():
    """The negative result, asserted rather than assumed.

    If this ever stops failing, either the simulator stopped generating
    sample-level variation or the per-cell test was quietly fixed, and
    either way the project's conclusion needs rechecking.
    """
    adata, _ = simulate.simulate(n_genes=300, n_de=0, seed=0)
    res = de.single_cell_de(adata, "condition", "stim", "ctrl")
    # it should be badly wrong, not marginally wrong
    assert de.n_significant(res, 0.05) > 50


def test_per_cell_is_calibrated_when_there_is_no_sample_variation():
    """The diagnosis.

    With sample_sd=0 every cell really is an independent draw from the same
    distribution, the per-cell test's assumption holds, and it behaves. This
    is what identifies unmodelled sample-level variance as the cause rather
    than something about cell counts or the rank test itself.
    """
    adata, _ = simulate.simulate(n_genes=300, n_de=0, sample_sd=0.0,
                                 donor_sd=0.5, seed=0)
    res = de.single_cell_de(adata, "condition", "stim", "ctrl")
    assert de.n_significant(res, 0.05) < 15


def test_pseudobulk_recovers_true_effects():
    adata, truth = simulate.simulate(n_genes=300, n_de=30, effect=1.0, seed=7)
    res = de.pseudobulk_de(adata, "condition", "donor", "stim", "ctrl",
                           paired=True, min_count=5, min_samples=3)
    m = simulate.confusion(res, truth, 0.05)
    assert m["power"] > 0.8
    assert m["fdp"] < 0.20


def test_log2fc_sign_matches_the_simulated_direction():
    adata, truth = simulate.simulate(n_genes=200, n_de=20, effect=1.5, seed=3)
    res = de.pseudobulk_de(adata, "condition", "donor", "stim", "ctrl",
                           paired=True, min_count=5, min_samples=3)
    sig = res[(res["qval"] < 0.05) & (res.index.isin(truth))]
    assert len(sig) > 5
    # effects were assigned +/- at random, so both directions must appear
    assert (sig["log2fc"] > 0).any() and (sig["log2fc"] < 0).any()


def test_pseudobulk_needs_enough_donors():
    """One donor per condition gives no way to estimate variance.

    Refusing is the right behaviour: with a single replicate the test has no
    denominator, and anything it returned would be fabricated confidence.
    """
    adata, _ = simulate.simulate(n_genes=100, n_donors=1, seed=0)
    with pytest.raises(ValueError, match="at least 2 donors"):
        de.pseudobulk_de(adata, "condition", "donor", "stim", "ctrl")


def test_unknown_normalisation_rejected():
    adata, _ = simulate.simulate(n_genes=100, seed=0)
    with pytest.raises(ValueError, match="unknown normalisation"):
        de.pseudobulk_de(adata, "condition", "donor", "stim", "ctrl",
                         normalisation="quantile")


def test_missing_condition_level_rejected():
    adata, _ = simulate.simulate(n_genes=100, seed=0)
    with pytest.raises(ValueError, match="no cells"):
        de.single_cell_de(adata, "condition", "stim", "nonexistent")
