"""The BH null reference is simulated by a committed script, not asserted."""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import bh_null_expectation as bne  # noqa: E402
from scpipe import de  # noqa: E402


def test_vectorised_bh_count_matches_the_library_bh():
    """The simulation's rejection count must be the same BH the analysis
    uses, otherwise the reference is for a different procedure."""
    rng = np.random.default_rng(0)
    p = rng.random((60, 200)) ** 3          # plenty of small values
    k = bne.bh_rejections(p, 0.05)
    for row, kk in zip(p, k):
        assert int((de._bh(row) < 0.05).sum()) == int(kk)


def test_bh_count_on_a_known_family():
    p = np.array([[0.001, 0.008, 0.039, 0.041, 0.6, 0.9]])
    # thresholds i*0.05/6: 0.0083, 0.0167, 0.025, 0.033, ...; p(2)=0.008
    # passes, p(3)=0.039 > 0.025 does not -> 2 rejections
    assert int(bne.bh_rejections(p, 0.05)[0]) == 2
    assert int(bne.bh_rejections(np.array([[0.5, 0.6, 0.7]]), 0.05)[0]) == 0


def test_simulated_reference_is_about_a_twentieth_not_25():
    """P(any rejection) is q = 0.05 for independent nulls, and the mean count
    is only slightly above that. 25 (alpha * m) is the uncorrected figure."""
    out = bne.simulate(m=500, draws=4000, q=0.05, seed=1)
    assert out["p_any"] == pytest.approx(0.05, abs=4 * out["p_any_se"] + 0.005)
    assert 0.03 < out["mean_rejections"] < 0.09
    assert out["mean_rejections"] < 1.0            # nowhere near alpha * m


def test_script_constant_agrees_with_the_simulation():
    import mixed_vs_others as mvo
    out = bne.simulate(m=500, draws=4000, q=0.05, seed=2)
    assert abs(mvo.BH_EXPECTED_FP - out["mean_rejections"]) < max(
        0.02, 4 * out["mean_rejections_se"])
