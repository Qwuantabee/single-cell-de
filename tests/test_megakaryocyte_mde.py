"""The megakaryocyte power calculation, checked against simulation."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

SCRIPT = (Path(__file__).resolve().parent.parent / "scripts"
          / "megakaryocyte_check.py")
spec = importlib.util.spec_from_file_location("megakaryocyte_check", SCRIPT)
mk = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mk)


def test_exact_mde_is_at_least_the_central_t_approximation():
    for n in (3, 4, 8):
        approx = mk.min_detectable_effect(0.8, n)
        exact = mk.exact_min_detectable_effect(0.8, n)
        assert exact >= approx - 1e-9
        assert exact < 1.5 * approx          # same ballpark, not a different test


def test_exact_mde_has_the_stated_power_in_simulation():
    """Draw paired differences with true mean = MDE and check ~80% reject."""
    from scipy import stats
    sd, n = 0.8, 4
    mde = mk.exact_min_detectable_effect(sd, n)
    rng = np.random.default_rng(0)
    x = rng.normal(mde, sd, size=(40000, n))
    p = stats.ttest_1samp(x, 0.0, axis=1).pvalue
    assert np.mean(p < 0.05) == pytest.approx(0.80, abs=0.01)


def test_bonferroni_mde_is_larger():
    base = mk.exact_min_detectable_effect(0.8, 4)
    bonf = mk.exact_min_detectable_effect(0.8, 4, alpha=0.05 / 85)
    assert bonf > 2 * base
    assert mk.exact_min_detectable_effect(float("nan"), 4) == float("inf")
