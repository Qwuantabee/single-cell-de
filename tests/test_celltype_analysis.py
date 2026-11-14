"""Small pure helpers of scripts/celltype_analysis.py."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "celltype_analysis.py"
spec = importlib.util.spec_from_file_location("celltype_analysis", SCRIPT)
ca = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ca)


def test_marker_share_counts_the_union_of_type_markers():
    novel = {"g1", "g2", "g3", "g4"}
    share, per_type = ca.marker_share(
        novel, {"A": ["g1", "g2", "zz"], "B": ["g2", "g3"]})
    assert per_type == {"A": 2, "B": 2}
    assert share == pytest.approx(3 / 4)       # g4 is nobody's marker


def test_marker_share_is_nan_not_zero_when_nothing_is_pooled_only():
    share, per_type = ca.marker_share(set(), {"A": ["g1"]})
    assert np.isnan(share) and per_type == {"A": 0}


def test_marker_share_is_zero_when_markers_do_not_overlap():
    share, _ = ca.marker_share({"g1"}, {"A": ["x"], "B": ["y"]})
    assert share == 0.0
