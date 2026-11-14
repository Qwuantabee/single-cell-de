"""The donor-swap null must permute labels, never cells.

A permutation null is only as good as its implementation, and the failure
mode here is silent: shuffling condition labels across *cells* also produces
a null-looking distribution, but it destroys the donor structure the methods
are being judged on, so the false-positive rate it reports is not the one the
comparison claims. These tests assert the property directly, on the module
the analysis script calls.
"""

import numpy as np
import pytest

from scpipe import permutation


def _design(n_donors=8, per_donor=30, seed=0):
    rng = np.random.default_rng(seed)
    donors, conds = [], []
    for i in range(n_donors):
        for c in ("ctrl", "stim"):
            donors.extend([f"d{i}"] * per_donor)
            conds.extend([c] * per_donor)
    return np.array(donors), np.array(conds)


def test_permutation_never_moves_a_cell_between_donors():
    """The donor column is unchanged by construction -- this is the core
    property, and it holds for the per-cell path, which is where a careless
    implementation would shuffle the cells themselves."""
    donors, conds = _design()
    for seed in range(10):
        rng = np.random.default_rng(seed)
        swapped, _ = permutation.swap_cells(donors, conds, rng)
        ok, problems = permutation.donor_structure_preserved(
            donors, donors, conds, swapped)
        assert ok, problems


def test_every_donor_keeps_exactly_one_aliquot_of_each_condition():
    donors, conds = _design()
    for seed in range(10):
        rng = np.random.default_rng(seed)
        swapped, _ = permutation.swap_cells(donors, conds, rng)
        ok, problems = permutation.donor_structure_preserved(
            donors, donors, conds, swapped)
        assert ok, problems
        for unit in sorted(set(donors)):
            rows = donors == unit
            before = {c: int((conds[rows] == c).sum()) for c in ("ctrl", "stim")}
            after = {c: int((swapped[rows] == c).sum()) for c in ("ctrl", "stim")}
            assert before == after, (unit, before, after)


def test_swapped_donor_has_its_conditions_exchanged_not_spread():
    """A donor marked for swapping must flip *all* of its cells, so that the
    within-donor contrast is reversed and nothing else changes."""
    donors, conds = _design(n_donors=3, per_donor=5)
    plan = permutation.sample_swap_plan(donors, np.random.default_rng(0))
    assert set(plan) == {"d0", "d1", "d2"}
    swapped, plan2 = permutation.swap_cells(donors, conds,
                                            np.random.default_rng(0))
    assert plan2 == plan
    for unit, do_swap in plan.items():
        rows = donors == unit
        if do_swap:
            flipped = np.where(conds[rows] == "ctrl", "stim", "ctrl")
        else:
            flipped = conds[rows]
        assert np.array_equal(swapped[rows], flipped), unit


def test_swap_samples_flips_a_pair_and_leaves_the_rest():
    donors = np.array(["d0", "d0", "d1", "d1"])
    conds = np.array(["ctrl", "stim", "ctrl", "stim"])
    rng = np.random.default_rng(1)
    out, plan = permutation.swap_samples(donors, conds, rng)
    for unit, do_swap in plan.items():
        rows = donors == unit
        if do_swap:
            assert out[rows][0] != conds[rows][0]
        else:
            assert np.array_equal(out[rows], conds[rows])


def test_swap_fixes_nothing_when_swap_probability_is_zero():
    donors, conds = _design()
    rng = np.random.default_rng(0)
    out, plan = permutation.swap_cells(donors, conds, rng,
                                       swap_probability=0.0)
    assert np.array_equal(out, conds)
    assert not any(plan.values())


def test_structure_check_detects_a_cell_shuffle():
    """The guard has to actually fail on the wrong implementation, otherwise
    the passing tests above prove nothing."""
    donors, conds = _design(n_donors=4, per_donor=5)
    rng = np.random.default_rng(3)
    bad_donors = rng.permutation(donors)
    ok, problems = permutation.donor_structure_preserved(
        donors, bad_donors, conds, conds)
    assert not ok
    assert "changed donor" in problems[0]


def test_structure_check_detects_an_unequal_relabelling():
    donors = np.array(["d0", "d0", "d1", "d1"])
    conds = np.array(["ctrl", "stim", "ctrl", "stim"])
    bad = np.array(["stim", "stim", "ctrl", "stim"])   # d0 lost its ctrl
    ok, problems = permutation.donor_structure_preserved(
        donors, donors, conds, bad)
    assert not ok
    assert "labels" in problems[0]


def test_swap_on_a_donor_with_one_aliquot_does_not_invent_one():
    """Swapping a donor that has no ctrl profile would move a cell into a
    condition its donor never had, so it must be left untouched."""
    donors = np.array(["d0", "d1", "d1"])
    conds = np.array(["stim", "ctrl", "stim"])
    rng = np.random.default_rng(0)
    out, _ = permutation.swap_samples(donors, conds, rng,
                                      swap_probability=1.0)
    assert out[0] == "stim"
    assert sorted(out[1:]) == ["ctrl", "stim"]
