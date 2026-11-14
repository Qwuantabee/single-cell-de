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


# --------------------------------------------------- distinct balanced plans


def test_distinct_plan_counts_match_the_combinatorics():
    # first unit fixed unflipped: C(n - 1, n / 2); a plan and its complement
    # are the same draw
    assert permutation.n_distinct_balanced_plans(4) == 3
    assert permutation.n_distinct_balanced_plans(6) == 10
    assert permutation.n_distinct_balanced_plans(8) == 35
    with pytest.raises(ValueError):
        permutation.n_distinct_balanced_plans(5)


def test_four_donors_have_only_three_distinct_balanced_draws():
    """The failure being fixed: sampling plans with replacement on a 4-donor
    panel gave repeated draws (21;21;21) that looked like replicates."""
    units = [f"d{i}" for i in range(4)]
    plans = permutation.distinct_balanced_plans(units)
    assert len(plans) == 3
    keys = {tuple(sorted(u for u, f in p.items() if f)) for p in plans}
    assert len(keys) == 3
    for p in plans:
        assert sum(p.values()) == 2
        assert p["d0"] is False        # complements collapsed onto one rep.
    # complements really are the same draw: no two plans are complementary
    for a in keys:
        comp = tuple(sorted(set(units) - set(a)))
        assert comp not in keys or comp == a


def test_distinct_plans_sampled_without_replacement_when_capped():
    units = [f"d{i}" for i in range(8)]
    plans = permutation.distinct_balanced_plans(
        units, np.random.default_rng(0), max_plans=20)
    keys = [tuple(sorted(u for u, f in p.items() if f)) for p in plans]
    assert len(keys) == 20 and len(set(keys)) == 20
    assert all(len(k) == 4 for k in keys)
    with pytest.raises(ValueError):
        permutation.distinct_balanced_plans(units, None, max_plans=5)


def test_balanced_plan_specs_enumerates_every_odd_panel_draw_once():
    """Five donors: hold one out (5 ways) x 3 plans on the other four = 15."""
    units = [f"d{i}" for i in range(5)]
    specs = permutation.balanced_plan_specs(
        units, np.random.default_rng(0), max_draws=200)
    assert permutation.n_plans_available(5) == 15
    assert len(specs) == 15
    keys = {(h, tuple(sorted(u for u, f in p.items() if f)))
            for h, _used, p in specs}
    assert len(keys) == 15                    # no repeats
    for held, used, plan in specs:
        assert held in units and held not in used and len(used) == 4
        assert set(plan) == set(used)         # the held-out donor is not in it
        assert sum(plan.values()) == 2


def test_balanced_plan_specs_cap_and_availability():
    assert permutation.n_plans_available(8) == 35
    assert permutation.n_plans_available(7) == 70
    assert permutation.n_plans_available(3) == 0
    units = [f"d{i}" for i in range(7)]
    specs = permutation.balanced_plan_specs(
        units, np.random.default_rng(1), max_draws=10)
    assert len(specs) == 10
    keys = {(h, tuple(sorted(u for u, f in p.items() if f)))
            for h, _u, p in specs}
    assert len(keys) == 10


def test_wilson_interval_and_bootstrap_ci_basic_properties():
    lo, hi = permutation.wilson_interval(0, 70)
    assert lo == 0.0 and 0.04 < hi < 0.07          # 0/70 does not rule out 5%
    lo, hi = permutation.wilson_interval(35, 70)
    assert lo < 0.5 < hi
    lo, hi = permutation.bootstrap_ci([0, 0, 0, 0], np.random.default_rng(0))
    assert (lo, hi) == (0.0, 0.0)
    lo, hi = permutation.bootstrap_ci([1, 2, 3, 4, 5],
                                      np.random.default_rng(0))
    assert lo < 3.0 < hi
