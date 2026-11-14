"""Permutation nulls that keep the donor structure intact.

The donor-swap null is the test that decides a method comparison on this
dataset. Its whole validity rests on one property: the permutation must move
*labels*, never *cells*. Each donor supplies one control and one stimulated
aliquot; swapping which of the two is called "stim" makes the true condition
effect exactly zero while leaving every cell inside its own donor, and leaving
the donor-to-donor spread -- the thing the mixed model is supposed to
absorb -- untouched.

The failure this module exists to prevent is the easy one to make by
accident: shuffling conditions across *cells* rather than across *samples*.
That destroys the donor structure as well as the condition effect, so the
null it produces is not the null the method will be judged on, and a
per-cell shuffle is precisely the mistake the rest of this repo is about.
It still looks like a null, which is why the property is pinned by a test
(tests/test_permutation.py) rather than by a comment.

Two functions, two granularities:

    swap_samples   -- one swap decision per donor, applied to all of that
                      donor's samples. This is the null for pseudobulk and
                      for the mixed model, both of which take one profile or
                      one group of cells per (donor, condition).
    swap_cells     -- per-cell version for models fitted to individual cells:
                      every cell is relabelled by its own donor's swap
                      decision. Cells never change donor.

Both return the permuted condition labels in the same order as the input, so
a caller can zip them straight back onto the rows.

A third decision rule, `balanced_swap_plan`, exists because the independent
0.5 coin is not a null when a real effect is present. Flipping each donor
with probability one half gives a Binomial(n, 1/2) number of flips, so a draw
can leave most donors unflipped; those donors keep their real condition
effect, the paired contrast no longer averages to zero, and a method that
correctly finds the effect gets scored as a false positive. Fixing the number
of flips at the balanced count makes the flipped and unflipped halves cancel
*exactly* in expectation, which is what "the true effect is zero" requires.
"""

import numpy as np
import pandas as pd


def sample_swap_plan(samples, rng, swap_probability=0.5):
    """Decide, per sample unit, whether to swap. Returns {unit: bool}.

    samples: array-like of the sample identifier (the donor). Each unit gets
        one decision, and the decision applies to every row belonging to it,
        so a two-aliquot donor ends up with its aliquots exchanged.
    """
    units = sorted(set(np.asarray(samples).astype(str)))
    draws = rng.random(len(units)) < swap_probability
    return dict(zip(units, draws.tolist()))


def balanced_swap_plan(samples, rng, n_swaps=None):
    """Decide exactly `n_swaps` units to flip, and leave the rest alone.

    Default n_swaps is len(units) // 2. When the number of units is even this
    halves the sample into flipped and unflipped, so a constant real effect of
    +e becomes +e on half the donors and -e on the other half and cancels
    exactly. With an odd count one unit is always left over, so a caller that
    needs an exact cancellation restricts to an even subset first (see
    scripts/mixed_vs_others.py, which holds one donor out per draw).

    This is not interchangeable with swap_probability=0.5: the coin version has
    a binomial spread, and an unlucky draw that flips two donors out of eight
    leaves six donors carrying the real effect.
    """
    units = sorted(set(np.asarray(samples).astype(str)))
    n = len(units)
    k = n // 2 if n_swaps is None else int(n_swaps)
    k = max(0, min(k, n))
    if n == 0:
        return {}
    if k == 0:
        flip = set()
    elif k == n:
        flip = set(range(n))
    else:
        flip = {int(i) for i in np.atleast_1d(
            rng.choice(n, size=k, replace=False))}
    return {unit: (i in flip) for i, unit in enumerate(units)}


def _apply_plan(samples, conditions, plan, group_a="stim", group_b="ctrl"):
    """Relabel rows whose unit was flipped, exchanging the two conditions."""
    samples = np.asarray(samples).astype(str)
    conditions = np.asarray(conditions).astype(str)
    out = conditions.copy()
    for unit, do_swap in plan.items():
        if not do_swap:
            continue
        rows = np.flatnonzero(samples == unit)
        labels = conditions[rows]
        counts = {lab: int((labels == lab).sum()) for lab in (group_a, group_b)}
        # exactly one of each: an unambiguous two-aliquot pair
        if counts[group_a] == 1 and counts[group_b] == 1:
            out[rows] = np.where(labels == group_a, group_b, group_a)
        # anything else (duplicated aliquots, a missing condition) is left
        # alone rather than half-swapped
    return out


def swap_samples(samples, conditions, rng, group_a="stim", group_b="ctrl",
                 other_labels=(), swap_probability=0.5, plan=None):
    """Swap the labels of each sample unit's two aliquots.

    samples: one entry per row (pseudobulk profile), the donor.
    conditions: the row's condition label.
    plan: an explicit {unit: flip} decision, for callers that need a balanced
        draw (see `balanced_swap_plan`). When omitted, each unit gets an
        independent coin with probability `swap_probability`.
    Returns a new array of labels, same length and order.

    A unit whose rows are not exactly one group_a and one group_b is left
    untouched: there is nothing to swap, and inventing a swap would move a
    donor into a condition it never had.
    """
    samples = np.asarray(samples).astype(str)
    conditions = np.asarray(conditions).astype(str)
    if plan is None:
        plan = sample_swap_plan(samples, rng, swap_probability)
    out = _apply_plan(samples, conditions, plan, group_a, group_b)
    return out, plan


def swap_cells(samples, conditions, rng, group_a="stim", group_b="ctrl",
               swap_probability=0.5, plan=None):
    """Per-cell labels, swapped by the cell's own donor's decision.

    Every cell keeps its donor. A donor that is swapped has all of its cells
    relabelled, which is what "the two aliquots were exchanged" means at
    single-cell resolution.
    """
    samples = np.asarray(samples).astype(str)
    conditions = np.asarray(conditions).astype(str)
    if plan is None:
        plan = sample_swap_plan(samples, rng, swap_probability)
    do_swap = np.array([plan[s] for s in samples], dtype=bool)
    is_a = conditions == group_a
    is_b = conditions == group_b
    out = conditions.copy()
    out[do_swap & is_a] = group_b
    out[do_swap & is_b] = group_a
    return out, plan


def donor_structure_preserved(samples_before, samples_after,
                              conditions_before, conditions_after):
    """Check the two properties a valid label permutation must have.

    1. No cell/sample changes donor: the donor column is identical.
    2. Within each donor, the multiset of condition labels is unchanged --
       a swap exchanges labels, it does not create or destroy an aliquot.

    Returns (ok, problems) where problems names the first violation found.
    """
    samples_before = np.asarray(samples_before).astype(str)
    samples_after = np.asarray(samples_after).astype(str)
    conditions_before = np.asarray(conditions_before).astype(str)
    conditions_after = np.asarray(conditions_after).astype(str)

    problems = []
    if not np.array_equal(samples_before, samples_after):
        moved = np.flatnonzero(samples_before != samples_after)
        problems.append(
            f"{len(moved)} rows changed donor (first: "
            f"{samples_before[moved[0]]} -> {samples_after[moved[0]]})")
    for unit in sorted(set(samples_before)):
        rows = samples_before == unit
        before = sorted(conditions_before[rows].tolist())
        after = sorted(conditions_after[rows].tolist())
        if before != after:
            problems.append(f"{unit}: labels {before} became {after}")
            break
    return (not problems), problems


def swap_table(samples, conditions, rng, **kwargs):
    """Convenience: swap and return one frame recording before/after.

    Kept small on purpose -- the audit trail for the null lives in the
    script's printed output and CSV, not in a new file format.
    """
    swapped, plan = swap_samples(samples, conditions, rng, **kwargs)
    return pd.DataFrame({
        "sample": np.asarray(samples).astype(str),
        "condition_before": np.asarray(conditions).astype(str),
        "condition_after": swapped,
    }), plan
