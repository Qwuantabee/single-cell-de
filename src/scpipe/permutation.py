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


def swap_samples(samples, conditions, rng, group_a="stim", group_b="ctrl",
                 other_labels=(), swap_probability=0.5):
    """Swap the labels of each sample unit's two aliquots.

    samples: one entry per row (pseudobulk profile), the donor.
    conditions: the row's condition label.
    Returns a new array of labels, same length and order.

    A unit whose rows are not exactly one group_a and one group_b is left
    untouched: there is nothing to swap, and inventing a swap would move a
    donor into a condition it never had.
    """
    samples = np.asarray(samples).astype(str)
    conditions = np.asarray(conditions).astype(str)
    plan = sample_swap_plan(samples, rng, swap_probability)

    out = conditions.copy()
    for unit, do_swap in plan.items():
        rows = np.flatnonzero(samples == unit)
        labels = conditions[rows]
        counts = {lab: int((labels == lab).sum()) for lab in (group_a, group_b)}
        if not do_swap:
            continue
        # exactly one of each: an unambiguous two-aliquot pair
        if counts[group_a] == 1 and counts[group_b] == 1:
            out[rows] = np.where(labels == group_a, group_b, group_a)
        # anything else (duplicated aliquots, a missing condition) is left
        # alone rather than half-swapped
    return out, plan


def swap_cells(samples, conditions, rng, group_a="stim", group_b="ctrl",
               swap_probability=0.5):
    """Per-cell labels, swapped by the cell's own donor's decision.

    Every cell keeps its donor. A donor that is swapped has all of its cells
    relabelled, which is what "the two aliquots were exchanged" means at
    single-cell resolution.
    """
    samples = np.asarray(samples).astype(str)
    conditions = np.asarray(conditions).astype(str)
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
