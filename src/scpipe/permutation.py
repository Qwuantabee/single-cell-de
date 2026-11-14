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
for a *constant* real effect, which is what "the true effect is zero"
requires.

What balance does NOT do: if the real effect differs between donors (a gene
that responds +2 in one donor and +0.5 in another), the flipped half and the
unflipped half no longer cancel in variance. The swapped contrast has mean
zero but carries the real effect's between-donor spread as extra variance, so
a test that is exact for iid noise sees *more* spread than it models only in
the direction that makes it conservative. The swap null therefore cannot show
that a procedure is calibrated or not over-conservative, only that it is not
anti-conservative. A caller that wants a cleaner null restricts the gene set
to genes with no measured effect, as scripts/nbglm_real.py does.
"""

import itertools

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
    halves the sample into flipped and unflipped, so a CONSTANT real effect of
    +e becomes +e on half the donors and -e on the other half and cancels in
    the mean. (An effect that varies between donors does not cancel in
    variance: see the module docstring. The swap null is conservative, not
    exact, for such genes.) With an odd count one unit is always left over,
    so a caller that needs the cancellation restricts to an even subset first
    (see scripts/mixed_vs_others.py, which holds one donor out per draw), and
    must record which donors were used.

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


def n_distinct_balanced_plans(n_units):
    """Distinct balanced plans on an even number of units, up to complement.

    A plan and its complement give the same |t| for any paired statistic (the
    contrast changes sign), so they are one draw, not two. Fixing the first
    unit unflipped leaves C(n - 1, n / 2) plans. n = 4 -> 3, 6 -> 10,
    8 -> 35. This is the ceiling on independent balanced draws: a 4-donor
    panel cannot yield 200 distinct ones.
    """
    if n_units < 2 or n_units % 2:
        raise ValueError(f"need an even number of units >= 2, got {n_units}")
    from math import comb
    return comb(n_units - 1, n_units // 2)


def distinct_balanced_plans(samples, rng=None, max_plans=None):
    """Distinct balanced swap plans, one per complement pair.

    Enumerates every plan (first unit unflipped, n/2 of the rest flipped)
    when there are at most `max_plans` of them; otherwise samples
    `max_plans` of them WITHOUT replacement. Sampling with replacement is
    what produced repeated draws (21;21;21) on a 4-donor panel, where only
    three distinct balanced plans exist.

    samples: the unit identifiers; must be an even count.
    Returns a list of {unit: bool} dicts, all with exactly n/2 flips.
    """
    units = sorted(set(np.asarray(samples).astype(str)))
    n = len(units)
    total = n_distinct_balanced_plans(n)      # raises on odd / tiny n
    k = n // 2
    rest = list(range(1, n))
    if max_plans is None or max_plans >= total:
        combos = list(itertools.combinations(rest, k))
    else:
        if rng is None:
            raise ValueError("rng is required to sample a subset of plans")
        seen, combos = set(), []
        while len(combos) < max_plans:
            pick = tuple(sorted(int(i) for i in rng.choice(
                rest, size=k, replace=False)))
            if pick not in seen:
                seen.add(pick)
                combos.append(pick)
    return [{u: (i in set(c)) for i, u in enumerate(units)} for c in combos]


def n_plans_available(n_units):
    """How many distinct balanced draws exist for n paired donors.

    Even n: C(n - 1, n / 2) (see n_distinct_balanced_plans). Odd n: one donor
    is held out per draw, so each (held-out donor, plan on the other n - 1)
    pair is a distinct draw: n * C(n - 2, (n - 1) / 2). n = 5 -> 15,
    7 -> 70, 8 -> 35. Fewer than four donors cannot be balanced: 0.
    """
    from math import comb
    if n_units < 4:
        return 0
    if n_units % 2 == 0:
        return comb(n_units - 1, n_units // 2)
    return n_units * comb(n_units - 2, (n_units - 1) // 2)


def balanced_plan_specs(units, rng, max_draws=200):
    """Distinct balanced draws as (held_out, used, plan) triples.

    Even panel: one draw per complement pair of balanced plans, held_out is
    None. Odd panel: every (held-out donor, plan on the remaining even set)
    pair is a distinct draw. At most `max_draws` are returned, sampled
    WITHOUT replacement when there are more; fewer when the plan space is
    smaller. A draw is never repeated, and each records the donors it used.
    """
    units = sorted(set(np.asarray(units).astype(str).tolist()))
    if len(units) < 4:
        raise ValueError(f"need >=4 donors for a balanced swap, got {len(units)}")
    if len(units) % 2 == 0:
        plans = distinct_balanced_plans(units, rng, max_draws)
        return [(None, list(units), p) for p in plans]
    specs = []
    for held in units:
        used = [u for u in units if u != held]
        for plan in distinct_balanced_plans(used):
            specs.append((held, used, plan))
    if len(specs) > max_draws:
        pick = rng.choice(len(specs), size=max_draws, replace=False)
        specs = [specs[int(i)] for i in sorted(pick)]
    else:
        order = rng.permutation(len(specs))
        specs = [specs[int(i)] for i in order]
    return specs


def wilson_interval(k, n, z=1.96):
    """Wilson score interval for a binomial proportion k / n."""
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return float(max(0.0, centre - half)), float(min(1.0, centre + half))


def bootstrap_ci(values, rng=None, n_boot=5000, level=0.95):
    """Percentile bootstrap CI for the mean over draws.

    Resamples the draws. It reflects variation across swap plans for a FIXED
    gene set; it does not include gene-sampling variation, and when the draws
    are an exhaustive enumeration of the plan space it is the spread of that
    space, not a sampling error.
    """
    v = np.asarray(values, dtype=float)
    if len(v) == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(0) if rng is None else rng
    means = v[rng.integers(0, len(v), size=(n_boot, len(v)))].mean(axis=1)
    lo = (1 - level) / 2
    return (float(np.quantile(means, lo)), float(np.quantile(means, 1 - lo)))


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