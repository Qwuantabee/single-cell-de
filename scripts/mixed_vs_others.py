"""Item 3: mixed-effects model vs per-cell and pseudobulk, on every cell type.

Three models are fitted to each cell-type subset of Kang 2018:

  per-cell     : Wilcoxon on cells, every cell a replicate (wrong, for scale)
  pseudobulk   : paired t-test on donor x condition sums (7-8 df)
  mixed        : expression ~ condition + (1|donor), fitted per gene

The first version of this script hardcoded CD4 T cells, and the CD4 result
was then described as a property of the method. It is a property of one
subset until the others are run, so this version loops over all seven
testable types (Megakaryocytes are excluded: 132 cells, and after the
min-cells filter they do not supply a donor panel -- see
scripts/megakaryocyte_check.py).

The mixed model is the one being evaluated, so it gets the harder questions:
does it converge, what does it cost, and does it agree with the method that
is already trusted on this data.

Cost is measured, not estimated. The per-cell fit is O(genes x iterative
optimisation) and is run on a deliberate subsample so the wall-clock is
recorded rather than extrapolated.

The donor-swap null is run BALANCED. Flipping each donor with a 0.5 coin
leaves a Binomial(n, 1/2) number of flips, so some draws leave most donors
carrying their real condition effect; a method that finds that real effect is
then scored as a false positive. A balanced draw flips exactly half the
donors, so the flipped and unflipped halves cancel for a CONSTANT real effect
and the truth is "nothing". The unbalanced version is still measured and
reported next to it, because the difference between the two is the whole
reason the balance matters -- but it is not a null and its numbers are not the
headline (its rows carry is_null = False and no expected-FP ratio).

What the swap null can and cannot show. A gene whose real effect differs
between donors does not cancel in variance, so the swapped contrast carries
that spread as extra noise and a test sees MORE variance than it models: the
swap null is conservative for such genes. It can therefore show that a method
is not anti-conservative, but "0 false positives in every draw" cannot show
that a method is calibrated (P(0 FP in 70 draws | exact calibration) is
0.95^70 = 2.8% only if every draw were an independent exact null, which they
are not). Two things are done about it: the gene set is the small-effect set
that scripts/nbglm_real.py uses (|paired log2FC| < 0.5, ISGs excluded) rather
than a random sample of all tested genes, and the raw per-gene rate of
p < 0.05 is reported next to the BH count, because a conservative null shows
there (well under 5%) while BH counts near zero are uninformative.

Draws are DISTINCT balanced plans (one per complement pair) enumerated or
sampled without replacement, never resampled: a 4-donor panel has three
distinct balanced plans, so repeated draws on it were repeats, not
replicates. With an odd donor count one donor is held out per draw and each
(held-out donor, plan) pair is a distinct draw. The number of plans
available is recorded; for the 5-donor types it is 15, so 200 independent
draws do not exist there. Every draw records the donors it used.

Run:
    python scripts/mixed_vs_others.py                      # all 7 types
    python scripts/mixed_vs_others.py --types "CD4 T cells,B cells"
    python scripts/mixed_vs_others.py --combine            # table from parts
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from scpipe import data, de, mixed, permutation, pseudobulk  # noqa: E402

RESULTS = ROOT / "results"
PARTS = RESULTS / "mixed_per_type"
RESULTS.mkdir(exist_ok=True)
PARTS.mkdir(exist_ok=True)

ALPHA = 0.05
MIN_CELLS = 25

# The seven testable types. Megakaryocytes are deliberately absent: too few
# cells per donor to build 8 paired profiles, so they cannot enter a paired
# or donor-random-effect comparison at all.
ALL_TYPES = [
    "CD4 T cells",
    "CD14+ Monocytes",
    "B cells",
    "NK cells",
    "CD8 T cells",
    "FCGR3A+ Monocytes",
    "Dendritic cells",
]

# expected false positives per gene set at BH q=0.05 on all-null data.
# NOT alpha * m (that is the uncorrected expectation). BH bounds the
# false-discovery proportion among rejections, so on a family of independent
# exact nulls P(any rejection) = q = 0.05 and the mean count is slightly
# above that. Simulated by scripts/bh_null_expectation.py (5,000 draws,
# m = 500: results/bh_null_expectation.csv); 0.05 is used as the round
# reference and tests/test_bh_null_expectation.py pins it to that simulation.
# It is the reference for independent exact nulls, not a derived bound for the
# swap null below (fixed gene set, correlated genes).
BH_EXPECTED_FP = 0.05

# per-cell fits are the expensive arm: cap them deliberately
PERCELL_GENES = 400
PERCELL_MAX_CELLS = 3200
NULL_GENES = 500
# upper bound on DISTINCT balanced draws; the plan space is smaller than this
# for every real panel here (n=5: 15, n=7: 70, n=8: 35), so it is exhaustive
NULL_PERMS = 200
UNBALANCED_PERMS = 10     # contrast only, not a null
NULL_FC_CUTOFF = 0.5      # same small-effect cutoff as scripts/nbglm_real.py
SEED = 0


def slug(cell_type):
    """File-safe stem for a cell type ('CD14+ Monocytes' -> 'CD14_Monocytes')."""
    return (cell_type.replace("+", "")
            .replace(" ", "_")
            .replace("-", "_")
            .replace("__", "_"))


def _subsample(adata, n_max, seed):
    """Take cells evenly across donor x condition, never at random globally.

    A uniform random subsample would let a donor with many cells dominate
    the fit; stratifying preserves the paired design, which is the only
    thing the random effect is estimated from.
    """
    rng = np.random.default_rng(seed)
    obs = adata.obs
    keys = obs["donor"].astype(str) + "||" + obs["condition"].astype(str)
    per = max(1, n_max // keys.nunique())
    idx = []
    for _, group in keys.groupby(keys).groups.items():
        g = np.asarray(group)
        idx.extend(rng.choice(g, size=min(per, len(g)), replace=False).tolist())
    return np.sort(np.asarray(idx))


def paired_units(meta):
    """Donors that contribute exactly one sample per condition.

    Only these can supply a paired contrast or a balanced swap; a donor with
    two control profiles has nothing to pair, and one with a single aliquot
    cannot be flipped without inventing a profile.
    """
    out = []
    for donor, rows in meta.groupby("donor", observed=True).groups.items():
        sub = meta.loc[rows, "condition"]
        if (sub == "stim").sum() == 1 and (sub == "ctrl").sum() == 1:
            out.append(str(donor))
    return sorted(out)


def balanced_swap(units, rng):
    """Balanced donor swap: flip exactly half of an even donor set.

    Odd donor counts cannot be balanced -- k = n/2 is not an integer, and
    rounding leaves one donor carrying its real effect on every draw. One
    donor is therefore held out at random, and the returned `used` list names
    the donors the null was actually computed on. Holding out a donor costs
    one degree of freedom and buys an exact cancellation, which is the
    property the null is supposed to have.
    """
    units = sorted(set(np.asarray(units).astype(str).tolist()))
    if len(units) < 4:
        raise ValueError(f"need >=4 donors for a balanced swap, got {len(units)}")
    held_out = None
    if len(units) % 2:
        held_out = str(rng.choice(units))
        units = [u for u in units if u != held_out]
    plan = permutation.balanced_swap_plan(units, rng)
    return plan, units, held_out


balanced_plan_specs = permutation.balanced_plan_specs
n_plans_available = permutation.n_plans_available
wilson_interval = permutation.wilson_interval
bootstrap_ci = permutation.bootstrap_ci


def null_gene_pool(pb, columns, markers=None, cutoff=NULL_FC_CUTOFF):
    """Genes the swap null can treat as null: small observed effect, no ISG.

    Same rule as nbglm_real.null_gene_set. The selection uses the observed
    paired log2FC of the very data being permuted, which is a mild selection
    effect (documented, and the reason the unrestricted random set is also
    reported for pseudobulk).
    """
    markers = set(data.ISG_MARKERS if markers is None else markers)
    common = [g for g in columns if g in pb.index]
    return [g for g in common
            if abs(pb.loc[g, "log2fc"]) < cutoff and g not in markers]


def paired_columns(logcpm, meta, labels, units, genes, group_a="stim",
                   group_b="ctrl"):
    """Pull the (condition A, condition B) rows of each donor, aligned by donor.

    Returns (A, B) as arrays with one row per unit, ordered identically in
    both, which is what a paired test requires.
    """
    labels = np.asarray(labels).astype(str)
    donor_arr = meta["donor"].astype(str).to_numpy()
    a_rows, b_rows = [], []
    for unit in units:
        own = np.flatnonzero(donor_arr == unit)
        if len(own) != 2:
            raise ValueError(f"{unit}: expected 2 profiles, got {len(own)}")
        ia = own[labels[own] == group_a]
        ib = own[labels[own] == group_b]
        if len(ia) != 1 or len(ib) != 1:
            raise ValueError(f"{unit}: swap left {len(ia)}/{len(ib)}")
        a_rows.append(int(ia[0]))
        b_rows.append(int(ib[0]))
    A = logcpm.iloc[a_rows][genes].to_numpy(dtype=float)
    B = logcpm.iloc[b_rows][genes].to_numpy(dtype=float)
    return A, B


def null_draw_detail(logcpm, meta, genes, rng, balanced=True, units=None,
                     spec=None, extra_genes=None):
    """One donor-swap null draw, with everything needed to audit it.

    The permutation moves labels only: `scpipe.permutation` re-checks the donor
    structure on every draw, and a violation raises rather than being printed
    and ignored, because a broken permutation silently invalidates the number
    it produces.

    spec: an explicit (held_out, used, plan) from balanced_plan_specs; when
        omitted a random balanced (or unbalanced) plan is drawn.
    extra_genes: optional second gene list scored with the pseudobulk arm only
        (cheap), used to contrast the small-effect set with a random set.

    Returns a dict: fp_pb, fp_mx (BH calls at ALPHA), raw_pb, raw_mx (the
    fraction of genes with raw p < ALPHA), used, held_out, flipped, boundary
    (fraction of mixed fits on the variance boundary), structure_ok.
    """
    donors = meta["donor"].astype(str).to_numpy()
    cond = meta["condition"].astype(str).to_numpy()
    units = paired_units(meta) if units is None else list(units)
    held_out = None
    if spec is not None:
        held_out, used, plan = spec
        used = list(used)
    elif balanced:
        plan, used, held_out = balanced_swap(units, rng)
    else:
        # unbalanced: independent coin on the paired donors only, so the two
        # arms differ in the balance of the flips and nothing else
        plan = permutation.sample_swap_plan(units, rng, 0.5)
        used = list(units)

    flipped, _ = permutation.swap_samples(donors, cond, rng, plan=plan)
    flipped_units = sorted(u for u in used if plan.get(u))
    if balanced and len(flipped_units) != len(used) // 2:
        raise RuntimeError(
            f"balanced draw flipped {len(flipped_units)} of {len(used)}")

    # ---- pseudobulk: paired t on the swapped labels ---------------------
    A, B = paired_columns(logcpm, meta, flipped, used, genes)
    _, p_t = stats.ttest_rel(A, B, axis=0)
    fp_pb = int((de._bh(p_t) < ALPHA).sum())
    raw_pb = float(np.mean(p_t < ALPHA))
    fp_pb_extra = None
    if extra_genes is not None:
        A2, B2 = paired_columns(logcpm, meta, flipped, used, extra_genes)
        _, p2 = stats.ttest_rel(A2, B2, axis=0)
        fp_pb_extra = int((de._bh(p2) < ALPHA).sum())

    # ---- mixed: same swapped labels, through the module -----------------
    rows = np.flatnonzero(np.isin(donors, used))
    permuted = logcpm.iloc[rows][genes].copy()
    permuted.index = pd.Index(
        [f"{d}||{c}" for d, c in zip(donors[rows], flipped[rows])])

    # donor-structure check on what the model is actually fed: the donor and
    # condition labels are parsed back out of the permuted frame's own index,
    # not reused from the arrays that built it, and the values must be the
    # original rows untouched (labels moved, never cells).
    after_d, after_c = mixed._split_index(permuted.index, "donor", "condition")
    ok, problems = permutation.donor_structure_preserved(
        donors[rows], after_d, cond[rows], after_c)
    if not ok:
        raise RuntimeError(f"donor structure broken: {problems}")
    if not np.array_equal(permuted.to_numpy(),
                          logcpm.iloc[rows][genes].to_numpy()):
        raise RuntimeError("expression values changed during the swap")

    pm = mixed.mixed_de(permuted, "donor", "condition", "stim", "ctrl",
                        donors=used)
    fp_mx = int((de._bh(pm["pval"].to_numpy()) < ALPHA).sum())
    raw_mx = float(np.nanmean(pm["pval"].to_numpy() < ALPHA))
    return {
        "fp_pb": fp_pb, "fp_mx": fp_mx, "raw_pb": raw_pb, "raw_mx": raw_mx,
        "fp_pb_extra": fp_pb_extra, "used": list(used), "held_out": held_out,
        "flipped": flipped_units,
        "boundary": float(pm["boundary"].mean()),
        "structure_ok": bool(ok),
    }


def null_draw(logcpm, meta, genes, rng, balanced=True, units=None):
    """One donor-swap null draw; returns (fp_pseudobulk, fp_mixed, used)."""
    d = null_draw_detail(logcpm, meta, genes, rng, balanced=balanced,
                         units=units)
    return d["fp_pb"], d["fp_mx"], d["used"]


def summarise_draws(draws, rng=None):
    """Mean FP with bootstrap CI, P(any FP) with Wilson CI, raw p<alpha rate.

    draws: list of per-draw detail dicts for ONE balance family.
    Returns {method: {...}} for "pseudobulk" and "mixed".
    """
    out = {}
    for method, fp_key, raw_key in (("pseudobulk", "fp_pb", "raw_pb"),
                                    ("mixed", "fp_mx", "raw_mx")):
        fp = np.array([d[fp_key] for d in draws], dtype=float)
        raw = np.array([d[raw_key] for d in draws], dtype=float)
        n_any = int((fp > 0).sum())
        lo, hi = bootstrap_ci(fp, rng)
        rlo, rhi = bootstrap_ci(raw, rng)
        alo, ahi = wilson_interval(n_any, len(fp))
        out[method] = {
            "draws": len(fp), "mean_false_positives": float(fp.mean()),
            "fp_ci_lo": lo, "fp_ci_hi": hi,
            "n_draws_with_any_fp": n_any,
            "p_any": n_any / len(fp), "p_any_lo": alo, "p_any_hi": ahi,
            "raw_p05_rate": float(raw.mean()),
            "raw_ci_lo": rlo, "raw_ci_hi": rhi,
            "per_draw": ";".join(str(int(v)) for v in fp),
        }
    return out


def build_null_rows(cell_type, gene_set, n_paired, n_genes_null,
                    n_distinct, n_avail, families):
    """Null summary rows and per-draw rows, for the CSVs.

    families: iterable of (balance, summary, draws, is_null). Only a real
    (balanced) null carries expected_fp_bh / ratio_observed_expected; the
    unbalanced contrast arm gets NaN there and is_null = False, so a ratio of
    up to hundreds can no longer be read off a row that is not a null.
    `permutation_preserved_donors` is the AND of every draw's own structure
    check, not a constant.
    """
    null_rows = []
    draw_rows = []
    for balance, summ, draws, is_null in families:
        struct_ok = all(dd["structure_ok"] for dd in draws)
        for method in ("pseudobulk", "mixed"):
            r = summ[method]
            mean_fp = r["mean_false_positives"]
            null_rows.append({
                "cell_type": cell_type, "method": method, "balance": balance,
                "is_null": is_null, "gene_set": gene_set,
                "n_paired_donors": n_paired,
                "n_donors_used": len(draws[0]["used"]) if draws else 0,
                "donors_used_per_draw": ";".join(
                    "|".join(dd["used"]) for dd in draws),
                "held_out_per_draw": ";".join(
                    str(dd["held_out"]) for dd in draws),
                "perms": r["draws"],
                "n_distinct_draws": (n_distinct if is_null else np.nan),
                "n_plans_available": (n_avail if is_null else np.nan),
                "genes_tested": n_genes_null,
                "mean_false_positives": mean_fp,
                "fp_ci_lo": r["fp_ci_lo"], "fp_ci_hi": r["fp_ci_hi"],
                "n_draws_with_any_fp": r["n_draws_with_any_fp"],
                "p_any": r["p_any"], "p_any_lo": r["p_any_lo"],
                "p_any_hi": r["p_any_hi"],
                "raw_p05_rate": r["raw_p05_rate"],
                "raw_ci_lo": r["raw_ci_lo"], "raw_ci_hi": r["raw_ci_hi"],
                "per_draw": r["per_draw"],
                "expected_fp_bh": (BH_EXPECTED_FP if is_null else np.nan),
                "ratio_observed_expected": (
                    mean_fp / BH_EXPECTED_FP if is_null else np.nan),
                "permutation_preserved_donors": struct_ok,
            })
        for i, dd in enumerate(draws):
            draw_rows.append({
                "cell_type": cell_type, "balance": balance, "draw": i,
                "is_null": is_null, "held_out": dd["held_out"],
                "donors_used": "|".join(dd["used"]),
                "donors_flipped": "|".join(dd["flipped"]),
                "fp_pseudobulk": dd["fp_pb"], "fp_mixed": dd["fp_mx"],
                "raw_p05_pseudobulk": dd["raw_pb"],
                "raw_p05_mixed": dd["raw_mx"],
                "fp_pseudobulk_random_all_genes": dd["fp_pb_extra"],
                "mixed_boundary_frac": dd["boundary"],
                "structure_ok": dd["structure_ok"],
            })
    return null_rows, draw_rows


def compare_frame(pb, mx, per_cell):
    """Agreement statistics shared by the script and its tests."""
    def called(frame, col="qval"):
        return set(frame.index[frame[col] < ALPHA])

    both = pd.concat([pb["log2fc"].rename("pb"), mx["log2fc"].rename("mx")],
                     axis=1).dropna()
    r_lfc = stats.pearsonr(both["pb"], both["mx"]) if len(both) > 2 else None
    rho = stats.spearmanr(both["pb"], both["mx"]) if len(both) > 2 else None

    def jaccard(a, b):
        return len(a & b) / len(a | b) if (a | b) else np.nan

    s_pc, s_pb, s_mx = called(per_cell), called(pb), called(mx)
    return {
        "percell_sig": len(s_pc),
        "pseudobulk_sig": len(s_pb),
        "mixed_sig": len(s_mx),
        "jaccard_pb_mx": jaccard(s_pb, s_mx),
        "jaccard_percell_mx": jaccard(s_pc, s_mx),
        "lfc_pearson_pb_mx": float(r_lfc[0]) if r_lfc else np.nan,
        "lfc_spearman_pb_mx": float(rho.statistic) if rho else np.nan,
        "n_genes_compared": len(both),
    }


def run_type(adata, cell_type, null_perms=NULL_PERMS, null_genes=NULL_GENES,
             seed=SEED, verbose=True, unbalanced_perms=UNBALANCED_PERMS,
             null_gene_set="small_effect"):
    """Every arm for one cell type. Returns the summary row (a dict)."""
    t_start = time.time()
    log = print if verbose else (lambda *a, **k: None)
    sub = adata[adata.obs["cell_type"].astype(str) == cell_type].copy()
    n_donors = int(sub.obs["donor"].nunique())
    log(f"\n=== {cell_type}: {sub.n_obs:,} cells x {sub.n_vars:,} genes, "
        f"{n_donors} donors")

    mat, meta = pseudobulk.aggregate(sub, "donor", "condition",
                                     min_cells=MIN_CELLS)
    mat = pseudobulk.filter_genes(mat, min_count=10, min_samples=3)
    logcpm = pseudobulk.normalise(mat, log=True)
    units = paired_units(meta)
    log(f"pseudobulk: {logcpm.shape[0]} profiles ({len(units)} paired "
        f"donors) x {logcpm.shape[1]:,} genes")

    # ---- arm 1: per-cell Wilcoxon --------------------------------------
    t0 = time.time()
    per_cell = de.single_cell_de(sub, "condition", "stim", "ctrl")
    t_percell = time.time() - t0
    per_cell = per_cell.loc[logcpm.columns]
    log(f"per-cell Wilcoxon: {de.n_significant(per_cell):,} sig, "
        f"{t_percell:.1f}s")

    # ---- arm 2: pseudobulk paired t ------------------------------------
    t0 = time.time()
    pb = de.pseudobulk_de(sub, "condition", "donor", "stim", "ctrl",
                          paired=True, min_cells=MIN_CELLS)
    pb = pb.loc[logcpm.columns]
    t_pb = time.time() - t0
    log(f"pseudobulk paired t: {de.n_significant(pb):,} sig, {t_pb:.3f}s")

    # the paired t on the ALREADY AGGREGATED matrix, which is the honest
    # denominator for "how much slower is a per-gene fit than the vectorised
    # test": t_pb above also contains cell aggregation. Averaged over many
    # repeats because one call is below timer resolution on the small types.
    A_all, B_all = paired_columns(logcpm, meta, meta["condition"].to_numpy(),
                                  units, list(logcpm.columns))
    reps = 20
    t0 = time.perf_counter()
    for _ in range(reps):
        stats.ttest_rel(A_all, B_all, axis=0)
    t_ttest_only = (time.perf_counter() - t0) / reps

    # ---- arm 3: mixed model per gene -----------------------------------
    # Fit on EXACTLY the rows the paired t uses. Before this restriction the
    # mixed arm also absorbed single-condition profiles from donors that
    # cannot be paired, so the two arms were not fitted on the same samples:
    # max |mixed_log2fc - pseudobulk_log2fc| was ~1e-15 on the 8-donor types
    # but 0.65 on CD8, 0.75 on FCGR3A+ and 0.45 on DC -- exactly the three
    # types with fewer than 8 paired donors.
    t0 = time.time()
    paired_rows = [f"{d}||{c}" for d in units for c in ("stim", "ctrl")]
    paired_rows = [r for r in paired_rows if r in logcpm.index]
    mx_input = logcpm.loc[sorted(paired_rows,
                                key=lambda r: (r.split("||")[0], r))]
    mx = mixed.mixed_de(mx_input, "donor", "condition", "stim", "ctrl",
                        donors=units)
    t_mx = time.time() - t0
    n_conv = int(mx["converged"].sum())
    log(f"mixedlm per gene: {de.n_significant(mx):,} sig, {t_mx:.1f}s, "
        f"converged {n_conv:,}/{len(mx):,} "
        f"({100 * n_conv / len(mx):.1f}%)")

    cmp = compare_frame(pb, mx, per_cell)

    # ---- what "converged" does and does not cover (C5) -----------------
    # converged = optimiser termination; boundary = random intercept collapsed
    # to ~0 (a pooled fit). The paired t on the same rows gives the reference
    # statistic: on non-boundary genes the mixed statistic should equal it,
    # on boundary genes it need not.
    t_ref, _ = stats.ttest_rel(A_all, B_all, axis=0)
    t_ref = pd.Series(t_ref, index=logcpm.columns)
    with np.errstate(divide="ignore", invalid="ignore"):
        se_ratio = (mx["stat"].abs() / t_ref.reindex(mx.index).abs()).replace(
            [np.inf, -np.inf], np.nan)
    is_bd = mx["boundary"].astype(bool)
    n_bd = int(is_bd.sum())
    n_conv_interior = int((mx["converged"] & ~is_bd).sum())
    mx_only = set(mx.index[mx["qval"] < ALPHA]) - set(pb.index[pb["qval"] < ALPHA])
    pb_only = set(pb.index[pb["qval"] < ALPHA]) - set(mx.index[mx["qval"] < ALPHA])
    lfc_max_diff = float((pb["log2fc"] - mx["log2fc"]).abs().max())
    log(f"boundary fits: {n_bd:,}/{len(mx):,} ({100 * n_bd / len(mx):.1f}%); "
        f"converged and interior {n_conv_interior:,}; max |stat/t| "
        f"interior {se_ratio[~is_bd].max():.3f} boundary "
        f"{se_ratio[is_bd].max() if n_bd else float('nan'):.3f}")

    isg = [g for g in data.ISG_MARKERS if g in logcpm.columns]
    isg_sig = {name: sum(frame.loc[g, "qval"] < ALPHA for g in isg)
               for name, frame in (("percell", per_cell),
                                   ("pseudobulk", pb), ("mixed", mx))}
    log(f"ISG recovery of {len(isg)}: {isg_sig}")

    # ---- null: balanced donor swap (the null) + unbalanced (for contrast)
    rng = np.random.default_rng(seed)
    all_genes = list(logcpm.columns)
    random_all = list(logcpm.columns[
        rng.choice(len(all_genes), size=min(null_genes, len(all_genes)),
                   replace=False)])
    pool = null_gene_pool(pb, all_genes)
    if null_gene_set == "small_effect":
        null_genes_list = list(np.asarray(pool, dtype=object)[
            rng.choice(len(pool), size=min(null_genes, len(pool)),
                       replace=False)])
    elif null_gene_set == "all":
        null_genes_list = random_all
    else:
        raise ValueError(f"unknown null_gene_set {null_gene_set!r}")
    log(f"null gene set: {null_gene_set}, {len(null_genes_list)} genes "
        f"({len(pool):,} small-effect genes available)")

    specs = balanced_plan_specs(units, rng, null_perms)
    avail = n_plans_available(len(units))
    bal_draws, unb_draws = [], []
    t0 = time.time()
    for held, used, plan in specs:
        draw_rng = np.random.default_rng(int(rng.integers(1 << 31)))
        bal_draws.append(null_draw_detail(
            logcpm, meta, null_genes_list, draw_rng, balanced=True,
            units=units, spec=(held, used, plan), extra_genes=random_all))
    t_bal = time.time() - t0
    t0 = time.time()
    for _ in range(unbalanced_perms):
        draw_rng = np.random.default_rng(int(rng.integers(1 << 31)))
        unb_draws.append(null_draw_detail(
            logcpm, meta, null_genes_list, draw_rng, balanced=False,
            units=units))
    t_unb = time.time() - t0
    n_genes_null = len(null_genes_list)
    bal = summarise_draws(bal_draws, np.random.default_rng(seed + 1))
    unb = summarise_draws(unb_draws, np.random.default_rng(seed + 2))
    n_distinct = len({(d["held_out"], tuple(d["flipped"]), tuple(d["used"]))
                      for d in bal_draws})
    fp_all_pb = float(np.mean([d["fp_pb_extra"] for d in bal_draws]))
    log(f"null: {len(bal_draws)} balanced draws ({n_distinct} distinct, "
        f"{avail} plans exist) in {t_bal:.1f}s, {len(unb_draws)} unbalanced "
        f"in {t_unb:.1f}s")
    for m in ("pseudobulk", "mixed"):
        r = bal[m]
        log(f"  balanced {m:<10} mean {r['mean_false_positives']:.3f} "
            f"[{r['fp_ci_lo']:.3f}, {r['fp_ci_hi']:.3f}] of {n_genes_null}; "
            f"P(any) {r['p_any']:.3f} [{r['p_any_lo']:.3f}, "
            f"{r['p_any_hi']:.3f}]; raw p<{ALPHA} rate "
            f"{r['raw_p05_rate']:.4f}")
    log(f"  unbalanced pseudobulk {unb['pseudobulk']['mean_false_positives']:.2f}"
        f" mixed {unb['mixed']['mean_false_positives']:.2f} (contrast, not a null)")

    # ---- per-cell mixed model: the expensive arm ----------------------
    idx = _subsample(sub, PERCELL_MAX_CELLS, seed)
    cell_sub = sub[idx]
    Xc = cell_sub.X
    Xc = Xc.toarray() if hasattr(Xc, "toarray") else np.asarray(Xc)
    genes_pc = list(logcpm.columns[
        np.random.default_rng(seed).choice(
            len(logcpm.columns),
            size=min(PERCELL_GENES, len(logcpm.columns)), replace=False)])
    col_idx = [cell_sub.var_names.get_loc(g) for g in genes_pc]
    vals = Xc[:, col_idx].astype(float)
    lib = Xc.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    vals = np.log2(vals / lib * 1e6 + 1.0)
    t0 = time.time()
    mx_cell = mixed.mixed_de_percell(
        vals, genes_pc, cell_sub.obs["donor"].to_numpy(),
        cell_sub.obs["condition"].to_numpy(), "stim", "ctrl")
    t_mx_cell = time.time() - t0
    per_gene_cell = t_mx_cell / len(genes_pc)
    log(f"per-cell mixed: {len(idx):,} cells x {len(genes_pc)} genes, "
        f"{t_mx_cell:.1f}s ({1000 * per_gene_cell:.1f} ms/gene), converged "
        f"{100 * mx_cell['converged'].mean():.0f}%")

    # ---- write the per-type artefacts -----------------------------------
    stem = slug(cell_type)
    models = pd.DataFrame({
        "percell_log2fc": per_cell["log2fc"],
        "percell_qval": per_cell["qval"],
        "pb_log2fc": pb["log2fc"],
        "pb_qval": pb["qval"],
        "mx_log2fc": mx["log2fc"],
        "mx_qval": mx["qval"],
        "mx_converged": mx["converged"],
        "mx_boundary": mx["boundary"],
        "mx_donor_var": mx["donor_var"],
        "mx_stat": mx["stat"],
        "paired_t_stat": t_ref.reindex(mx.index),
    })
    models.to_csv(PARTS / f"{stem}_models.csv")

    null_rows, draw_rows = build_null_rows(
        cell_type, null_gene_set, len(units), n_genes_null, n_distinct, avail,
        (("balanced", bal, bal_draws, True),
         ("unbalanced", unb, unb_draws, False)))
    pd.DataFrame(null_rows).to_csv(PARTS / f"{stem}_null.csv", index=False)
    pd.DataFrame(draw_rows).to_csv(PARTS / f"{stem}_null_draws.csv",
                                   index=False)

    row = {
        "cell_type": cell_type,
        "n_cells": int(sub.n_obs),
        "n_donors": n_donors,
        "n_paired_donors": len(units),
        "n_genes_tested": int(logcpm.shape[1]),
        "percell_sig": cmp["percell_sig"],
        "pseudobulk_sig": cmp["pseudobulk_sig"],
        "mixed_sig": cmp["mixed_sig"],
        "mixed_converged": n_conv,
        "mixed_pct_converged": round(100 * n_conv / len(mx), 2),
        "mixed_boundary": n_bd,
        "mixed_pct_boundary": round(100 * n_bd / len(mx), 2),
        "mixed_converged_interior": n_conv_interior,
        "mixed_pct_converged_interior": round(
            100 * n_conv_interior / len(mx), 2),
        "mixed_median_donor_var": float(mx["donor_var"].median()),
        "mixed_stat_over_t_max_interior": float(se_ratio[~is_bd].max()),
        "mixed_stat_over_t_max_boundary": (
            float(se_ratio[is_bd].max()) if n_bd else np.nan),
        "mixed_only_calls": len(mx_only),
        "mixed_only_calls_boundary": int(sum(bool(is_bd.get(g, False))
                                             for g in mx_only)),
        "pseudobulk_only_calls": len(pb_only),
        "lfc_max_abs_diff_pb_mx": lfc_max_diff,
        "jaccard_pb_mx": round(cmp["jaccard_pb_mx"], 4),
        "jaccard_percell_mx": round(cmp["jaccard_percell_mx"], 4),
        "lfc_pearson_pb_mx": round(cmp["lfc_pearson_pb_mx"], 4),
        "lfc_spearman_pb_mx": round(cmp["lfc_spearman_pb_mx"], 4),
        "mixed_s": round(t_mx, 1),
        "mixed_ms_per_gene_pseudobulk": round(1000 * t_mx / len(mx), 2),
        "pseudobulk_s": round(t_pb, 4),
        "paired_t_only_s": round(t_ttest_only, 5),
        "percell_s": round(t_percell, 1),
        "isg_present": len(isg),
        "isg_percell": isg_sig["percell"],
        "isg_pseudobulk": isg_sig["pseudobulk"],
        "isg_mixed": isg_sig["mixed"],
        "null_gene_set": null_gene_set,
        "null_genes": n_genes_null,
        "null_genes_available": len(pool),
        "null_perms": bal["mixed"]["draws"],
        "null_distinct_draws": n_distinct,
        "null_plans_available": avail,
        "null_n_donors_used": len(bal_draws[0]["used"]),
        "null_fp_pb_balanced": round(bal["pseudobulk"]["mean_false_positives"], 3),
        "null_fp_pb_ci_lo": round(bal["pseudobulk"]["fp_ci_lo"], 3),
        "null_fp_pb_ci_hi": round(bal["pseudobulk"]["fp_ci_hi"], 3),
        "null_fp_mx_balanced": round(bal["mixed"]["mean_false_positives"], 3),
        "null_fp_mx_ci_lo": round(bal["mixed"]["fp_ci_lo"], 3),
        "null_fp_mx_ci_hi": round(bal["mixed"]["fp_ci_hi"], 3),
        "null_pany_pb": round(bal["pseudobulk"]["p_any"], 4),
        "null_pany_pb_lo": round(bal["pseudobulk"]["p_any_lo"], 4),
        "null_pany_pb_hi": round(bal["pseudobulk"]["p_any_hi"], 4),
        "null_pany_mx": round(bal["mixed"]["p_any"], 4),
        "null_pany_mx_lo": round(bal["mixed"]["p_any_lo"], 4),
        "null_pany_mx_hi": round(bal["mixed"]["p_any_hi"], 4),
        "null_raw_p05_pb": round(bal["pseudobulk"]["raw_p05_rate"], 4),
        "null_raw_p05_mx": round(bal["mixed"]["raw_p05_rate"], 4),
        "null_fp_pb_random_all_genes": round(fp_all_pb, 3),
        "null_ratio_pb_balanced": round(
            bal["pseudobulk"]["mean_false_positives"] / BH_EXPECTED_FP, 2),
        "null_ratio_mx_balanced": round(
            bal["mixed"]["mean_false_positives"] / BH_EXPECTED_FP, 2),
        "null_fp_pb_unbalanced": round(
            unb["pseudobulk"]["mean_false_positives"], 2),
        "null_fp_mx_unbalanced": round(
            unb["mixed"]["mean_false_positives"], 2),
        "null_unbalanced_perms": len(unb_draws),
        "null_seconds": round(t_bal + t_unb, 1),
        "percell_mixed_cells_fitted": int(len(idx)),
        "percell_mixed_genes_fitted": int(len(genes_pc)),
        "percell_mixed_s": round(t_mx_cell, 1),
        "percell_mixed_ms_per_gene": round(1000 * per_gene_cell, 2),
        "percell_mixed_projected_full_s": round(
            per_gene_cell * len(logcpm.columns), 1),
        "percell_mixed_pct_converged": round(
            100 * mx_cell["converged"].mean(), 2),
        "elapsed_s": round(time.time() - t_start, 1),
        "timestamp": pd.Timestamp.now().isoformat(timespec="seconds"),
    }
    pd.DataFrame([row]).to_csv(PARTS / f"{stem}_summary.csv", index=False)
    return row


def combine_results(parts=PARTS, results=RESULTS):
    """Rebuild the two headline CSVs from whatever per-type parts exist."""
    sums = []
    nulls = []
    draws = []
    for ct in ALL_TYPES:
        stem = slug(ct)
        sp = Path(parts) / f"{stem}_summary.csv"
        if sp.exists():
            sums.append(pd.read_csv(sp))
        npth = Path(parts) / f"{stem}_null.csv"
        if npth.exists():
            nulls.append(pd.read_csv(npth))
        dpth = Path(parts) / f"{stem}_null_draws.csv"
        if dpth.exists():
            draws.append(pd.read_csv(dpth))
    if not sums:
        raise SystemExit(f"no per-type summaries in {parts}")
    table = pd.concat(sums, ignore_index=True)
    table["_order"] = table["cell_type"].map({c: i for i, c in enumerate(ALL_TYPES)})
    table = table.sort_values("_order").drop(columns="_order")
    table.to_csv(Path(results) / "mixed_vs_others_summary.csv", index=False)
    if nulls:
        pd.concat(nulls, ignore_index=True).to_csv(
            Path(results) / "mixed_null_donor_swap.csv", index=False)
    if draws:
        pd.concat(draws, ignore_index=True).to_csv(
            Path(results) / "mixed_null_draws.csv", index=False)
    return table


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--types", default=None,
                    help="comma-separated subset of ALL_TYPES")
    ap.add_argument("--null-perms", type=int, default=NULL_PERMS,
                    help="max DISTINCT balanced draws (plan space is smaller "
                         "than 200 for every real panel, so this is "
                         "exhaustive)")
    ap.add_argument("--unbalanced-perms", type=int, default=UNBALANCED_PERMS)
    ap.add_argument("--null-gene-set", default="small_effect",
                    choices=["small_effect", "all"])
    ap.add_argument("--null-genes", type=int, default=NULL_GENES)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--combine", action="store_true",
                    help="only rebuild the combined tables from the parts")
    args = ap.parse_args(argv)

    if args.combine:
        table = combine_results()
        print(table.to_string(index=False))
        return

    types = ALL_TYPES
    if args.types:
        types = [t.strip() for t in args.types.split(",") if t.strip()]
        unknown = [t for t in types if t not in ALL_TYPES]
        if unknown:
            raise SystemExit(
                f"unknown cell type(s) {unknown}; choose from {ALL_TYPES}")

    adata = data.load_kang()
    t_all = time.time()
    rows = [run_type(adata, ct, null_perms=args.null_perms,
                     null_genes=args.null_genes, seed=args.seed,
                     unbalanced_perms=args.unbalanced_perms,
                     null_gene_set=args.null_gene_set)
            for ct in types]
    table = pd.DataFrame(rows)
    print("\n--- per-type table ---")
    cols = ["cell_type", "n_genes_tested", "percell_sig", "pseudobulk_sig",
            "mixed_sig", "mixed_pct_converged", "lfc_pearson_pb_mx",
            "null_fp_pb_balanced", "null_fp_mx_balanced"]
    print(table[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    if len(types) == len(ALL_TYPES):
        combine_results()
        print(f"\nwrote combined tables to {RESULTS}")
    print(f"total {time.time() - t_all:.1f}s")


if __name__ == "__main__":
    main()