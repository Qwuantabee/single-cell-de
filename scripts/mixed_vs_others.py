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
donors, so the flipped and unflipped halves cancel and the truth really is
"nothing". The unbalanced version is still measured and reported next to it,
because the difference between the two is the whole reason the balance
matters -- but it is not a null and its numbers are not the headline.

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
# false-discovery proportion among rejections, so on all-null data a
# calibrated procedure rejects a mean of about q * 1 / (1 + ...) ~ 0.05.
# Measured by simulation in scripts/false_positive_test.py: 0.055.
BH_EXPECTED_FP = 0.05

# per-cell fits are the expensive arm: cap them deliberately
PERCELL_GENES = 400
PERCELL_MAX_CELLS = 3200
NULL_GENES = 500
NULL_PERMS = 10
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


def null_draw(logcpm, meta, genes, rng, balanced=True, units=None):
    """One donor-swap null draw; returns (fp_pseudobulk, fp_mixed).

    The permutation moves labels only: `scpipe.permutation` re-checks the donor
    structure on every draw, and a violation raises rather than being printed
    and ignored, because a broken permutation silently invalidates the number
    it produces.
    """
    donors = meta["donor"].astype(str).to_numpy()
    cond = meta["condition"].astype(str).to_numpy()
    units = paired_units(meta) if units is None else list(units)
    if balanced:
        plan, used, _held = balanced_swap(units, rng)
    else:
        # unbalanced: independent coin on the paired donors only, so the two
        # arms differ in the balance of the flips and nothing else
        plan = permutation.sample_swap_plan(units, rng, 0.5)
        used = list(units)

    flipped, _ = permutation.swap_samples(donors, cond, rng, plan=plan)
    ok, problems = permutation.donor_structure_preserved(
        donors, donors, cond, flipped)
    if not ok:
        raise RuntimeError(f"donor structure broken: {problems}")

    # ---- pseudobulk: paired t on the swapped labels ---------------------
    A, B = paired_columns(logcpm, meta, flipped, used, genes)
    _, p_t = stats.ttest_rel(A, B, axis=0)
    fp_pb = int((de._bh(p_t) < ALPHA).sum())

    # ---- mixed: same swapped labels, through the module -----------------
    rows = np.flatnonzero(np.isin(donors, used))
    permuted = logcpm.iloc[rows][genes].copy()
    permuted.index = pd.Index(
        [f"{d}||{c}" for d, c in zip(donors[rows], flipped[rows])])
    pm = mixed.mixed_de(permuted, "donor", "condition", "stim", "ctrl")
    fp_mx = int((de._bh(pm["pval"].to_numpy()) < ALPHA).sum())
    return fp_pb, fp_mx, list(used)


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
             seed=SEED, verbose=True):
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
    log(f"pseudobulk paired t: {de.n_significant(pb):,} sig, {t_pb:.1f}s")

    # ---- arm 3: mixed model per gene -----------------------------------
    t0 = time.time()
    mx = mixed.mixed_de(logcpm, "donor", "condition", "stim", "ctrl")
    t_mx = time.time() - t0
    n_conv = int(mx["converged"].sum())
    log(f"mixedlm per gene: {de.n_significant(mx):,} sig, {t_mx:.1f}s, "
        f"converged {n_conv:,}/{len(mx):,} "
        f"({100 * n_conv / len(mx):.1f}%)")

    cmp = compare_frame(pb, mx, per_cell)

    isg = [g for g in data.ISG_MARKERS if g in logcpm.columns]
    isg_sig = {name: sum(frame.loc[g, "qval"] < ALPHA for g in isg)
               for name, frame in (("percell", per_cell),
                                   ("pseudobulk", pb), ("mixed", mx))}
    log(f"ISG recovery of {len(isg)}: {isg_sig}")

    # ---- null: balanced donor swap (the null) + unbalanced (for contrast)
    rng = np.random.default_rng(seed)
    null_genes_list = list(logcpm.columns[
        rng.choice(len(logcpm.columns),
                   size=min(null_genes, len(logcpm.columns)), replace=False)])
    fp_bal = {"pseudobulk": [], "mixed": []}
    fp_unb = {"pseudobulk": [], "mixed": []}
    held_outs = []
    t0 = time.time()
    for i in range(null_perms):
        draw_rng = np.random.default_rng(int(rng.integers(1 << 31)))
        a, b, used = null_draw(logcpm, meta, null_genes_list, draw_rng,
                               balanced=True, units=units)
        if i == 0:
            held_outs.append(len(units) - len(used))
        fp_bal["pseudobulk"].append(a)
        fp_bal["mixed"].append(b)
    t_bal = time.time() - t0
    t0 = time.time()
    for _ in range(null_perms):
        draw_rng = np.random.default_rng(int(rng.integers(1 << 31)))
        a, b, _used = null_draw(logcpm, meta, null_genes_list, draw_rng,
                                balanced=False, units=units)
        fp_unb["pseudobulk"].append(a)
        fp_unb["mixed"].append(b)
    t_unb = time.time() - t0
    n_genes_null = len(null_genes_list)
    log(f"null: {null_perms} balanced draws in {t_bal:.1f}s, unbalanced in "
        f"{t_unb:.1f}s")
    log(f"  balanced   pseudobulk {np.mean(fp_bal['pseudobulk']):.2f} "
        f"mixed {np.mean(fp_bal['mixed']):.2f} (of {n_genes_null})")
    log(f"  unbalanced pseudobulk {np.mean(fp_unb['pseudobulk']):.2f} "
        f"mixed {np.mean(fp_unb['mixed']):.2f} (of {n_genes_null})")

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
        "mx_donor_var": mx["donor_var"],
    })
    models.to_csv(PARTS / f"{stem}_models.csv")

    null_rows = []
    for method in ("pseudobulk", "mixed"):
        null_rows.append({
            "cell_type": cell_type, "method": method, "balance": "balanced",
            "mean_false_positives": float(np.mean(fp_bal[method])),
            "perms": null_perms, "genes_tested": n_genes_null,
            "per_draw": ";".join(str(v) for v in fp_bal[method]),
            "expected_fp_bh": BH_EXPECTED_FP,
            "ratio_observed_expected": (
                float(np.mean(fp_bal[method])) / BH_EXPECTED_FP),
        })
        null_rows.append({
            "cell_type": cell_type, "method": method, "balance": "unbalanced",
            "mean_false_positives": float(np.mean(fp_unb[method])),
            "perms": null_perms, "genes_tested": n_genes_null,
            "per_draw": ";".join(str(v) for v in fp_unb[method]),
            "expected_fp_bh": BH_EXPECTED_FP,
            "ratio_observed_expected": (
                float(np.mean(fp_unb[method])) / BH_EXPECTED_FP),
        })
    null_df = pd.DataFrame(null_rows)
    null_df["permutation_preserved_donors"] = True
    null_df.to_csv(PARTS / f"{stem}_null.csv", index=False)

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
        "jaccard_pb_mx": round(cmp["jaccard_pb_mx"], 4),
        "jaccard_percell_mx": round(cmp["jaccard_percell_mx"], 4),
        "lfc_pearson_pb_mx": round(cmp["lfc_pearson_pb_mx"], 4),
        "lfc_spearman_pb_mx": round(cmp["lfc_spearman_pb_mx"], 4),
        "mixed_s": round(t_mx, 1),
        "mixed_ms_per_gene_pseudobulk": round(1000 * t_mx / len(mx), 2),
        "pseudobulk_s": round(t_pb, 1),
        "percell_s": round(t_percell, 1),
        "isg_present": len(isg),
        "isg_percell": isg_sig["percell"],
        "isg_pseudobulk": isg_sig["pseudobulk"],
        "isg_mixed": isg_sig["mixed"],
        "null_fp_pb_balanced": round(float(np.mean(fp_bal["pseudobulk"])), 2),
        "null_fp_mx_balanced": round(float(np.mean(fp_bal["mixed"])), 2),
        "null_ratio_pb_balanced": round(
            float(np.mean(fp_bal["pseudobulk"])) / BH_EXPECTED_FP, 2),
        "null_ratio_mx_balanced": round(
            float(np.mean(fp_bal["mixed"])) / BH_EXPECTED_FP, 2),
        "null_fp_pb_unbalanced": round(
            float(np.mean(fp_unb["pseudobulk"])), 2),
        "null_fp_mx_unbalanced": round(
            float(np.mean(fp_unb["mixed"])), 2),
        "null_genes": n_genes_null,
        "null_perms": null_perms,
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
    for ct in ALL_TYPES:
        stem = slug(ct)
        sp = Path(parts) / f"{stem}_summary.csv"
        if sp.exists():
            sums.append(pd.read_csv(sp))
        npth = Path(parts) / f"{stem}_null.csv"
        if npth.exists():
            nulls.append(pd.read_csv(npth))
    if not sums:
        raise SystemExit(f"no per-type summaries in {parts}")
    table = pd.concat(sums, ignore_index=True)
    table["_order"] = table["cell_type"].map({c: i for i, c in enumerate(ALL_TYPES)})
    table = table.sort_values("_order").drop(columns="_order")
    table.to_csv(Path(results) / "mixed_vs_others_summary.csv", index=False)
    if nulls:
        pd.concat(nulls, ignore_index=True).to_csv(
            Path(results) / "mixed_null_donor_swap.csv", index=False)
    return table


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--types", default=None,
                    help="comma-separated subset of ALL_TYPES")
    ap.add_argument("--null-perms", type=int, default=NULL_PERMS)
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
                     null_genes=args.null_genes, seed=args.seed)
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
