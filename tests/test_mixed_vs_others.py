"""Tests for the multi-cell-type mixed-model driver.

The driver's job is to make a claim about *every* testable cell type, and the
way that goes wrong is quiet: a cell type silently dropped from the loop, a
null that is not actually null because the draws were unbalanced, or a
per-type row written without the numbers the README quotes. These tests pin
the mechanical properties that make the per-type table trustworthy, on small
synthetic inputs where the right answer is known.

The statistical claims themselves are measured by running the script; what is
tested here is that the machinery around them is not lying.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import mixed_vs_others as mvo  # noqa: E402
from scpipe import permutation  # noqa: E402


def _profiles(n_donors=8, n_genes=30, effect=0.0, seed=0, noise=0.1):
    """Pseudobulk-like frame: one ctrl and one stim profile per donor."""
    rng = np.random.default_rng(seed)
    donors = [f"p{i}" for i in range(n_donors)]
    rows, index = [], []
    for d in donors:
        base = rng.normal(0, 1.0, n_genes)
        for cond in ("ctrl", "stim"):
            shift = effect if cond == "stim" else 0.0
            rows.append(base + shift + rng.normal(0, noise, n_genes))
            index.append(f"{d}||{cond}")
    genes = [f"g{i:02d}" for i in range(n_genes)]
    logcpm = pd.DataFrame(rows, index=index, columns=genes, dtype=float)
    meta = pd.DataFrame({
        "donor": [i.split("||")[0] for i in index],
        "condition": [i.split("||")[1] for i in index],
    }, index=index)
    return logcpm, meta


def test_all_seven_testable_types_are_in_the_loop():
    """Six untested types is the exact error this script exists to fix, so
    the roster is asserted, not assumed."""
    assert len(mvo.ALL_TYPES) == 7
    assert "Megakaryocytes" not in mvo.ALL_TYPES
    for ct in ("CD4 T cells", "CD14+ Monocytes", "B cells", "NK cells",
               "CD8 T cells", "FCGR3A+ Monocytes", "Dendritic cells"):
        assert ct in mvo.ALL_TYPES


def test_bh_reference_is_the_corrected_one():
    """0.05 per gene set, NOT alpha*500 = 25. That confusion is what the
    README correction was about, so the constant is pinned."""
    assert mvo.BH_EXPECTED_FP == pytest.approx(0.05)
    assert mvo.BH_EXPECTED_FP != pytest.approx(0.05 * 500)


def test_balanced_plan_flips_exactly_half():
    units = [f"p{i}" for i in range(8)]
    for seed in range(20):
        plan = permutation.balanced_swap_plan(
            units, np.random.default_rng(seed))
        assert set(plan) == set(units)
        assert sum(plan.values()) == 4, plan


def test_balanced_plan_has_no_binomial_spread():
    """The coin version can flip 0 or 8 donors; the balanced one cannot.
    This is the whole reason the null changed."""
    units = [f"p{i}" for i in range(8)]
    coin_counts = {sum(permutation.sample_swap_plan(
        units, np.random.default_rng(s), 0.5).values()) for s in range(60)}
    assert min(coin_counts) <= 1 or max(coin_counts) >= 7
    bal_counts = {sum(permutation.balanced_swap_plan(
        units, np.random.default_rng(s)).values()) for s in range(60)}
    assert bal_counts == {4}


def test_balanced_swap_holds_out_a_donor_when_the_count_is_odd():
    """Five paired donors cannot be halved: rounding leaves one donor with
    its real effect on every draw, so one is dropped instead."""
    units = [f"p{i}" for i in range(5)]
    plan, used, held = mvo.balanced_swap(units, np.random.default_rng(0))
    assert len(used) == 4 and held in units and held not in used
    assert sum(plan.values()) == 2


def test_balanced_swap_refuses_a_panel_too_small_to_balance():
    with pytest.raises(ValueError):
        mvo.balanced_swap(["p0", "p1"], np.random.default_rng(0))


def test_balanced_null_kills_a_constant_effect():
    """The property the balanced null is for: a real, constant condition
    effect must produce no calls, because half the donors are inverted."""
    logcpm, meta = _profiles(effect=3.0, noise=0.05)
    genes = list(logcpm.columns)
    fp_pb, fp_mx = [], []
    for seed in range(6):
        a, b, _u = mvo.null_draw(logcpm, meta, genes,
                                 np.random.default_rng(seed), balanced=True)
        fp_pb.append(a)
        fp_mx.append(b)
    assert max(fp_pb) <= 1, fp_pb
    assert np.mean(fp_mx) <= 2, fp_mx


def test_unbalanced_null_keeps_the_effect_and_is_not_a_null():
    """Documents the defect being fixed rather than just claiming it: with
    p=0.5 per donor, some draws leave the effect largely intact, so this arm
    calls genes even though it is presented as a null. If this ever stops
    failing to be a null, the README's contrast between the arms is stale."""
    logcpm, meta = _profiles(effect=3.0, noise=0.05)
    genes = list(logcpm.columns)
    counts = [mvo.null_draw(logcpm, meta, genes, np.random.default_rng(s),
                            balanced=False)[0] for s in range(12)]
    assert np.mean(counts) > 0.5, counts


def test_null_draw_rejects_a_broken_permutation(monkeypatch):
    """The donor-structure guard must raise, not annotate. A silent warning
    would let an invalid null ship with a valid-looking number."""
    logcpm, meta = _profiles()
    genes = list(logcpm.columns)

    def broken(samples_before, samples_after, *a, **k):
        return False, ["planted violation"]

    monkeypatch.setattr(permutation, "donor_structure_preserved", broken)
    with pytest.raises(RuntimeError, match="planted violation"):
        mvo.null_draw(logcpm, meta, genes, np.random.default_rng(0),
                      balanced=True)


def test_paired_units_excludes_donors_missing_an_aliquot():
    index = ["p0||ctrl", "p0||stim", "p1||ctrl", "p1||stim", "p2||ctrl"]
    meta = pd.DataFrame({
        "donor": [i.split("||")[0] for i in index],
        "condition": [i.split("||")[1] for i in index],
    }, index=index)
    assert mvo.paired_units(meta) == ["p0", "p1"]


def test_compare_frame_reports_identity_for_identical_arms():
    rng = np.random.default_rng(0)
    df = pd.DataFrame({
        "log2fc": rng.normal(size=50),
        "qval": rng.random(50) * 0.5,
    }, index=[f"g{i}" for i in range(50)])
    cmp = mvo.compare_frame(df, df.copy(), df.copy())
    assert cmp["lfc_pearson_pb_mx"] == pytest.approx(1.0)
    assert cmp["jaccard_pb_mx"] == pytest.approx(1.0)
    assert cmp["percell_sig"] == cmp["pseudobulk_sig"] == cmp["mixed_sig"]


def test_slug_is_stable_and_file_safe():
    assert mvo.slug("CD4 T cells") == "CD4_T_cells"
    assert mvo.slug("CD14+ Monocytes") == "CD14_Monocytes"
    assert mvo.slug("FCGR3A+ Monocytes") == "FCGR3A_Monocytes"
    assert "/" not in mvo.slug("Dendritic cells")


def test_combine_results_orders_types_and_requires_parts(tmp_path):
    parts = tmp_path / "parts"
    parts.mkdir()
    for ct in reversed(mvo.ALL_TYPES):      # write out of order on purpose
        pd.DataFrame([{"cell_type": ct, "n_genes_tested": 10}] * 1).to_csv(
            parts / f"{mvo.slug(ct)}_summary.csv", index=False)
    table = mvo.combine_results(parts=parts, results=tmp_path)
    assert list(table["cell_type"]) == mvo.ALL_TYPES
    assert (tmp_path / "mixed_vs_others_summary.csv").exists()

    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(SystemExit):
        mvo.combine_results(parts=empty, results=tmp_path)
