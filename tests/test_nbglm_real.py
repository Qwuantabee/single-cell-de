"""Tests for the real-data NB GLM comparison script.

The script's job is to report numbers, so most of it is I/O. What can be got
wrong quietly, and is therefore tested here, is the bookkeeping: a call set
that silently drops NaN q-values, an overlap that double-counts, an agreement
statistic computed on genes the methods did not both test, an ISG rank that
ignores q = NaN, and -- the one that matters most -- a permutation that moves
cells between donors or fails to swap a donor's labels exactly.

The permutation tests are written against the same failure the module itself
guards: a label swap that is really a cell shuffle still looks like a null.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "nbglm_real.py"


def _load():
    spec = importlib.util.spec_from_file_location("nbglm_real", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


nr = _load()


def _frame(rows):
    """Result frame from {gene: (log2fc, qval)}, the columns the helpers use."""
    idx = list(rows)
    return pd.DataFrame({
        "log2fc": [rows[g][0] for g in idx],
        "qval": [rows[g][1] for g in idx],
    }, index=idx)


def _adata(labels):
    """Minimal object with the two obs columns the swap helpers read."""
    import anndata as ad

    obs = pd.DataFrame({
        "donor": [d for d, _ in labels],
        "condition": [c for _, c in labels],
    })
    return ad.AnnData(X=np.zeros((len(labels), 2)), obs=obs)


# ------------------------------------------------------------- call sets


def test_called_genes_ignores_nan_and_ties_alpha_strictly():
    r = _frame({"a": (1.0, 0.01), "b": (1.0, 0.05), "c": (1.0, 0.049),
                "d": (1.0, np.nan)})
    # q == alpha is not a call: the threshold is strict
    assert nr.called_genes(r) == {"a", "c"}


def test_called_genes_on_empty_frame():
    assert nr.called_genes(pd.DataFrame(columns=["log2fc", "qval"])) == set()
    assert nr.called_genes(None) == set()


def test_overlap_and_jaccard():
    pb = _frame({"a": (1, 0.01), "b": (1, 0.01), "c": (1, 0.9)})
    nb = _frame({"b": (1, 0.01), "c": (1, 0.01), "d": (1, 0.01)})
    st = nr.overlap_stats(pb, nb)
    assert (st["n_pb"], st["n_nb"]) == (2, 3)
    assert st["n_both"] == 1
    assert st["n_pb_only"] == 1 and st["n_nb_only"] == 2
    assert st["jaccard"] == pytest.approx(1 / 4)   # union a,b,c,d
    assert st["pb_only_genes"] == ["a"]
    assert st["nb_only_genes"] == ["c", "d"]


def test_jaccard_is_nan_when_nobody_calls_anything():
    empty = _frame({"a": (1, 0.9)})
    assert np.isnan(nr.overlap_stats(empty, empty)["jaccard"])


def test_identical_call_sets_give_jaccard_one():
    same = _frame({"a": (1, 0.01), "b": (1, 0.02)})
    assert nr.overlap_stats(same, same)["jaccard"] == pytest.approx(1.0)


# -------------------------------------------------------------- agreement


def test_log2fc_agreement_uses_only_shared_genes():
    pb = _frame({"a": (1.0, 0.01), "b": (2.0, 0.01)})
    nb = _frame({"b": (2.0, 0.01), "z": (99.0, 0.5)})
    agr = nr.log2fc_agreement(pb, nb)
    assert agr["n_genes"] == 1          # z is not shared, so it cannot dilute
    assert np.isnan(agr["pearson"])     # one point has no correlation


def test_log2fc_agreement_perfect_and_reversed():
    pb = _frame({"a": (1.0, 0.1), "b": (2.0, 0.1), "c": (3.0, 0.1),
                 "d": (4.0, 0.1)})
    same = nr.log2fc_agreement(pb, pb)
    assert same["pearson"] == pytest.approx(1.0)
    assert same["spearman"] == pytest.approx(1.0)
    assert same["mean_abs_diff"] == pytest.approx(0.0)

    flipped = _frame({"a": (-1.0, 0.1), "b": (-2.0, 0.1), "c": (-3.0, 0.1),
                      "d": (-4.0, 0.1)})
    rev = nr.log2fc_agreement(pb, flipped)
    assert rev["pearson"] == pytest.approx(-1.0)
    assert rev["mean_abs_diff"] == pytest.approx(5.0)


def test_log2fc_agreement_drops_nan_rows():
    pb = _frame({"a": (1.0, 0.1), "b": (2.0, 0.1), "c": (3.0, 0.1),
                 "d": (4.0, 0.1)})
    nb = _frame({"a": (1.0, 0.1), "b": (2.0, 0.1), "c": (3.0, 0.1),
                 "d": (np.nan, 0.1)})
    assert nr.log2fc_agreement(pb, nb)["n_genes"] == 3


# --------------------------------------------------------------- ISG rank


def test_rank_by_q_puts_nan_last_and_shares_no_rank():
    r = _frame({"a": (1, 0.05), "b": (1, 0.001), "c": (1, np.nan)})
    ranks = nr.rank_by_q(r)
    assert ranks["b"] == 1 and ranks["a"] == 2 and ranks["c"] == 3


def test_isg_table_marks_unmeasured_marker_absent():
    pb = _frame({"ISG15": (2.0, 0.001), "MX1": (1.0, 0.9)})
    nb = _frame({"ISG15": (2.0, 0.001), "MX1": (1.0, 0.001)})
    tab = nr.isg_table(pb, nb, markers=["ISG15", "MX1", "IFIT1"])
    assert tab.loc["ISG15", "called_pb"] and tab.loc["ISG15", "called_nb"]
    assert not tab.loc["MX1", "called_pb"] and tab.loc["MX1", "called_nb"]
    assert tab.loc["IFIT1", "rank_pb"] == -1      # not in the frame at all
    assert not tab.loc["IFIT1", "called_pb"]


def test_isg_median_rank_ignores_markers_absent_from_one_method():
    pb = _frame({"ISG15": (2, 0.001), "MX1": (2, 0.02), "OAS1": (2, 0.5)})
    nb = _frame({"ISG15": (2, 0.001), "MX1": (2, 0.5), "OAS1": (2, 0.001)})
    tab = nr.isg_table(pb, nb, markers=["ISG15", "MX1", "OAS1", "IFIT3"])
    med = nr.isg_median_rank(tab)
    assert med["n_isg"] == 3                    # IFIT3 was never measured
    assert med["median_rank_pb"] == 2.0         # ranks 1, 2, 3 -> median 2
    assert med["median_rank_nb"] == 2.0
    assert med["n_called_pb"] == 2 and med["n_called_nb"] == 2


# ------------------------------------------------------------ null gene set


def test_null_gene_set_excludes_markers_and_large_effects():
    pb = _frame({"small": (0.1, 0.5), "big": (2.0, 0.001),
                 "ISG15": (0.1, 0.5), "other": (0.4, 0.5)})
    nb = _frame({"small": (0.1, 0.5), "big": (2.0, 0.001),
                 "ISG15": (0.1, 0.5), "other": (0.4, 0.5),
                 "only_nb": (0.0, 0.5)})
    got = nr.null_gene_set(pb, nb)
    assert "small" in got and "other" in got
    assert "big" not in got                     # real effect, not null
    assert "ISG15" not in got                   # marker, excluded by name
    assert "only_nb" not in got                 # not tested by both


# ---------------------------------------------------------- swap invariants


def test_check_cell_swap_accepts_a_real_exchange():
    labels = ([("d1", "ctrl")] * 3 + [("d1", "stim")] * 2
              + [("d2", "ctrl")] * 4 + [("d2", "stim")] * 1)
    ad = _adata(labels)
    a = ad.copy()
    cond = ad.obs["condition"].tolist()
    a.obs["condition"] = ["stim"] * 3 + ["ctrl"] * 2 + cond[5:]
    ok, problems = nr.check_cell_swap(
        ad.obs["donor"], a.obs["donor"], ad.obs["condition"],
        a.obs["condition"])
    assert ok, problems


def test_check_cell_swap_rejects_a_cell_moving_donor():
    ad = _adata([("d1", "ctrl"), ("d1", "stim"), ("d2", "ctrl"),
                 ("d2", "stim")])
    moved = ad.copy()
    moved.obs["donor"] = ["d2", "d1", "d1", "d2"]
    ok, problems = nr.check_cell_swap(
        ad.obs["donor"], moved.obs["donor"], ad.obs["condition"],
        ad.obs["condition"])
    assert not ok and "changed donor" in problems[0]


def test_check_cell_swap_rejects_a_partial_flip():
    """Only some of a donor's cells relabelled is not a swap.

    The donor here is deliberately unequal (3 ctrl, 1 stim): with equal
    counts an exchange and a partial relabelling produce the same numbers, so
    an equal-count donor cannot detect this failure at all.
    """
    ad = _adata([("d1", "ctrl"), ("d1", "ctrl"), ("d1", "ctrl"),
                 ("d1", "stim"), ("d2", "ctrl"), ("d2", "stim")])
    bad = ad.copy()
    bad.obs["condition"] = ["stim", "ctrl", "ctrl", "stim",
                            "stim", "ctrl"]
    ok, problems = nr.check_cell_swap(
        ad.obs["donor"], bad.obs["donor"], ad.obs["condition"],
        bad.obs["condition"])
    assert not ok and "is not an exchange" in problems[0]


def test_n_donors_flipped_counts_only_changed_donors():
    ad = _adata([("d1", "ctrl"), ("d1", "stim"), ("d2", "ctrl"),
                 ("d2", "stim"), ("d3", "ctrl"), ("d3", "stim")])
    after = ad.copy()
    after.obs["condition"] = ["stim", "ctrl", "ctrl", "stim", "stim", "ctrl"]
    assert nr._n_donors_flipped(ad, after) == 2


def test_swap_preserves_donor_column_and_nothing_else():
    labels = []
    for d in ("d1", "d2", "d3", "d4"):
        labels += [(d, "ctrl")] * 7 + [(d, "stim")] * 5
    ad = _adata(labels)
    rng = np.random.default_rng(0)
    perm = nr.permute_adata(ad, rng)
    assert np.array_equal(ad.obs["donor"].to_numpy(),
                          perm.obs["donor"].to_numpy())
    ok, problems = nr.check_cell_swap(
        ad.obs["donor"], perm.obs["donor"], ad.obs["condition"],
        perm.obs["condition"])
    assert ok, problems


def test_balanced_plan_flips_exactly_half_the_donors():
    labels = []
    for d in [f"d{i}" for i in range(8)]:
        labels += [(d, "ctrl")] * 4 + [(d, "stim")] * 4
    ad = _adata(labels)
    for seed in range(5):
        rng = np.random.default_rng(seed)
        perm = nr.permute_adata(ad, np.random.default_rng(seed + 100),
                                balanced=True)
        assert nr._n_donors_flipped(ad, perm) == 4
        ok, problems = nr.check_cell_swap(
            ad.obs["donor"], perm.obs["donor"], ad.obs["condition"],
            perm.obs["condition"])
        assert ok, problems
        del rng


def test_fair_coin_plan_is_not_forced_to_half():
    """The native draw must stay a coin flip: forcing k = n/2 would hide it."""
    labels = []
    for d in [f"d{i}" for i in range(8)]:
        labels += [(d, "ctrl")] * 4 + [(d, "stim")] * 4
    ad = _adata(labels)
    ks = {nr._n_donors_flipped(ad, nr.permute_adata(
        ad, np.random.default_rng(s))) for s in range(20)}
    assert len(ks) > 1


# ------------------------------------------------------------- null summary


def test_null_summary_reports_balanced_separately():
    recs = [
        {"fp_pseudobulk": 0, "fp_nb_glm": 0, "fp_pseudobulk_null": 0,
         "fp_nb_glm_null": 0, "n_pseudobulk": 10, "n_nb_glm": 10,
         "balanced": True},
        {"fp_pseudobulk": 0, "fp_nb_glm": 9, "fp_pseudobulk_null": 0,
         "fp_nb_glm_null": 1, "n_pseudobulk": 10, "n_nb_glm": 10,
         "balanced": False},
    ]
    tab = nr.null_summary(recs).set_index("method")
    assert tab.loc["nb_glm", "mean_fp"] == pytest.approx(4.5)
    assert tab.loc["nb_glm", "mean_fp_balanced"] == pytest.approx(0.0)
    assert tab.loc["nb_glm", "max_fp"] == 9
    assert tab.loc["pseudobulk", "mean_fp"] == pytest.approx(0.0)


def test_null_summary_marks_no_balanced_draws():
    recs = [{"fp_pseudobulk": 1, "fp_nb_glm": 2, "fp_pseudobulk_null": 1,
             "fp_nb_glm_null": 2, "n_pseudobulk": 10, "n_nb_glm": 10,
             "balanced": False}]
    tab = nr.null_summary(recs).set_index("method")
    assert np.isnan(tab.loc["nb_glm", "mean_fp_balanced"])
    assert tab.loc["nb_glm", "n_balanced"] == 0
