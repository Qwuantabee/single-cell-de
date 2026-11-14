"""Tests for cell-type-stratified pseudobulk.

The behaviour that matters here is a specific failure mode: if the cell
composition differs between conditions, a pooled all-cells test reports the
composition difference as expression change. These tests build that
situation on purpose and check that the stratified version does not fall for
it, because a test that only checks shapes would pass on a broken
implementation.
"""

import numpy as np
import pandas as pd
import pytest

from scpipe import pseudobulk


def _mixture(n_donors=6, cells=120, seed=0):
    """Two cell types, one of which dominates one condition.

    Gene A is a marker of type_x and gene B of type_y; both are constant per
    cell within their type, in both conditions. So there is zero expression
    change. The only difference between conditions is the number of cells of
    each type recovered -- which is what a pooled test will mistake for a
    change in both genes.
    """
    import anndata as ad
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    rows, obs = [], []
    for d in range(n_donors):
        for cond in ("ctrl", "stim"):
            # ctrl is mostly type_y, stim is mostly type_x
            n_x = int(cells * (0.15 if cond == "ctrl" else 0.85))
            n_y = cells - n_x
            for ct, n in (("type_x", n_x), ("type_y", n_y)):
                block = np.zeros((n, 4), dtype=np.float32)
                block[:, 0] = 200 if ct == "type_x" else 5
                block[:, 1] = 5 if ct == "type_x" else 200
                block[:, 2] = 100
                block[:, 3] = 50
                block += rng.poisson(2, size=block.shape)
                rows.append(block)
                obs += [(f"donor_{d}", cond, ct)] * n

    X = sp.csr_matrix(np.vstack(rows))
    o = pd.DataFrame(obs, columns=["donor", "condition", "cell_type"])
    o.index = [f"c{i}" for i in range(len(o))]
    var = pd.DataFrame(index=["marker_x", "marker_y", "shared", "shared2"])
    return ad.AnnData(X=X, obs=o, var=var)


def test_composition_reports_the_imbalance():
    adata = _mixture()
    props, shift = pseudobulk.cell_type_composition(adata)
    assert np.allclose(props.sum(axis=1).to_numpy(), 1.0)
    # type_x is 6x more common in whichever condition is listed first
    assert abs(shift.loc["type_x", "log2_ratio"]) > 2.0
    assert abs(shift.loc["type_y", "log2_ratio"]) > 2.0
    # the two types move in opposite directions
    assert (shift.loc["type_x", "log2_ratio"]
            * shift.loc["type_y", "log2_ratio"]) < 0
    assert shift.loc["type_x", "n_a"] > shift.loc["type_x", "n_b"] or \
        shift.loc["type_x", "n_a"] < shift.loc["type_x", "n_b"]


def test_composition_rejects_single_condition():
    adata = _mixture()
    adata = adata[adata.obs["condition"] == "ctrl"].copy()
    with pytest.raises(ValueError, match="expected 2 conditions"):
        pseudobulk.cell_type_composition(adata)


def test_stratified_test_sees_no_expression_change():
    """The whole point: the pooled shift is composition, not biology.

    Pooling calls both markers; stratifying by cell type does not, because
    within a cell type nothing changed.
    """
    adata = _mixture()

    # pooled: aggregate all cells of a donor together
    pooled, _ = pseudobulk.aggregate(adata, "donor", "condition", min_cells=10)
    pooled = pooled.loc[:, ["marker_x", "marker_y"]]
    a = pooled.loc[pooled.index.str.endswith("ctrl")]
    b = pooled.loc[pooled.index.str.endswith("stim")]
    pooled_shift = (b.mean() - a.mean()).abs()
    assert (pooled_shift > 1.0).all(), "pooled test should be fooled"

    results, composition, skipped = pseudobulk.pseudobulk_by_celltype(
        adata, group_a="ctrl", group_b="stim", min_cells=10, min_donors=4)
    assert not skipped
    assert set(results.index.get_level_values("cell_type")) == {
        "type_x", "type_y"}
    # within a cell type the contrast is exact, so the fold change is ~0
    for gene in ("marker_x", "marker_y"):
        sub = results.xs(gene, level="gene")
        assert (sub["log2fc"].abs() < 1.0).all(), (
            f"{gene} still moves within cell type: {sub['log2fc'].to_dict()}")


def test_stratified_index_and_columns():
    adata = _mixture()
    results, composition, skipped = pseudobulk.pseudobulk_by_celltype(
        adata, group_a="ctrl", group_b="stim", min_cells=10, min_donors=4)
    assert results.index.names == ["cell_type", "gene"]
    for col in ("log2fc", "pval", "qval", "qval_global", "n_donors"):
        assert col in results.columns
    assert (results["n_donors"] == 6).all()
    # global BH can never be smaller than the within-type correction
    assert (results["qval_global"] >= results["qval"] - 1e-12).all()


def test_skips_cell_types_with_too_few_donors():
    adata = _mixture(n_donors=6)
    # donor_0 contributes no type_x cells in ctrl, so type_x has 5 usable
    # donors and must be skipped rather than tested on a subset of donors
    # while type_y is tested on all 6. Silently dropping donors would make
    # the two cell types incomparable.
    drop = ((adata.obs["donor"] == "donor_0")
            & (adata.obs["cell_type"] == "type_x")
            & (adata.obs["condition"] == "ctrl"))
    assert drop.sum() > 0
    trimmed = adata[~drop.to_numpy()].copy()
    results, _, skipped = pseudobulk.pseudobulk_by_celltype(
        trimmed, group_a="ctrl", group_b="stim", min_cells=10, min_donors=6)
    assert "type_y" in results.index.get_level_values("cell_type")
    assert "type_x" in skipped
    assert "donors" in skipped["type_x"]


def test_raises_when_nothing_is_testable():
    adata = _mixture(n_donors=2)
    with pytest.raises(ValueError, match="no cell type had enough donors"):
        pseudobulk.pseudobulk_by_celltype(
            adata, group_a="ctrl", group_b="stim", min_cells=10,
            min_donors=4)


def test_significant_by_celltype_counts_match_filter():
    adata = _mixture()
    results, _, _ = pseudobulk.pseudobulk_by_celltype(
        adata, group_a="ctrl", group_b="stim", min_cells=10, min_donors=4)
    counts = pseudobulk.significant_by_celltype(results, alpha=0.05)
    for ct in counts.index:
        expected = int((results.loc[ct, "qval"] < 0.05).sum())
        assert counts.loc[ct] == expected
