import numpy as np
import pandas as pd
import pytest

from scpipe import pseudobulk, simulate


@pytest.fixture(scope="module")
def small():
    adata, truth = simulate.simulate(n_genes=60, n_donors=4,
                                     cells_per_donor=50, n_de=5, seed=0)
    return adata, truth


def test_aggregate_shape(small):
    adata, _ = small
    mat, meta = pseudobulk.aggregate(adata, "donor", "condition")
    assert mat.shape == (8, 60)          # 4 donors x 2 conditions
    assert list(meta.columns) == ["donor", "condition", "n_cells"]


def test_aggregate_sums_counts(small):
    """The pseudobulk profile must be the literal sum of its cells."""
    adata, _ = small
    mat, meta = pseudobulk.aggregate(adata, "donor", "condition")
    key = mat.index[0]
    donor, cond = key.split("||")
    mask = ((adata.obs["donor"] == donor) &
            (adata.obs["condition"] == cond)).to_numpy()
    expected = np.asarray(adata.X[mask].sum(axis=0)).ravel()
    assert np.allclose(mat.loc[key].to_numpy(), expected)


def test_aggregate_drops_tiny_groups(small):
    """A group with too few cells is dropped, the rest are kept.

    Pseudobulk from a handful of cells is noise dressed as a replicate, and
    leaving it in widens the variance estimate for every gene.
    """
    adata, _ = small
    # keep all 50 cells of each group except one, cut down to 4
    donor0 = adata.obs["donor"] == "donor_0"
    ctrl = adata.obs["condition"] == "ctrl"
    victim = np.flatnonzero((donor0 & ctrl).to_numpy())
    keep = np.ones(adata.n_obs, dtype=bool)
    keep[victim[4:]] = False
    trimmed = adata[keep].copy()

    mat, meta = pseudobulk.aggregate(trimmed, "donor", "condition",
                                     min_cells=10)
    assert mat.shape[0] == 7                      # 8 groups, 1 dropped
    assert "donor_0||ctrl" not in mat.index
    assert (meta["n_cells"] >= 10).all()


def test_aggregate_raises_when_everything_filtered(small):
    adata, _ = small
    with pytest.raises(ValueError, match="no .* combination"):
        pseudobulk.aggregate(adata, "donor", "condition", min_cells=10_000)


def test_size_factors_centred_and_positive(small):
    adata, _ = small
    mat, _ = pseudobulk.aggregate(adata, "donor", "condition")
    sf = pseudobulk.size_factors(mat)
    assert (sf > 0).all()
    # median of ratios is centred on 1 by construction
    assert 0.5 < np.median(sf) < 2.0


def test_size_factors_track_sequencing_depth():
    """Doubling a sample's counts must double its size factor, so that
    normalisation undoes it exactly."""
    rng = np.random.default_rng(0)
    base = rng.poisson(50, size=(4, 100)).astype(float)
    mat = pd.DataFrame(base, index=[f"s{i}" for i in range(4)],
                       columns=[f"g{j}" for j in range(100)])
    mat.iloc[0] *= 2
    sf = pseudobulk.size_factors(mat)
    assert sf.iloc[0] / sf.iloc[1:].mean() == pytest.approx(2.0, rel=0.1)


def test_median_ratio_resists_composition_bias():
    """The property that total-count normalisation lacks.

    A few genes go up 10x. Total counts rise, so CPM makes every other gene
    look like it went down. Median-of-ratios should not move.
    """
    rng = np.random.default_rng(1)
    counts = rng.poisson(100, size=(6, 200)).astype(float)
    mat = pd.DataFrame(counts, index=[f"s{i}" for i in range(6)],
                       columns=[f"g{j}" for j in range(200)])
    treated = mat.index[3:]
    mat.loc[treated, mat.columns[:10]] *= 10      # 5% of genes, 10x up

    null_genes = mat.columns[10:]
    cpm_vals = pseudobulk.cpm(mat, log=True)
    med_vals = pseudobulk.normalise(mat, log=True)

    cpm_shift = (cpm_vals.loc[treated, null_genes].mean().mean()
                 - cpm_vals.loc[mat.index[:3], null_genes].mean().mean())
    med_shift = (med_vals.loc[treated, null_genes].mean().mean()
                 - med_vals.loc[mat.index[:3], null_genes].mean().mean())

    assert abs(cpm_shift) > 0.3                 # CPM is pushed down
    assert abs(med_shift) < 0.05                # median-of-ratios is not
    assert abs(med_shift) < abs(cpm_shift)


def test_filter_genes_removes_undetected():
    mat = pd.DataFrame({"on": [100, 120, 90, 110],
                        "off": [0, 0, 1, 0],
                        "sparse": [500, 0, 0, 0]})
    kept = pseudobulk.filter_genes(mat, min_count=10, min_samples=3)
    assert list(kept.columns) == ["on"]


def test_cpm_rows_sum_to_a_million():
    mat = pd.DataFrame(np.array([[1.0, 3.0], [10.0, 10.0]]))
    out = pseudobulk.cpm(mat, log=False)
    assert np.allclose(out.sum(axis=1).to_numpy(), 1e6)
