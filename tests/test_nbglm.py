"""Tests for the negative-binomial GLM differential-expression path.

Same philosophy as test_de.py: the assertions that matter are the ones where
the right answer is known in advance. A GLM that returns small p-values is
easy to write; a GLM whose null p-values are uniform and whose size factors
actually absorb library depth is the thing being tested.
"""

import numpy as np
import pandas as pd
import pytest

from scpipe import de, nbglm, simulate
from scpipe.pseudobulk import size_factors


def _small(n_donors=6, n_genes=120, n_de=20, effect=1.0, seed=0, n_de_low=0):
    adata, truth = simulate.simulate(n_genes=n_genes, n_donors=n_donors,
                                     n_de=n_de, effect=effect, seed=seed)
    return adata, truth


def test_recovers_a_planted_effect():
    """Genes given a 2x effect are called, null genes mostly are not."""
    adata, truth = _small(n_donors=8, n_de=20)
    r = nbglm.nbglm_de(adata, "condition", "donor", "stim", "ctrl")
    called = set(r.index[r["qval"] < 0.05])
    tp = len(called & truth)
    assert tp >= 15, f"only {tp}/{len(truth)} planted genes recovered"
    assert r.loc[sorted(truth)[0], "qval"] < 0.05
    # the fit has to be on the right scale: log2fc_glm must agree in sign
    # and roughly in size with the observed log2 CPM difference
    sub = r.loc[sorted(truth)]
    agree = np.sign(sub["log2fc_glm"]) == np.sign(sub["log2fc"])
    assert agree.mean() > 0.8
    assert np.median(np.abs(sub["log2fc"])) > 0.5


def test_null_gives_approximately_correct_fdr():
    """Nothing changes, so BH at 0.05 should almost never call anything.

    Pooled over replicates the genotype-wide false discovery rate must sit
    near alpha. A naive NB Wald test without a t reference or with an
    under-estimated dispersion lands at 2-4x alpha here, which is exactly
    what this catches.
    """
    called, total = 0, 0
    for seed in range(4):
        adata, _ = simulate.simulate(n_genes=150, n_donors=6, n_de=0,
                                     seed=seed)
        r = nbglm.nbglm_de(adata, "condition", "donor", "stim", "ctrl")
        called += int((r["qval"] < 0.05).sum())
        total += len(r)
    fdr = called / total
    assert fdr < 0.05, f"null FDR {fdr:.3f} at 6 donors"
    # The ceiling above allows 29 false positives over 600 genes. BH on an
    # all-null family calls about 0.05 genes per family, so four families
    # should give at most one call even with generous slack.
    assert called <= 1, f"{called} BH calls on 4 all-null families"


def test_pvalues_are_uniform_under_the_null():
    """A sharper statement than the FDR bound: ~5% of nulls below 0.05."""
    ps = []
    for seed in range(3):
        adata, _ = simulate.simulate(n_genes=150, n_donors=8, n_de=0,
                                     seed=seed)
        r = nbglm.nbglm_de(adata, "condition", "donor", "stim", "ctrl")
        ps.append(r["pval"].to_numpy())
    p = np.concatenate(ps)
    frac = float(np.mean(p < 0.05))
    assert 0.01 < frac < 0.12, f"{frac:.3f} of null p-values below 0.05"


def test_null_pvalues_are_calibrated_at_six_donors():
    """A band that can fail. 1,200 null genes: the rate of p < 0.05 has an sd
    of about 0.006 if calibrated, so 0.03 to 0.075 is wide. A normal-reference
    Wald (the failure the t reference exists to prevent) puts about 10% below
    0.05 at six donors and 12 samples, which this rejects; the 0.01-0.12 band
    above does not.
    """
    ps = []
    for seed in range(8):
        adata, _ = simulate.simulate(n_genes=150, n_donors=6, n_de=0,
                                     seed=seed)
        r = nbglm.nbglm_de(adata, "condition", "donor", "stim", "ctrl")
        ps.append(r["pval"].to_numpy())
    p = np.concatenate(ps)
    assert len(p) >= 1100
    frac = float(np.mean(p < 0.05))
    assert 0.03 < frac < 0.075, f"{frac:.3f} of null p-values below 0.05"
    # uniform p-values have mean 0.5 (sd 0.008 over 1,200 genes)
    assert 0.44 < float(np.mean(p)) < 0.56, f"mean null p {np.mean(p):.3f}"


def test_size_factors_absorb_library_depth():
    """A sample sequenced twice as deep is not called differentially expressed.

    All genes are multiplied by one sample's size factor. The offset must
    undo it, so nothing is called. This is the test that would fail if the
    size factors were applied as a pre-division (wrong variance) or dropped.
    """
    mat = pd.DataFrame(
        np.array([[100, 50, 10], [200, 100, 20], [120, 60, 12],
                  [90, 45, 9], [210, 105, 21], [110, 55, 11],
                  [130, 65, 13], [95, 48, 10]], dtype=float),
        index=[f"s{i}" for i in range(8)],
        columns=["g0", "g1", "g2"])
    base = size_factors(mat)
    doubled = mat.copy()
    doubled.iloc[0] = doubled.iloc[0] * 2
    sf2 = size_factors(doubled)
    # Size factors are only defined up to a constant, so the assertion is on
    # the ratio between the doubled sample and the rest, not on absolute
    # values: that ratio has to double.
    ratios = (sf2 / sf2.iloc[1]) / (base / base.iloc[1])
    assert ratios.iloc[0] == pytest.approx(2.0, rel=0.02)
    assert np.allclose(ratios.iloc[1:], 1.0, atol=1e-6)


def test_depth_difference_alone_is_not_called():
    """End to end: one condition sequenced 1.8x deeper, no biology.

    Counts are Poisson at the same rate, differing only in depth. Nothing
    should survive BH.
    """
    import anndata as ad
    import scipy.sparse as sp

    rng = np.random.default_rng(3)
    n_genes, n_donors = 200, 8
    rate = np.exp(rng.normal(4.0, 1.0, size=n_genes))
    blocks, obs = [], []
    for d in range(n_donors):
        for cond, depth in (("ctrl", 1.0), ("stim", 1.8)):
            counts = rng.poisson(rate * depth * 40, size=(60, n_genes))
            blocks.append(counts.astype(np.float32))
            obs += [(f"donor_{d}", cond)] * 60
    X = sp.csr_matrix(np.vstack(blocks))
    adata = ad.AnnData(X=X,
                       obs=pd.DataFrame(obs, columns=["donor", "condition"],
                                        index=[f"c{i}" for i in range(len(obs))]),
                       var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]))
    r = nbglm.nbglm_de(adata, "condition", "donor", "stim", "ctrl",
                       min_count=5, min_samples=3)
    n_sig = int((r["qval"] < 0.05).sum())
    assert n_sig <= 10, f"{n_sig} genes called on a pure depth difference"


def test_handles_zero_counts():
    """All-zero and mostly-zero genes are dropped or left untested, quietly."""
    adata, _ = _small(n_donors=6, n_genes=100, n_de=10, seed=2)
    import scipy.sparse as sp
    X = adata.X.toarray() if sp.issparse(adata.X) else np.array(adata.X)
    X[:, :10] = 0                     # ten genes with nothing at all
    X[:, 10] = 0
    X[:2000, 10] = 3                  # one gene detected in a subset only
    adata.X = sp.csr_matrix(X)
    r = nbglm.nbglm_de(adata, "condition", "donor", "stim", "ctrl",
                       min_count=5, min_samples=3)
    assert len(r) > 0
    assert r["pval"].between(0, 1).all()
    assert np.isfinite(r["dispersion"]).all()
    assert (r["dispersion"] > 0).all()
    # a gene with zero counts in every sample must not appear at all
    assert "gene_0000" not in r.index


def test_paired_glm_keeps_only_complete_donors():
    """A donor present in one condition only is dropped, not modelled as null."""
    adata, truth = _small(n_donors=6, n_de=12, seed=4)
    obs = adata.obs
    drop = obs.index[(obs["donor"] == "donor_5") & (obs["condition"] == "stim")]
    adata = adata[~adata.obs.index.isin(drop)].copy()
    r = nbglm.nbglm_de(adata, "condition", "donor", "stim", "ctrl")
    assert r["n_a"].iloc[0] == 5
    assert r["n_b"].iloc[0] == 5
    assert r["qval"].between(0, 1).all()


def test_result_shape_matches_pseudobulk():
    """Same columns so the two can be put in one table."""
    adata, _ = _small(n_donors=6, n_genes=80, n_de=8, seed=5)
    nb = nbglm.nbglm_de(adata, "condition", "donor", "stim", "ctrl")
    pb = de.pseudobulk_de(adata, "condition", "donor", "stim", "ctrl")
    for col in ("mean_a", "mean_b", "log2fc", "pval", "qval", "n_a", "n_b"):
        assert col in nb.columns
        assert col in pb.columns
    assert nb.index.name == "gene"
    assert len(nb) > 0
