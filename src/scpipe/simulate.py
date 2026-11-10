"""Simulated single-cell data with a known answer.

The real-data null experiment (shuffling donors) shows that the per-cell test
produces false positives, but it cannot show how many true positives each
method finds, because in real data the truth is unknown.

This simulator makes both measurable. Counts are drawn from a negative
binomial, which is the standard model for RNA-seq overdispersion, with a
per-donor random effect so that cells from the same donor are correlated.
That correlation is the entire reason the per-cell test fails: it is what
makes 1000 cells from one donor worth far less than 1000 cells from 1000
donors, and the per-cell test assumes they are worth the same.

Setting donor_sd=0 removes the random effect and the per-cell test becomes
valid again, which is a useful way to confirm the diagnosis rather than just
observing the symptom.
"""

import numpy as np


def simulate(n_genes=500, n_donors=8, cells_per_donor=400, n_de=50,
             effect=1.0, donor_sd=0.5, sample_sd=0.3, dispersion=0.3,
             seed=0, base_mean=5.0, cell_count_variation=0.0):
    """Two conditions, paired within donor.

    n_de genes are given a real effect; the rest are null. Returns an
    AnnData with raw counts and obs columns donor / condition, plus the
    set of truly-changed gene names.

    donor_sd: per-donor shift applied identically to both of that donor's
        samples. Captures stable between-individual differences.

    sample_sd: per-(donor, condition) shift, independent between a donor's
        two samples. This is the one that matters. It represents everything
        that differs between two aliquots of the same donor's cells --
        processing day, capture efficiency, how that individual's cells
        happened to respond -- and it is the variance component that a
        per-cell test cannot see.

        The first version of this simulator had sample_sd=0, so a donor's
        two samples differed only by Poisson noise over hundreds of cells.
        Under that model the per-cell test looked perfectly calibrated,
        which is not what happens on real data. The simulator was wrong,
        not the literature.

    cell_count_variation: log-scale spread in cells recovered per sample.
        Real experiments vary several-fold; equal counts make the pooled
        test look better behaved than it is.

    dispersion: negative binomial dispersion; variance = mu + dispersion*mu^2.
    """
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)

    gene_names = [f"gene_{i:04d}" for i in range(n_genes)]
    de_idx = rng.choice(n_genes, size=n_de, replace=False)
    de_genes = {gene_names[i] for i in de_idx}

    # baseline expression spans a realistic dynamic range
    base = np.exp(rng.normal(np.log(base_mean), 1.2, size=n_genes))

    log_fc = np.zeros(n_genes)
    log_fc[de_idx] = effect * rng.choice([-1, 1], size=n_de)

    blocks, obs_rows = [], []
    for d in range(n_donors):
        donor_shift = rng.normal(0, donor_sd, size=n_genes)
        for cond in ("ctrl", "stim"):
            # independent between this donor's two samples -- see docstring
            sample_shift = rng.normal(0, sample_sd, size=n_genes)
            mu = base * np.exp(donor_shift + sample_shift
                               + (log_fc if cond == "stim" else 0))
            mu = np.clip(mu, 1e-3, None)
            n = cells_per_donor
            if cell_count_variation > 0:
                n = int(cells_per_donor
                        * np.exp(rng.normal(0, cell_count_variation)))
                n = max(n, 25)

            # negative binomial via gamma-poisson
            shape = 1.0 / dispersion
            gamma = rng.gamma(shape, mu / shape, size=(n, n_genes))
            counts = rng.poisson(gamma)

            blocks.append(counts)
            obs_rows += [(f"donor_{d}", cond)] * n

    X = sp.csr_matrix(np.vstack(blocks).astype(np.float32))
    obs = pd.DataFrame(obs_rows, columns=["donor", "condition"])
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    var = pd.DataFrame(index=gene_names)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.uns["true_de_genes"] = sorted(de_genes)
    return adata, de_genes


def confusion(result, true_genes, alpha=0.05):
    """True/false positives and negatives against the known answer."""
    called = set(result.index[result["qval"] < alpha])
    tested = set(result.index)
    truth = set(true_genes) & tested

    tp = len(called & truth)
    fp = len(called - truth)
    fn = len(truth - called)
    tn = len(tested) - tp - fp - fn

    return {
        "n_called": len(called),
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "power": tp / len(truth) if truth else np.nan,
        # observed false discovery proportion: of everything we called,
        # what fraction was wrong. Should sit near alpha.
        "fdp": fp / len(called) if called else 0.0,
    }
