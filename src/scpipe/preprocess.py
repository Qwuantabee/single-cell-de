"""Quality control and the standard clustering pipeline.

Nothing here is novel -- it is the scanpy workflow that most single-cell
analyses run. The parts worth attention are the QC thresholds, because they
are the step most often done by copying numbers out of a tutorial, and they
are dataset dependent.

The functions return the object plus a record of what was removed, so the
filtering is auditable instead of being a silent step that changed the cell
count by 30%.
"""

import numpy as np


def qc_metrics(adata, mito_prefix=("MT-", "mt-")):
    """Per-cell counts, genes detected, and mitochondrial fraction.

    High mitochondrial fraction usually means a cell whose membrane broke
    during dissociation and lost cytoplasmic RNA, leaving the mitochondria
    behind. It is a dying-cell signal, not a biological one, which is why
    it gets filtered -- but the threshold depends on the tissue, and in
    cardiac or muscle samples a hard 5% cut removes real cells.
    """
    import scanpy as sc

    adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True,
                               percent_top=None, log1p=False)
    return adata


def filter_cells(adata, min_genes=200, max_mt_pct=15, max_counts=None,
                 verbose=True):
    """Drop low-quality cells, reporting what went.

    min_genes=200 removes empty or near-empty droplets.
    max_mt_pct=15 is deliberately looser than the 5% often used; PBMCs
    tolerate it and a tight cut preferentially removes monocytes, which
    have higher baseline mitochondrial content than lymphocytes. Filtering
    on a threshold that correlates with cell type silently changes the
    composition of the experiment.
    """
    before = adata.n_obs
    keep = np.ones(adata.n_obs, dtype=bool)
    reasons = {}

    if "n_genes_by_counts" in adata.obs:
        bad = adata.obs["n_genes_by_counts"].to_numpy() < min_genes
        reasons[f"fewer than {min_genes} genes"] = int(bad.sum())
        keep &= ~bad

    if "pct_counts_mt" in adata.obs:
        bad = adata.obs["pct_counts_mt"].to_numpy() > max_mt_pct
        reasons[f"over {max_mt_pct}% mitochondrial"] = int(bad.sum())
        keep &= ~bad

    if max_counts is not None and "total_counts" in adata.obs:
        bad = adata.obs["total_counts"].to_numpy() > max_counts
        reasons[f"over {max_counts} counts (likely doublets)"] = int(bad.sum())
        keep &= ~bad

    out = adata[keep].copy()
    if verbose:
        print(f"QC: {before:,} -> {out.n_obs:,} cells "
              f"({100 * (1 - out.n_obs / before):.1f}% removed)")
        for reason, n in reasons.items():
            print(f"    {n:>6,} {reason}")
    return out


def standard_pipeline(adata, n_hvg=2000, n_pcs=30, resolution=1.0,
                      batch_key=None, seed=0):
    """Normalise, select variable genes, PCA, neighbours, Leiden, UMAP.

    Raw counts are kept in .layers['counts'] because the differential
    expression in de.py needs them -- aggregating normalised values and then
    normalising again is a mistake that is easy to make and hard to see.
    """
    import scanpy as sc

    adata = adata.copy()
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg,
                                batch_key=batch_key)
    adata = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs, random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata, random_state=seed)
    return adata


# Canonical PBMC markers. Used for annotating clusters by hand, which is
# still the most reliable approach when the cell types are well known.
PBMC_MARKERS = {
    "CD4 T": ["IL7R", "CD3D", "CCR7"],
    "CD8 T": ["CD8A", "CD8B", "GZMK"],
    "NK": ["GNLY", "NKG7", "KLRD1"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "CD14 Mono": ["CD14", "LYZ", "S100A9"],
    "FCGR3A Mono": ["FCGR3A", "MS4A7"],
    "DC": ["FCER1A", "CST3"],
    "Megakaryocyte": ["PPBP", "PF4"],
}


def score_markers(adata, markers=None):
    """Mean expression of each marker set per cluster.

    Returns a cluster x cell-type table. Annotating from this rather than
    from a single marker gene matters because most individual markers are
    expressed in more than one lineage.
    """
    import pandas as pd
    import scanpy as sc

    markers = markers or PBMC_MARKERS
    source = adata.raw.to_adata() if adata.raw is not None else adata

    rows = {}
    for label, genes in markers.items():
        present = [g for g in genes if g in source.var_names]
        if not present:
            continue
        sc.tl.score_genes(source, present, score_name=f"_score_{label}")
        rows[label] = source.obs[f"_score_{label}"].to_numpy()

    scores = pd.DataFrame(rows, index=adata.obs_names)
    scores["cluster"] = adata.obs["leiden"].to_numpy()
    return scores.groupby("cluster").mean()
