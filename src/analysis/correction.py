import scanpy as sc


def correct_suspicious_cells(adata):
    """
    Create a corrected dataset by excluding cells that have
    multiple independent artifact signals.

    The original dataset is preserved.
    """

    if "requires_investigation" not in adata.obs:
        raise ValueError(
            "Run artifact evidence aggregation first."
        )

    corrected = adata[
        ~adata.obs["requires_investigation"]
    ].copy()

    return corrected


def reanalyze(adata):
    """
    Re-run normalization, clustering, and UMAP after correction.
    """

    adata = adata.copy()

    sc.pp.normalize_total(
        adata,
        target_sum=1e4
    )

    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        flavor="seurat"
    )

    adata = adata[
        :,
        adata.var["highly_variable"]
    ].copy()

    sc.pp.scale(
        adata,
        max_value=10
    )

    sc.tl.pca(
        adata,
        svd_solver="arpack"
    )

    sc.pp.neighbors(
        adata,
        n_neighbors=10,
        n_pcs=30
    )

    sc.tl.umap(adata)

    sc.tl.leiden(
        adata,
        resolution=0.5
    )

    return adata