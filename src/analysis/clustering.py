import scanpy as sc


def cluster_cells(adata):
    """
    Perform PCA, neighborhood analysis, UMAP, and Leiden clustering.
    """

    sc.tl.pca(adata, svd_solver="arpack")

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