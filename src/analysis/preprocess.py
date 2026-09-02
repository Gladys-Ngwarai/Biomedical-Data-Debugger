import scanpy as sc


def preprocess_scrna(adata):
    """
    Normalize and prepare scRNA-seq data for clustering.
    """

    adata = adata.copy()

    # Normalize each cell to the same total count
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Log-transform expression values
    sc.pp.log1p(adata)

    # Identify highly variable genes
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        flavor="seurat"
    )

    # Keep highly variable genes
    adata = adata[:, adata.var["highly_variable"]].copy()

    # Scale expression
    sc.pp.scale(adata, max_value=10)

    return adata