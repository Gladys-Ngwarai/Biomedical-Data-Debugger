def flag_low_quality_cells(
    adata,
    min_genes=200,
    max_genes=None,
    max_mito_percent=20.0,
):
    """
    Flag potentially low-quality cells using basic QC metrics.

    Cells are flagged but NOT removed.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing QC metrics.
    min_genes : int
        Minimum number of detected genes expected.
    max_genes : int or None
        Optional maximum number of detected genes.
    max_mito_percent : float
        Maximum mitochondrial percentage allowed.

    Returns
    -------
    anndata.AnnData
        AnnData object with a low_quality_flag column.
    """

    flag = adata.obs["n_genes_by_counts"] < min_genes

    if max_genes is not None:
        flag = flag | (
            adata.obs["n_genes_by_counts"] > max_genes
        )

    flag = flag | (
        adata.obs["pct_counts_mt"] > max_mito_percent
    )

    adata.obs["low_quality_flag"] = flag

    return adata