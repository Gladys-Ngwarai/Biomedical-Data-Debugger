import scanpy as sc


def plot_qc_metrics(adata, output_path=None):
    """
    Generate basic QC plots for an scRNA-seq dataset.

    Parameters
    ----------
    adata : anndata.AnnData
        Dataset with QC metrics already calculated.
    output_path : str, optional
        Directory where plots should be saved.
    """

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False
    )

    if output_path:
        sc.settings.figdir = output_path

    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        show=False
    )

    sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_mt",
        show=False
    )