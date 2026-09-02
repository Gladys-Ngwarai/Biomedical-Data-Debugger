import matplotlib.pyplot as plt
import scanpy as sc


def create_qc_figure(adata):
    """Create basic QC visualizations."""

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(15, 4),
    )

    axes[0].hist(
        adata.obs["n_genes_by_counts"],
        bins=40,
    )

    axes[0].set_title(
        "Genes per cell"
    )

    axes[0].set_xlabel(
        "Genes"
    )

    axes[0].set_ylabel(
        "Cells"
    )

    axes[1].hist(
        adata.obs["total_counts"],
        bins=40,
    )

    axes[1].set_title(
        "Total counts"
    )

    axes[1].set_xlabel(
        "Counts"
    )

    axes[2].hist(
        adata.obs["pct_counts_mt"],
        bins=40,
    )

    axes[2].set_title(
        "Mitochondrial percentage"
    )

    axes[2].set_xlabel(
        "Percent"
    )

    fig.tight_layout()

    return fig


def create_umap_figure(
    adata,
    color="leiden",
):
    """Create a UMAP figure."""

    sc.pl.umap(
        adata,
        color=color,
        show=False,
    )

    fig = plt.gcf()

    fig.tight_layout()

    return fig