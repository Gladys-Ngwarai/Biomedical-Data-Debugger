import scanpy as sc


def plot_clusters(adata, output_path="figures/umap_clusters.png"):
    """
    Save a UMAP showing Leiden clusters.
    """

    sc.pl.umap(
        adata,
        color="leiden",
        show=False
    )

    import matplotlib.pyplot as plt

    plt.savefig(output_path, bbox_inches="tight")
    plt.close()