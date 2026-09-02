import scanpy as sc


def find_cluster_markers(adata, top_n=5):
    """
    Identify marker genes for each Leiden cluster.
    """

    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon"
    )

    for cluster in adata.obs["leiden"].cat.categories:
        genes = adata.uns["rank_genes_groups"]["names"][cluster][:top_n]

        print(f"\nCluster {cluster}:")
        print(", ".join(genes))

    return adata