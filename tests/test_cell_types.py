from src.analysis.load_data import load_pbmc3k
from src.analysis.preprocess import preprocess_scrna
from src.analysis.clustering import cluster_cells
from src.analysis.cell_types import annotate_clusters


def main():
    adata = load_pbmc3k()

    adata = preprocess_scrna(adata)
    adata = cluster_cells(adata)
    adata = annotate_clusters(adata)

    print("\nCluster annotation complete")

    for cluster, cell_type in adata.uns["cluster_cell_types"].items():
        print(f"Cluster {cluster}: {cell_type}")


if __name__ == "__main__":
    main()