from src.analysis.load_data import load_pbmc3k
from src.analysis.preprocess import preprocess_scrna
from src.analysis.clustering import cluster_cells
from src.analysis.annotation import find_cluster_markers


def main():
    adata = load_pbmc3k()

    adata = preprocess_scrna(adata)
    adata = cluster_cells(adata)
    adata = find_cluster_markers(adata)

    print("\nCluster marker analysis complete")


if __name__ == "__main__":
    main()