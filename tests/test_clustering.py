from src.analysis.load_data import load_pbmc3k
from src.analysis.preprocess import preprocess_scrna
from src.analysis.clustering import cluster_cells
from src.analysis.visualization import plot_clusters

def main():
    adata = load_pbmc3k()

    print("Dataset loaded")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")

    adata = preprocess_scrna(adata)

    print("\nPreprocessing complete")
    print(f"Genes after HVG selection: {adata.n_vars}")

    adata = cluster_cells(adata)
    plot_clusters(adata)
    print("\nClustering complete")
    print(f"Number of clusters: {adata.obs['leiden'].nunique()}")
    print("Clusters:", sorted(adata.obs["leiden"].unique()))


if __name__ == "__main__":
    main()