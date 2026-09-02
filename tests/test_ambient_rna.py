from src.analysis.load_data import load_pbmc3k
from src.analysis.preprocess import preprocess_scrna
from src.analysis.clustering import cluster_cells
from src.artifacts.ambient_rna import calculate_ambient_rna_evidence


def main():
    adata = load_pbmc3k()

    adata = preprocess_scrna(adata)
    adata = cluster_cells(adata)
    adata = calculate_ambient_rna_evidence(adata)

    print("\nAmbient RNA evidence analysis complete")
    print("Scores calculated:")

    for cell_type, score_name in adata.uns["ambient_marker_scores"].items():
        print(f"- {cell_type}: {score_name}")


if __name__ == "__main__":
    main()