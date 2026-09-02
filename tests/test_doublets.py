from src.analysis.load_data import load_pbmc3k
from src.artifacts.doublets import detect_doublets


def main():
    # Load raw-count dataset
    adata = load_pbmc3k()

    print("Dataset loaded")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")

    # Detect potential doublets
    adata = detect_doublets(adata)

    total_cells = adata.n_obs
    doublet_cells = adata.obs["predicted_doublet"].sum()

    percentage = (doublet_cells / total_cells) * 100

    print("\nDoublet detection complete")
    print(f"Potential doublets: {doublet_cells}")
    print(f"Percentage flagged: {percentage:.2f}%")


if __name__ == "__main__":
    main()