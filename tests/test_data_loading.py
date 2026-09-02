from src.analysis.load_data import load_pbmc3k


def main():
    adata = load_pbmc3k()

    print("Dataset loaded successfully")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")


if __name__ == "__main__":
    main()