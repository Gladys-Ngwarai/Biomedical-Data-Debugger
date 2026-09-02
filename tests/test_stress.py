from src.analysis.load_data import load_pbmc3k
from src.artifacts.stress import calculate_stress_score


def main():
    adata = load_pbmc3k()

    adata = calculate_stress_score(adata)

    print("Stress analysis complete")
    print(f"Mean stress score: {adata.obs['stress_score'].mean():.4f}")
    print(f"Maximum stress score: {adata.obs['stress_score'].max():.4f}")


if __name__ == "__main__":
    main()