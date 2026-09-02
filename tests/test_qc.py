from pathlib import Path

from src.analysis.load_data import load_pbmc3k
from src.qc.basic_qc import calculate_basic_qc
from src.qc.qc_plots import plot_qc_metrics


def main():
    # Load dataset
    adata = load_pbmc3k()

    print("Dataset loaded")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")

    # Calculate QC metrics
    adata = calculate_basic_qc(adata)

    print("\nQC metrics calculated")

    # Create output directory
    output_dir = Path("data/processed/qc")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate QC plots
    plot_qc_metrics(
        adata,
        output_path=str(output_dir)
    )

    print(f"\nQC plots generated in: {output_dir}")


if __name__ == "__main__":
    main()