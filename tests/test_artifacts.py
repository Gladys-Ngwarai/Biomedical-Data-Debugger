from src.analysis.load_data import load_pbmc3k
from src.qc.basic_qc import calculate_basic_qc
from src.artifacts.low_quality import flag_low_quality_cells


def main():
    # Load data
    adata = load_pbmc3k()

    # Calculate QC metrics
    adata = calculate_basic_qc(adata)

    # Flag suspicious cells
    adata = flag_low_quality_cells(adata)

    total_cells = adata.n_obs
    suspicious_cells = adata.obs["low_quality_flag"].sum()

    print("Artifact detection complete")
    print(f"Total cells: {total_cells}")
    print(f"Potentially low-quality cells: {suspicious_cells}")

    percentage = (suspicious_cells / total_cells) * 100

    print(f"Percentage flagged: {percentage:.2f}%")


if __name__ == "__main__":
    main()