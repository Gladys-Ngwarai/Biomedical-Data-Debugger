from src.analysis.load_data import load_pbmc3k
from src.qc.basic_qc import calculate_basic_qc
from src.artifacts.low_quality import flag_low_quality_cells
from src.artifacts.doublets import detect_doublets
from src.artifacts.stress import calculate_stress_score
from src.artifacts.evidence import aggregate_artifact_evidence


def main():
    adata = load_pbmc3k()

    adata = calculate_basic_qc(adata)
    adata = flag_low_quality_cells(adata)
    adata = detect_doublets(adata)
    adata = calculate_stress_score(adata)

    adata = aggregate_artifact_evidence(adata)

    print("\nEvidence aggregation complete")

    print(
        "Cells requiring investigation:",
        adata.obs["requires_investigation"].sum()
    )

    print("\nEvidence count:")
    print(
        adata.obs["artifact_evidence_count"]
        .value_counts()
        .sort_index()
    )

    print("\nEvidence combinations:")

    flagged = adata.obs[adata.obs["requires_investigation"]]

    print(
        flagged[
            [
                "low_quality_flag",
                "predicted_doublet",
                "stress_score",
                "artifact_evidence_count",
            ]
        ].to_string()
    )


if __name__ == "__main__":
    main()