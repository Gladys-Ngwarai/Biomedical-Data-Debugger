from src.analysis.load_data import load_pbmc3k
from src.qc.basic_qc import calculate_basic_qc
from src.analysis.preprocess import preprocess_scrna
from src.analysis.clustering import cluster_cells
from src.analysis.cell_types import annotate_clusters

from src.artifacts.low_quality import flag_low_quality_cells
from src.artifacts.doublets import detect_doublets
from src.artifacts.stress import calculate_stress_score
from src.artifacts.evidence import aggregate_artifact_evidence
from src.artifacts.debug_report import generate_debug_report


def main():
    adata = load_pbmc3k()

    # Artifact evidence
    adata = calculate_basic_qc(adata)
    adata = flag_low_quality_cells(adata)
    adata = detect_doublets(adata)
    adata = calculate_stress_score(adata)

    # Biological context
    analysis_adata = preprocess_scrna(adata)
    analysis_adata = cluster_cells(analysis_adata)
    analysis_adata = annotate_clusters(analysis_adata)

    # Transfer cluster identities back
    adata.obs["leiden"] = analysis_adata.obs["leiden"]
    adata.uns["cluster_cell_types"] = (
        analysis_adata.uns["cluster_cell_types"]
    )

    # Combine evidence
    adata = aggregate_artifact_evidence(adata)

    reports = generate_debug_report(adata)

    print("\nBiomedical Data Debugger Report")
    print("=" * 40)

    for report in reports:
        cell_type = adata.uns["cluster_cell_types"].get(
            report["cluster"],
            "Unknown"
        )

        print(f"\nCell: {report['cell']}")
        print(f"Cluster: {report['cluster']}")
        print(f"Likely cell type: {cell_type}")
        print(
            f"Evidence count: "
            f"{report['artifact_evidence_count']}"
        )
        print("Reasons: " + ", ".join(report["reasons"]))
        print(
            f"Recommendation: "
            f"{report['recommendation']}"
        )


if __name__ == "__main__":
    main()