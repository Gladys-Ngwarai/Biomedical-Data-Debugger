from src.analysis.load_data import load_pbmc3k
from src.qc.basic_qc import calculate_basic_qc
from src.analysis.preprocess import preprocess_scrna
from src.analysis.clustering import cluster_cells
from src.analysis.cell_types import annotate_clusters

from src.artifacts.low_quality import flag_low_quality_cells
from src.artifacts.doublets import detect_doublets
from src.artifacts.stress import calculate_stress_score
from src.artifacts.ambient_rna import calculate_ambient_rna_evidence
from src.artifacts.evidence import aggregate_artifact_evidence

from src.analysis.correction import (
    correct_suspicious_cells,
    reanalyze,
)

from src.analysis.comparison import compare_before_after
from src.analysis.trust_score import (
    calculate_trust_score,
    interpret_trust_score,
)


def main():

    # -------------------------
    # LOAD DATA
    # -------------------------

    adata = load_pbmc3k()

    # -------------------------
    # QC + ARTIFACT DETECTION
    # -------------------------

    adata = calculate_basic_qc(adata)

    adata = flag_low_quality_cells(adata)

    adata = detect_doublets(adata)

    adata = calculate_stress_score(adata)

    # -------------------------
    # BEFORE ANALYSIS
    # -------------------------

    analysis_adata = preprocess_scrna(adata)

    analysis_adata = cluster_cells(
        analysis_adata
    )

    analysis_adata = annotate_clusters(
        analysis_adata
    )

    # Copy cluster information back.
    adata.obs["leiden"] = (
        analysis_adata.obs["leiden"]
    )

    adata.uns["cluster_cell_types"] = (
        analysis_adata.uns["cluster_cell_types"]
    )

    # -------------------------
    # AMBIENT EVIDENCE
    # -------------------------

    analysis_adata = calculate_ambient_rna_evidence(
        analysis_adata
    )

    for column in [
        "ambient_T_cell",
        "ambient_NK_cell",
        "ambient_B_cell",
        "ambient_Monocyte",
        "ambient_Platelet",
        "ambient_Dendritic",
    ]:

        if column in analysis_adata.obs:

            adata.obs[column] = (
                analysis_adata.obs[column]
            )

    adata.uns["ambient_marker_scores"] = (
        analysis_adata.uns[
            "ambient_marker_scores"
        ]
    )

    # -------------------------
    # EVIDENCE AGGREGATION
    # -------------------------

    adata = aggregate_artifact_evidence(
        adata
    )

    # -------------------------
    # BEFORE
    # -------------------------

    before_cells = adata.n_obs

    flagged_cells = int(
        adata.obs[
            "requires_investigation"
        ].sum()
    )

    before_clusters = (
        adata.obs["leiden"].nunique()
    )

    # -------------------------
    # CORRECTION
    # -------------------------

    corrected = correct_suspicious_cells(
        adata
    )

    # -------------------------
    # AFTER
    # -------------------------

    corrected = reanalyze(
        corrected
    )

    after_cells = corrected.n_obs

    after_clusters = (
        corrected.obs["leiden"].nunique()
    )

    # -------------------------
    # COMPARISON
    # -------------------------

    comparison = compare_before_after(
        adata,
        corrected
    )

    # -------------------------
    # TRUST SCORE
    # -------------------------

    trust_score = calculate_trust_score(
        before_cells=before_cells,
        after_cells=after_cells,
        before_clusters=before_clusters,
        after_clusters=after_clusters,
        flagged_cells=flagged_cells,
    )

    interpretation = interpret_trust_score(
        trust_score
    )

    # -------------------------
    # REPORT
    # -------------------------

    print()
    print("Biomedical Data Debugger")
    print("=" * 50)

    print()
    print("BEFORE DEBUGGING")
    print("-" * 30)
    print(f"Cells: {comparison['before_cells']}")
    print(
        f"Clusters: "
        f"{comparison['before_clusters']}"
    )

    print()
    print("ARTIFACT INVESTIGATION")
    print("-" * 30)
    print(
        f"Suspicious cells: "
        f"{flagged_cells}"
    )

    print()
    print("CORRECTION")
    print("-" * 30)
    print(
        f"Cells removed: "
        f"{comparison['cells_removed']}"
    )

    print()
    print("AFTER DEBUGGING")
    print("-" * 30)
    print(f"Cells: {comparison['after_cells']}")
    print(
        f"Clusters: "
        f"{comparison['after_clusters']}"
    )

    print()
    print("BIOLOGICAL TRUST")
    print("-" * 30)
    print(
        f"Trust score: "
        f"{trust_score}/100"
    )

    print()
    print("Conclusion:")
    print(interpretation)

    print()
    print("Re-analysis complete.")