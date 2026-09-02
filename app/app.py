import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
    
import streamlit as st

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

from src.analysis.trust_score import (
    calculate_trust_score,
    interpret_trust_score,
)


st.set_page_config(
    page_title="Biomedical Data Debugger",
    page_icon="🧬",
    layout="wide",
)


# -------------------------
# TITLE
# -------------------------

st.title("🧬 Biomedical Data Debugger")

st.subheader(
    "Artifact vs. True Biology"
)

st.write(
    "Don't trust the biological conclusion "
    "until you debug the data."
)


# -------------------------
# SIDEBAR
# -------------------------

st.sidebar.header("Analysis")

run_analysis = st.sidebar.button(
    "Run scRNA-seq Debugger"
)


# -------------------------
# ANALYSIS
# -------------------------

if run_analysis:

    with st.spinner(
        "Running biomedical data debugger..."
    ):

        # Load
        adata = load_pbmc3k()

        # QC
        adata = calculate_basic_qc(
            adata
        )

        adata = flag_low_quality_cells(
            adata
        )

        # Doublets
        adata = detect_doublets(
            adata
        )

        # Stress
        adata = calculate_stress_score(
            adata
        )

        # -------------------------
        # BEFORE
        # -------------------------

        analysis_adata = preprocess_scrna(
            adata
        )

        analysis_adata = cluster_cells(
            analysis_adata
        )

        analysis_adata = annotate_clusters(
            analysis_adata
        )

        adata.obs["leiden"] = (
            analysis_adata.obs["leiden"]
        )

        adata.uns["cluster_cell_types"] = (
            analysis_adata.uns[
                "cluster_cell_types"
            ]
        )

        # Ambient evidence
        analysis_adata = (
            calculate_ambient_rna_evidence(
                analysis_adata
            )
        )

        for column in analysis_adata.obs.columns:

            if column.startswith(
                "ambient_"
            ):

                adata.obs[column] = (
                    analysis_adata.obs[column]
                )

        adata.uns["ambient_marker_scores"] = (
            analysis_adata.uns[
                "ambient_marker_scores"
            ]
        )

        # Aggregate evidence
        adata = aggregate_artifact_evidence(
            adata
        )

        # -------------------------
        # METRICS
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

        corrected = (
            correct_suspicious_cells(
                adata
            )
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
        # TRUST
        # -------------------------

        trust_score = (
            calculate_trust_score(
                before_cells,
                after_cells,
                before_clusters,
                after_clusters,
                flagged_cells,
            )
        )

        interpretation = (
            interpret_trust_score(
                trust_score
            )
        )

    # -------------------------
    # DASHBOARD
    # -------------------------

    st.success(
        "Analysis completed successfully."
    )

    st.header("Debugging Summary")

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric(
            "Cells Before",
            before_cells
        )

    with col2:
        st.metric(
            "Suspicious Cells",
            flagged_cells
        )

    with col3:
        st.metric(
            "Cells After",
            after_cells
        )

    with col4:
        st.metric(
            "Trust Score",
            f"{trust_score}/100"
        )

    # -------------------------
    # BEFORE / AFTER
    # -------------------------

    st.header(
        "Before vs. After Debugging"
    )

    col1, col2 = st.columns(2)

    with col1:

        st.subheader(
            "Before"
        )

        st.write(
            f"Cells: {before_cells}"
        )

        st.write(
            f"Clusters: {before_clusters}"
        )

    with col2:

        st.subheader(
            "After"
        )

        st.write(
            f"Cells: {after_cells}"
        )

        st.write(
            f"Clusters: {after_clusters}"
        )

    # -------------------------
    # CONCLUSION
    # -------------------------

    st.header(
        "Biological Conclusion"
    )

    st.info(
        interpretation
    )

    # -------------------------
    # ARTIFACT BREAKDOWN
    # -------------------------

    st.header(
        "Artifact Evidence"
    )

    low_quality = int(
        adata.obs[
            "low_quality_flag"
        ].sum()
    )

    doublets = int(
        adata.obs[
            "predicted_doublet"
        ].sum()
    )

    stress = int(
        adata.obs[
            "stress_flag"
        ].sum()
    )

    ambient = int(
        adata.obs[
            "ambient_evidence"
        ].sum()
    )

    st.write(
        f"**Low-quality cells:** {low_quality}"
    )

    st.write(
        f"**Potential doublets:** {doublets}"
    )

    st.write(
        f"**High-stress cells:** {stress}"
    )

    st.write(
        f"**Lineage-marker evidence:** {ambient}"
    )

    st.caption(
        "Lineage-marker evidence is a screening "
        "signal and does not prove ambient RNA."
    )

else:

    st.info(
        "Click **Run scRNA-seq Debugger** "
        "to analyze the demonstration dataset."
    )
