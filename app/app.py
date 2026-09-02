import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


import streamlit as st
import scanpy as sc

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

from src.variants.load_vcf import load_vcf
from src.variants.variant_debugger import debug_variants


# =========================================================
# PAGE
# =========================================================

st.set_page_config(
    page_title="Biomedical Data Debugger",
    page_icon="🧬",
    layout="wide",
)


# =========================================================
# HEADER
# =========================================================

st.title("🧬 Biomedical Data Debugger")

st.subheader("Artifact vs. True Biology")

st.write(
    "Don't trust the biological conclusion "
    "until you debug the data."
)


# =========================================================
# MODULE SELECTION
# =========================================================

st.sidebar.header("1. Choose Module")

module = st.sidebar.radio(
    "Debugger module",
    [
        "Single-cell RNA-seq",
        "DNA Variant Debugger",
    ],
)


# =========================================================
# SINGLE-CELL RNA-SEQ
# =========================================================

if module == "Single-cell RNA-seq":

    st.sidebar.header("2. Input Data")

    input_method = st.sidebar.radio(
        "Choose dataset",
        [
            "Demo dataset",
            "Upload scRNA-seq (.h5ad)",
        ],
    )

    adata = None

    # -----------------------------------------------------
    # DEMO DATASET
    # -----------------------------------------------------

    if input_method == "Demo dataset":

        st.sidebar.success(
            "PBMC3k demonstration dataset selected."
        )

        if st.sidebar.button(
            "Run scRNA-seq Debugger"
        ):
            adata = load_pbmc3k()

    # -----------------------------------------------------
    # USER DATASET
    # -----------------------------------------------------

    else:

        uploaded_file = st.sidebar.file_uploader(
            "Upload an AnnData scRNA-seq dataset",
            type=["h5ad"],
        )

        if uploaded_file is not None:

            try:

                adata = sc.read_h5ad(
                    uploaded_file
                )

                st.sidebar.success(
                    "Dataset loaded successfully."
                )

                st.sidebar.write(
                    f"Cells: {adata.n_obs}"
                )

                st.sidebar.write(
                    f"Genes: {adata.n_vars}"
                )

            except Exception as e:

                st.sidebar.error(
                    "Could not read this .h5ad file."
                )

                st.sidebar.exception(e)

        run_uploaded = st.sidebar.button(
            "Run scRNA-seq Debugger"
        )

        if not run_uploaded:
            adata = None

    # =====================================================
    # RUN scRNA ANALYSIS
    # =====================================================

    if adata is not None:

        with st.spinner(
            "Running scRNA-seq debugger..."
        ):

            # Basic QC
            adata = calculate_basic_qc(adata)

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

            # Before analysis
            analysis_adata = preprocess_scrna(
                adata
            )

            analysis_adata = cluster_cells(
                analysis_adata
            )

            analysis_adata = annotate_clusters(
                analysis_adata
            )

            # Copy clustering information
            adata.obs["leiden"] = (
                analysis_adata.obs["leiden"]
            )

            adata.uns["cluster_cell_types"] = (
                analysis_adata.uns[
                    "cluster_cell_types"
                ]
            )

            # Ambient / lineage evidence
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

            # Evidence aggregation
            adata = aggregate_artifact_evidence(
                adata
            )

            # Before
            before_cells = adata.n_obs

            flagged_cells = int(
                adata.obs[
                    "requires_investigation"
                ].sum()
            )

            before_clusters = (
                adata.obs["leiden"].nunique()
            )

            # Correction
            corrected = (
                correct_suspicious_cells(
                    adata
                )
            )

            # Re-analysis
            corrected = reanalyze(
                corrected
            )

            after_cells = corrected.n_obs

            after_clusters = (
                corrected.obs["leiden"].nunique()
            )

            # Trust score
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

        # =================================================
        # RESULTS
        # =================================================

        st.success(
            "scRNA-seq analysis completed successfully."
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

        # =================================================
        # BEFORE / AFTER
        # =================================================

        st.header(
            "Before vs. After Debugging"
        )

        col1, col2 = st.columns(2)

        with col1:

            st.subheader("Before")

            st.write(
                f"Cells: {before_cells}"
            )

            st.write(
                f"Clusters: {before_clusters}"
            )

        with col2:

            st.subheader("After")

            st.write(
                f"Cells: {after_cells}"
            )

            st.write(
                f"Clusters: {after_clusters}"
            )

        # =================================================
        # BIOLOGICAL CONCLUSION
        # =================================================

        st.header(
            "Biological Conclusion"
        )

        st.info(
            interpretation
        )

        # =================================================
        # ARTIFACT EVIDENCE
        # =================================================

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

        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric(
                "Low Quality",
                low_quality
            )

        with col2:
            st.metric(
                "Potential Doublets",
                doublets
            )

        with col3:
            st.metric(
                "High Stress",
                stress
            )

        with col4:
            st.metric(
                "Lineage Evidence",
                ambient
            )

        st.caption(
            "Lineage-marker evidence is a screening "
            "signal and does not prove ambient RNA."
        )

        # =================================================
        # SUSPICIOUS CELLS
        # =================================================

        st.header(
            "Cells Requiring Investigation"
        )

        suspicious = adata.obs[
            adata.obs[
                "requires_investigation"
            ]
        ].copy()

        if len(suspicious) > 0:

            display_columns = [
                "leiden",
                "artifact_evidence_count",
                "low_quality_flag",
                "predicted_doublet",
                "stress_flag",
                "ambient_evidence",
            ]

            display_columns = [
                column
                for column in display_columns
                if column in suspicious.columns
            ]

            st.dataframe(
                suspicious[
                    display_columns
                ]
            )

        else:

            st.success(
                "No cells currently require investigation."
            )

    else:

        st.info(
            "Choose a dataset from the sidebar "
            "to begin debugging."
        )

        st.markdown(
            """
### Supported input

**`.h5ad` — AnnData single-cell RNA-seq dataset**

The debugger evaluates:

- Quality control
- Potential doublets
- Stress response
- Lineage-marker contamination evidence
- Clustering
- Cell-type signals
- Suspicious-cell evidence
- Correction
- Re-analysis
- Biological trust
            """
        )


# =========================================================
# DNA VARIANT DEBUGGER
# =========================================================

else:

    st.sidebar.header("2. Input Data")

    input_method = st.sidebar.radio(
        "Choose variant dataset",
        [
            "Demo VCF",
            "Upload VCF",
        ],
    )

    vcf = None

    # -----------------------------------------------------
    # DEMO VCF
    # -----------------------------------------------------

    if input_method == "Demo VCF":

        st.sidebar.success(
            "Synthetic demonstration VCF selected."
        )

        if st.sidebar.button(
            "Run DNA Variant Debugger"
        ):

            vcf = load_vcf(
                "data/raw/test_variants.vcf"
            )

    # -----------------------------------------------------
    # USER VCF
    # -----------------------------------------------------

    else:

        uploaded_vcf = st.sidebar.file_uploader(
            "Upload a VCF file",
            type=["vcf"],
        )

        if uploaded_vcf is not None:

            try:

                vcf = load_vcf(
                    uploaded_vcf
                )

                st.sidebar.success(
                    "VCF loaded successfully."
                )

                st.sidebar.write(
                    f"Variants: {len(vcf)}"
                )

            except Exception as e:

                st.sidebar.error(
                    "Could not read this VCF file."
                )

                st.sidebar.exception(e)

        run_vcf = st.sidebar.button(
            "Run DNA Variant Debugger"
        )

        if not run_vcf:
            vcf = None

    # =====================================================
    # RUN VARIANT ANALYSIS
    # =====================================================

    if vcf is not None:

        with st.spinner(
            "Running DNA variant debugger..."
        ):

            reports = debug_variants(
                vcf
            )

        st.success(
            "DNA variant analysis completed successfully."
        )

        # =================================================
        # SUMMARY
        # =================================================

        total_variants = len(reports)

        suspicious_variants = sum(
            report["requires_investigation"]
            for report in reports
        )

        clean_variants = (
            total_variants
            - suspicious_variants
        )

        st.header("Variant Debugging Summary")

        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric(
                "Variants",
                total_variants
            )

        with col2:
            st.metric(
                "Require Investigation",
                suspicious_variants
            )

        with col3:
            st.metric(
                "No Major Evidence",
                clean_variants
            )

        # =================================================
        # VARIANT EVIDENCE
        # =================================================

        st.header(
            "Variant Evidence"
        )

        table_rows = []

        for report in reports:

            table_rows.append(
                {
                    "Variant": report["variant"],
                    "Depth": report["depth"],
                    "Allele Fraction": report[
                        "allele_fraction"
                    ],
                    "Mapping Quality": report[
                        "mapping_quality"
                    ],
                    "Evidence Count": report[
                        "artifact_evidence"
                    ],
                    "Status": (
                        "Investigate"
                        if report[
                            "requires_investigation"
                        ]
                        else "No major evidence"
                    ),
                }
            )

        st.dataframe(
            table_rows,
            use_container_width=True
        )

        # =================================================
        # SUSPICIOUS VARIANTS
        # =================================================

        st.header(
            "Variants Requiring Investigation"
        )

        suspicious_reports = [
            report
            for report in reports
            if report[
                "requires_investigation"
            ]
        ]

        if suspicious_reports:

            for report in suspicious_reports:

                st.subheader(
                    report["variant"]
                )

                st.write(
                    f"Depth: {report['depth']}"
                )

                st.write(
                    "Allele fraction: "
                    f"{report['allele_fraction']}"
                )

                st.write(
                    "Mapping quality: "
                    f"{report['mapping_quality']}"
                )

                st.write(
                    "Evidence:"
                )

                for reason in report["reasons"]:

                    st.write(
                        f"• {reason}"
                    )

                st.divider()

        else:

            st.success(
                "No variants currently require investigation."
            )

        st.caption(
            "Variant evidence indicates potential technical "
            "problems; it does not by itself prove that a "
            "variant is false."
        )

    else:

        st.info(
            "Choose a VCF dataset from the sidebar "
            "to begin debugging."
        )

        st.markdown(
            """
### Supported input

**`.vcf` — Variant Call Format**

The debugger evaluates available variant evidence such as:

- Sequencing depth
- Alternate-allele support
- Allele fraction
- Mapping quality
- Suspicious evidence patterns

Variants are flagged for investigation rather than
automatically declared false.
            """
        )