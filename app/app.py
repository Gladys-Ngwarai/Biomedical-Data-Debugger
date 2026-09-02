import sys
from pathlib import Path
import tempfile

import pandas as pd
import streamlit as st

# Make project root importable when running:
# streamlit run app/app.py
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

# -------------------------------------------------------------------
# Imports
# -------------------------------------------------------------------

from src.analysis.load_data import (
    load_pbmc3k,
    load_scrna_file,
    load_10x_mtx,
    validate_scrna,
)
from src.qc.basic_qc import calculate_basic_qc

from src.artifacts.low_quality import flag_low_quality_cells
from src.artifacts.doublets import detect_doublets
from src.artifacts.ambient_rna import calculate_ambient_rna_evidence
from src.artifacts.stress import calculate_stress_score
from src.artifacts.evidence import aggregate_artifact_evidence
from src.artifacts.debug_report import generate_debug_report

from src.analysis.preprocess import preprocess_scrna
from src.analysis.clustering import cluster_cells
from src.analysis.cell_types import annotate_clusters

from src.analysis.correction import (
    correct_suspicious_cells,
    reanalyze,
)
from src.analysis.comparison import compare_before_after
from src.analysis.trust_score import (
    calculate_trust_score,
    interpret_trust_score,
)

from src.analysis.reporting import (
    build_scrna_summary,
    reports_to_dataframe,
)

from src.analysis.plots import (
    create_qc_figure,
    create_umap_figure,
)

from src.variants.load_vcf import load_vcf
from src.variants.variant_debugger import debug_variants


# -------------------------------------------------------------------
# Page configuration
# -------------------------------------------------------------------

st.set_page_config(
    page_title="Biomedical Data Debugger",
    page_icon="🧬",
    layout="wide",
)


# -------------------------------------------------------------------
# Header
# -------------------------------------------------------------------

st.title("🧬 Biomedical Data Debugger")

st.caption(
    "Don't trust the biological conclusion until you debug the data."
)

st.markdown(
    """
**Artifact vs. True Biology**

A modular biomedical data debugger that investigates whether
apparent biological findings are supported by the data or may
be influenced by technical artifacts.
"""
)


# -------------------------------------------------------------------
# Sidebar
# -------------------------------------------------------------------

st.sidebar.title("Debugger")

module = st.sidebar.selectbox(
    "Analysis module",
    [
        "Single-cell RNA-seq",
        "DNA Variant Debugger",
    ],
)

st.sidebar.markdown("---")

st.sidebar.caption(
    "Evidence is flagged for investigation. "
    "The system does not automatically declare biological findings false."
)


# ===================================================================
# SINGLE-CELL RNA-SEQ
# ===================================================================

if module == "Single-cell RNA-seq":

    st.header("Single-cell RNA-seq Debugger")

    st.write(
        "Investigate whether apparent cellular structure is biologically "
        "coherent or potentially driven by technical artifacts."
    )

    # ---------------------------------------------------------------
    # Dataset input
    # ---------------------------------------------------------------

    st.subheader("1. Dataset")

    uploaded_file = st.file_uploader(
        "Upload an scRNA-seq dataset",
        type=["h5ad", "h5"],
        help=(
            "Supported formats: AnnData .h5ad and "
            "10x Genomics .h5."
        ),
    )

    st.caption(
        "You can also provide a 10x Matrix Market dataset below."
    )

    with st.expander("Advanced: 10x Matrix Market input"):

        matrix_file = st.file_uploader(
            "matrix.mtx or matrix.mtx.gz",
            type=["mtx", "gz"],
            key="matrix_file",
        )

        barcodes_file = st.file_uploader(
            "barcodes.tsv or barcodes.tsv.gz",
            type=["tsv", "gz"],
            key="barcodes_file",
        )

        features_file = st.file_uploader(
            "features.tsv, genes.tsv, or compressed equivalent",
            type=["tsv", "gz"],
            key="features_file",
        )

    # ---------------------------------------------------------------
    # Debugging settings
    # ---------------------------------------------------------------

    st.subheader("2. Debugging Settings")

    col1, col2, col3 = st.columns(3)

    with col1:
        min_genes = st.number_input(
            "Minimum genes per cell",
            min_value=0,
            value=200,
            step=25,
        )

    with col2:
        max_mito = st.number_input(
            "Maximum mitochondrial %",
            min_value=0.0,
            max_value=100.0,
            value=20.0,
            step=1.0,
        )

    with col3:
        doublet_rate = st.number_input(
            "Expected doublet rate",
            min_value=0.0,
            max_value=0.50,
            value=0.05,
            step=0.01,
        )

    # ---------------------------------------------------------------
    # Run
    # ---------------------------------------------------------------

    run_debugger = st.button(
        "Run scRNA-seq Debugger",
        type="primary",
        use_container_width=True,
    )

    if run_debugger:

        try:

            # =======================================================
            # LOAD DATA
            # =======================================================

            with st.spinner("Loading dataset..."):

                if uploaded_file is not None:

                    suffix = Path(
                        uploaded_file.name
                    ).suffix

                    with tempfile.NamedTemporaryFile(
                        delete=False,
                        suffix=suffix,
                    ) as temp:

                        temp.write(
                            uploaded_file.getbuffer()
                        )

                        temp_path = temp.name

                    adata = load_scrna_file(temp_path)

                elif (
                    matrix_file is not None
                    and barcodes_file is not None
                    and features_file is not None
                ):

                    with tempfile.TemporaryDirectory() as temp_dir:

                        temp_dir = Path(temp_dir)

                        matrix_path = (
                            temp_dir / matrix_file.name
                        )

                        barcodes_path = (
                            temp_dir / barcodes_file.name
                        )

                        features_path = (
                            temp_dir / features_file.name
                        )

                        matrix_path.write_bytes(
                            matrix_file.getbuffer()
                        )

                        barcodes_path.write_bytes(
                            barcodes_file.getbuffer()
                        )

                        features_path.write_bytes(
                            features_file.getbuffer()
                        )

                        adata = load_10x_mtx(
                            temp_dir
                        )

                else:

                    st.info(
                        "No dataset uploaded — using the built-in "
                        "PBMC3k dataset for demonstration."
                    )

                    adata = load_pbmc3k()

            # =======================================================
            # VALIDATE
            # =======================================================

            validation = validate_scrna(adata)

            st.success(
                "Dataset loaded successfully."
            )

            # =======================================================
            # DATASET SUMMARY
            # =======================================================

            st.subheader("Dataset Summary")

            c1, c2, c3 = st.columns(3)

            c1.metric(
                "Cells",
                validation["cells"],
            )

            c2.metric(
                "Genes",
                validation["genes"],
            )

            c3.metric(
                "Layers",
                len(validation["layers"]),
            )

            # =======================================================
            # QC
            # =======================================================

            with st.spinner(
                "Calculating quality-control metrics..."
            ):

                adata = calculate_basic_qc(
                    adata
                )

                adata = flag_low_quality_cells(
                    adata,
                    min_genes=min_genes,
                    max_mito_percent=max_mito,
                )

            # =======================================================
            # QC VISUALIZATION
            # =======================================================

            st.subheader(
                "Quality-Control Overview"
            )

            try:

                qc_fig = create_qc_figure(
                    adata
                )

                st.pyplot(
                    qc_fig,
                    use_container_width=True,
                )

            except Exception as error:

                st.warning(
                    f"QC visualization unavailable: {error}"
                )

            # =======================================================
            # ARTIFACT INVESTIGATION
            # =======================================================

            with st.spinner(
                "Investigating technical artifacts..."
            ):

                adata = detect_doublets(
                    adata,
                    expected_doublet_rate=doublet_rate,
                )

                adata = calculate_ambient_rna_evidence(
                    adata
                )

                adata = calculate_stress_score(
                    adata
                )

            # =======================================================
            # BIOLOGICAL STRUCTURE
            # =======================================================

            with st.spinner(
                "Analyzing cellular structure..."
            ):

                processed = preprocess_scrna(
                    adata
                )

                processed = cluster_cells(
                    processed
                )

                processed = annotate_clusters(
                    processed
                )

            # -------------------------------------------------------
            # Carry QC/artifact information onto processed object
            # -------------------------------------------------------

            for column in adata.obs.columns:

                if column not in processed.obs:

                    processed.obs[column] = (
                        adata.obs.loc[
                            processed.obs_names,
                            column,
                        ]
                    )

            processed.uns.update(
                adata.uns
            )

            # =======================================================
            # AGGREGATE EVIDENCE
            # =======================================================

            processed = aggregate_artifact_evidence(
                processed
            )

            reports = generate_debug_report(
                processed
            )

            summary = build_scrna_summary(
                processed,
                reports,
            )

            # =======================================================
            # DEBUGGING OVERVIEW
            # =======================================================

            st.subheader(
                "Debugging Overview"
            )

            c1, c2, c3, c4 = st.columns(4)

            c1.metric(
                "Cells analyzed",
                summary["cells_analyzed"],
            )

            c2.metric(
                "Cells requiring investigation",
                summary[
                    "cells_requiring_investigation"
                ],
            )

            c3.metric(
                "Clusters",
                summary["clusters"],
            )

            doublet_status = processed.uns.get(
                "doublet_detection_status",
                "unknown",
            )

            c4.metric(
                "Doublet detection",
                (
                    "Completed"
                    if doublet_status == "completed"
                    else "Limited"
                ),
            )

            # =======================================================
            # ARTIFACT BREAKDOWN
            # =======================================================

            st.subheader(
                "Artifact Evidence Breakdown"
            )

            low_quality_count = int(
                processed.obs[
                    "low_quality_flag"
                ].sum()
            )

            doublet_count = int(
                processed.obs[
                    "predicted_doublet"
                ].sum()
            )

            stress_count = int(
                processed.obs[
                    "stress_flag"
                ].sum()
            )

            ambient_count = int(
                processed.obs[
                    "ambient_evidence"
                ].sum()
            )

            artifact_df = pd.DataFrame(
                {
                    "Evidence type": [
                        "Low quality",
                        "Potential doublet",
                        "High stress",
                        "Unexpected lineage-marker signal",
                    ],
                    "Cells": [
                        low_quality_count,
                        doublet_count,
                        stress_count,
                        ambient_count,
                    ],
                }
            )

            st.dataframe(
                artifact_df,
                use_container_width=True,
                hide_index=True,
            )

            # =======================================================
            # UMAP
            # =======================================================

            st.subheader(
                "Detected Cellular Structure"
            )

            try:

                umap_fig = create_umap_figure(
                    processed,
                    color="leiden",
                )

                st.pyplot(
                    umap_fig,
                    use_container_width=True,
                )

            except Exception as error:

                st.warning(
                    f"UMAP visualization unavailable: {error}"
                )

            # =======================================================
            # CELL TYPE STRUCTURE
            # =======================================================

            cell_types = processed.uns.get(
                "cluster_cell_types",
                {},
            )

            if cell_types:

                st.subheader(
                    "Detected Biological Structure"
                )

                type_df = pd.DataFrame(
                    [
                        {
                            "Cluster": cluster,
                            "Likely cell type": cell_type,
                        }
                        for cluster, cell_type
                        in cell_types.items()
                    ]
                )

                st.dataframe(
                    type_df,
                    use_container_width=True,
                    hide_index=True,
                )

            # =======================================================
            # INVESTIGATION REPORT
            # =======================================================

            st.subheader(
                "Cells Requiring Investigation"
            )

            if reports:

                report_df = reports_to_dataframe(
                    reports
                )

                st.dataframe(
                    report_df,
                    use_container_width=True,
                    hide_index=True,
                )

                # Downloadable report

                csv_data = report_df.to_csv(
                    index=False
                )

                st.download_button(
                    "Download investigation report",
                    data=csv_data,
                    file_name=(
                        "scrna_debug_report.csv"
                    ),
                    mime="text/csv",
                )

            else:

                st.success(
                    "No cells currently meet the "
                    "multi-evidence investigation criterion."
                )

            # =======================================================
            # OPTIONAL CORRECTION
            # =======================================================

            st.subheader(
                "Optional Correction & Re-analysis"
            )

            if reports:

                st.write(
                    "Cells meeting the investigation criterion "
                    "are excluded from a separate corrected copy. "
                    "The original dataset remains unchanged."
                )

                run_correction = st.button(
                    "Apply correction and re-analyze",
                    use_container_width=True,
                )

                if run_correction:

                    try:

                        with st.spinner(
                            "Applying correction and re-analyzing..."
                        ):

                            corrected = (
                                correct_suspicious_cells(
                                    processed
                                )
                            )

                            corrected = reanalyze(
                                corrected
                            )

                            comparison = (
                                compare_before_after(
                                    processed,
                                    corrected,
                                )
                            )

                            score = (
                                calculate_trust_score(
                                    before_cells=(
                                        processed.n_obs
                                    ),
                                    after_cells=(
                                        corrected.n_obs
                                    ),
                                    before_clusters=(
                                        processed.obs[
                                            "leiden"
                                        ].nunique()
                                    ),
                                    after_clusters=(
                                        corrected.obs[
                                            "leiden"
                                        ].nunique()
                                    ),
                                    flagged_cells=len(
                                        reports
                                    ),
                                )
                            )

                            interpretation = (
                                interpret_trust_score(
                                    score
                                )
                            )

                        # ------------------------------------------------
                        # Before / after metrics
                        # ------------------------------------------------

                        c1, c2, c3 = st.columns(3)

                        c1.metric(
                            "Cells before",
                            comparison[
                                "before_cells"
                            ],
                        )

                        c2.metric(
                            "Cells after",
                            comparison[
                                "after_cells"
                            ],
                        )

                        c3.metric(
                            "Clusters after",
                            comparison[
                                "after_clusters"
                            ],
                        )

                        st.info(
                            f"**Debugging Stability Score: "
                            f"{score}/100**\n\n"
                            f"{interpretation}"
                        )

                        before_after = pd.DataFrame(
                            {
                                "Metric": [
                                    "Cells",
                                    "Clusters",
                                ],
                                "Before": [
                                    comparison[
                                        "before_cells"
                                    ],
                                    comparison[
                                        "before_clusters"
                                    ],
                                ],
                                "After": [
                                    comparison[
                                        "after_cells"
                                    ],
                                    comparison[
                                        "after_clusters"
                                    ],
                                ],
                            }
                        )

                        st.dataframe(
                            before_after,
                            use_container_width=True,
                            hide_index=True,
                        )

                        # ------------------------------------------------
                        # Interpretation
                        # ------------------------------------------------

                        if (
                            comparison[
                                "before_clusters"
                            ]
                            == comparison[
                                "after_clusters"
                            ]
                        ):

                            st.success(
                                "The number of detected clusters "
                                "was preserved after debugging."
                            )

                        else:

                            st.warning(
                                "The detected cluster structure "
                                "changed after debugging and should "
                                "be reviewed."
                            )

                        # ------------------------------------------------
                        # Corrected UMAP
                        # ------------------------------------------------

                        st.subheader(
                            "Post-debugging Cellular Structure"
                        )

                        try:

                            corrected_umap = (
                                create_umap_figure(
                                    corrected,
                                    color="leiden",
                                )
                            )

                            st.pyplot(
                                corrected_umap,
                                use_container_width=True,
                            )

                        except Exception as error:

                            st.warning(
                                "Post-debugging UMAP "
                                f"unavailable: {error}"
                            )

                    except Exception as error:

                        st.error(
                            "Correction could not be completed."
                        )

                        st.warning(
                            f"Reason: {error}"
                        )

                        st.info(
                            "The original dataset was preserved "
                            "and was not overwritten."
                        )

            else:

                st.info(
                    "No correction was required because no cells "
                    "met the investigation criterion."
                )

            # =======================================================
            # FINAL INTERPRETATION
            # =======================================================

            st.subheader(
                "Biological Interpretation"
            )

            st.write(
                "The debugger identifies evidence that may indicate "
                "technical artifacts affecting cellular structure. "
                "A flagged cell is not automatically considered false."
            )

            st.write(
                "The biological conclusion should be considered "
                "more stable when the major cellular structure "
                "remains consistent after debugging."
            )

            # =======================================================
            # DISCLAIMER
            # =======================================================

            st.markdown("---")

            st.caption(
                "Biomedical Data Debugger flags evidence requiring "
                "investigation. It does not establish biological truth "
                "or clinical validity, and it does not automatically "
                "declare cells false."
            )


# ===================================================================
# DNA VARIANT DEBUGGER
# ===================================================================

else:

    st.header(
        "DNA Variant Debugger"
    )

    st.write(
        "Evaluate available sequencing evidence and identify "
        "variants requiring investigation."
    )

    # ---------------------------------------------------------------
    # Upload
    # ---------------------------------------------------------------

    uploaded_file = st.file_uploader(
        "Upload a VCF or VCF.GZ file",
        type=["vcf", "gz"],
    )

    if uploaded_file is not None:

        run_dna = st.button(
            "Run DNA Variant Debugger",
            type="primary",
            use_container_width=True,
        )

        if run_dna:

            temp_path = (
                PROJECT_ROOT
                / "data"
                / "raw"
                / "uploaded_variants.vcf"
            )

            temp_path.parent.mkdir(
                parents=True,
                exist_ok=True,
            )

            temp_path.write_bytes(
                uploaded_file.getbuffer()
            )

            try:

                # ===================================================
                # LOAD
                # ===================================================

                records = load_vcf(
                    temp_path
                )

                # ===================================================
                # DEBUG
                # ===================================================

                reports, summary = (
                    debug_variants(records)
                )

                st.success(
                    "DNA variant debugging complete."
                )

                # ===================================================
                # SUMMARY
                # ===================================================

                c1, c2, c3 = st.columns(3)

                c1.metric(
                    "Variants analyzed",
                    summary[
                        "variants_analyzed"
                    ],
                )

                c2.metric(
                    "Requires investigation",
                    summary[
                        "suspicious_records"
                    ],
                )

                c3.metric(
                    "Samples detected",
                    len(
                        summary["samples"]
                    ),
                )

                # ===================================================
                # VARIANT TABLE
                # ===================================================

                st.subheader(
                    "Variant Evidence"
                )

                display_rows = []

                for report in reports:

                    display_rows.append(
                        {
                            "Variant": (
                                report["variant"]
                            ),
                            "Sample": (
                                report["sample"]
                            ),
                            "Genotype": (
                                report["genotype"]
                            ),
                            "Depth": (
                                report["depth"]
                            ),
                            "Alt reads": (
                                report["alt_reads"]
                            ),
                            "Allele fraction": (
                                report[
                                    "allele_fraction"
                                ]
                            ),
                            "Mapping quality": (
                                report[
                                    "mapping_quality"
                                ]
                            ),
                            "Variant quality": (
                                report["quality"]
                            ),
                            "Evidence": (
                                report[
                                    "artifact_evidence"
                                ]
                            ),
                            "Status": (
                                "Requires investigation"
                                if report[
                                    "requires_investigation"
                                ]
                                else
                                "No major artifact evidence"
                            ),
                        }
                    )

                if display_rows:

                    variant_df = pd.DataFrame(
                        display_rows
                    )

                    st.dataframe(
                        variant_df,
                        use_container_width=True,
                        hide_index=True,
                    )

                    st.download_button(
                        "Download variant report",
                        data=variant_df.to_csv(
                            index=False
                        ),
                        file_name=(
                            "variant_debug_report.csv"
                        ),
                        mime="text/csv",
                    )

                # ===================================================
                # FLAGGED VARIANTS
                # ===================================================

                suspicious = [
                    report
                    for report in reports
                    if report[
                        "requires_investigation"
                    ]
                ]

                if suspicious:

                    st.subheader(
                        "Why These Variants Were Flagged"
                    )

                    for report in suspicious:

                        with st.expander(
                            report["variant"]
                        ):

                            for reason in report[
                                "reasons"
                            ]:

                                st.write(
                                    f"- {reason}"
                                )

                            st.caption(
                                "Evidence requiring investigation "
                                "does not by itself prove that "
                                "the variant is false."
                            )

                else:

                    st.success(
                        "No major artifact evidence was detected "
                        "under the current thresholds."
                    )

            except Exception as error:

                st.error(
                    "The VCF could not be analyzed."
                )

                st.warning(
                    f"Reason: {error}"
                )

    else:

        st.info(
            "Upload a VCF file to begin."
        )


# -------------------------------------------------------------------
# Footer
# -------------------------------------------------------------------

st.markdown("---")

st.caption(
    "Biomedical Data Debugger • Evidence-driven biomedical data quality analysis"
)