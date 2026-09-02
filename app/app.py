import sys
from pathlib import Path
import tempfile

import pandas as pd
import streamlit as st


PROJECT_ROOT = Path(
    __file__
).resolve().parents[1]

sys.path.insert(
    0,
    str(PROJECT_ROOT),
)


from src.analysis.load_data import (
    load_scrna_file,
    load_10x_mtx,
    validate_scrna,
)

from src.qc.basic_qc import (
    calculate_basic_qc,
)

from src.artifacts.low_quality import (
    flag_low_quality_cells,
)

from src.artifacts.doublets import (
    detect_doublets,
)

from src.artifacts.ambient_rna import (
    calculate_ambient_rna_evidence,
)

from src.artifacts.stress import (
    calculate_stress_score,
)

from src.artifacts.evidence import (
    aggregate_artifact_evidence,
)

from src.artifacts.debug_report import (
    generate_debug_report,
)

from src.analysis.preprocess import (
    preprocess_scrna,
)

from src.analysis.clustering import (
    cluster_cells,
)

from src.analysis.cell_types import (
    annotate_clusters,
)

from src.analysis.correction import (
    correct_suspicious_cells,
    reanalyze,
)

from src.analysis.comparison import (
    compare_before_after,
)

from src.analysis.trust_score import (
    calculate_trust_score,
    interpret_trust_score,
)

from src.variants.load_vcf import (
    load_vcf,
)

from src.variants.variant_debugger import (
    debug_variants,
)


st.set_page_config(
    page_title="Biomedical Data Debugger",
    page_icon="🧬",
    layout="wide",
)


st.title(
    "🧬 Biomedical Data Debugger"
)

st.caption(
    "Don't trust the biological conclusion "
    "until you debug the data."
)


module = st.sidebar.selectbox(
    "Analysis module",
    [
        "Single-cell RNA-seq",
        "DNA Variant Debugger",
    ],
)


# ============================================================
# SINGLE-CELL RNA-SEQ
# ============================================================

if module == "Single-cell RNA-seq":

    st.header(
        "Single-cell RNA-seq Debugger"
    )

    st.write(
        "Investigate whether apparent cellular "
        "structure is biologically coherent or "
        "potentially driven by technical artifacts."
    )

    st.subheader(
        "1. Dataset"
    )

    uploaded_file = st.file_uploader(
        "Upload an scRNA-seq dataset",
        type=[
            "h5ad",
            "h5",
        ],
        help=(
            "Supported: AnnData .h5ad and "
            "10x Genomics .h5."
        ),
    )

    st.caption(
        "10x Matrix Market datasets can be added "
        "in the advanced input section below."
    )

    # --------------------------------------------------------
    # Advanced Matrix Market upload
    # --------------------------------------------------------

    with st.expander(
        "Advanced: 10x Matrix Market input"
    ):

        matrix_file = st.file_uploader(
            "matrix.mtx or matrix.mtx.gz",
            type=[
                "mtx",
                "gz",
            ],
            key="matrix_file",
        )

        barcodes_file = st.file_uploader(
            "barcodes.tsv or barcodes.tsv.gz",
            type=[
                "tsv",
                "gz",
            ],
            key="barcodes_file",
        )

        features_file = st.file_uploader(
            "features.tsv, genes.tsv, or compressed equivalent",
            type=[
                "tsv",
                "gz",
            ],
            key="features_file",
        )

    st.subheader(
        "2. Debugging settings"
    )

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

    run_debugger = st.button(
        "Run scRNA-seq Debugger",
        type="primary",
        use_container_width=True,
    )

    if run_debugger:

        try:

            # ==================================================
            # LOAD
            # ==================================================

            with st.spinner(
                "Loading dataset..."
            ):

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

                    adata = load_scrna_file(
                        temp_path
                    )

                elif (
                    matrix_file is not None
                    and barcodes_file is not None
                    and features_file is not None
                ):

                    with tempfile.TemporaryDirectory() as temp_dir:

                        temp_dir = Path(
                            temp_dir
                        )

                        matrix_path = (
                            temp_dir
                            / matrix_file.name
                        )

                        barcodes_path = (
                            temp_dir
                            / barcodes_file.name
                        )

                        features_path = (
                            temp_dir
                            / features_file.name
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

                    # Demo fallback
                    from src.analysis.load_data import (
                        load_pbmc3k
                    )

                    adata = load_pbmc3k()

            validation = validate_scrna(
                adata
            )

            st.success(
                "Dataset loaded successfully."
            )

            # ==================================================
            # DATASET SUMMARY
            # ==================================================

            st.subheader(
                "Dataset Summary"
            )

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
                len(
                    validation["layers"]
                ),
            )

            # ==================================================
            # QC
            # ==================================================

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

            # ==================================================
            # ARTIFACT DETECTION
            # ==================================================

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

            # ==================================================
            # BIOLOGICAL STRUCTURE
            # ==================================================

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

            # Transfer artifact information
            # onto clustered dataset.
            for column in adata.obs.columns:

                if column not in processed.obs:

                    processed.obs[
                        column
                    ] = adata.obs.loc[
                        processed.obs_names,
                        column,
                    ]

            processed.uns.update(
                adata.uns
            )

            # ==================================================
            # EVIDENCE AGGREGATION
            # ==================================================

            processed = aggregate_artifact_evidence(
                processed
            )

            reports = generate_debug_report(
                processed
            )

            # ==================================================
            # OVERVIEW
            # ==================================================

            st.subheader(
                "Debugging Overview"
            )

            c1, c2, c3, c4 = st.columns(4)

            c1.metric(
                "Cells analyzed",
                processed.n_obs,
            )

            c2.metric(
                "Cells requiring investigation",
                len(reports),
            )

            c3.metric(
                "Clusters",
                processed.obs[
                    "leiden"
                ].nunique(),
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

            # ==================================================
            # ARTIFACT BREAKDOWN
            # ==================================================

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

            # ==================================================
            # SUSPICIOUS CELLS
            # ==================================================

            st.subheader(
                "Cells Requiring Investigation"
            )

            if reports:

                report_df = pd.DataFrame(
                    reports
                )

                st.dataframe(
                    report_df,
                    use_container_width=True,
                    hide_index=True,
                )

            else:

                st.success(
                    "No cells currently meet the "
                    "multi-evidence investigation criterion."
                )

            # ==================================================
            # CELL TYPES
            # ==================================================

            st.subheader(
                "Detected Biological Structure"
            )

            cell_types = processed.uns.get(
                "cluster_cell_types",
                {},
            )

            if cell_types:

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

            # ==================================================
            # CORRECTION
            # ==================================================

            st.subheader(
                "Correction & Re-analysis"
            )

            if reports:

                with st.spinner(
                    "Applying debugging correction and re-analyzing..."
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

                    score = calculate_trust_score(
                        before_cells=processed.n_obs,
                        after_cells=corrected.n_obs,
                        before_clusters=processed.obs[
                            "leiden"
                        ].nunique(),
                        after_clusters=corrected.obs[
                            "leiden"
                        ].nunique(),
                        flagged_cells=len(
                            reports
                        ),
                    )

                    interpretation = (
                        interpret_trust_score(
                            score
                        )
                    )

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

            else:

                st.info(
                    "No correction was required because "
                    "no cells met the investigation criterion."
                )

            # ==================================================
            # SCIENTIFIC LIMITATION
            # ==================================================

            st.caption(
                "Biomedical Data Debugger flags evidence "
                "requiring investigation. It does not establish "
                "biological truth or clinical validity, and it "
                "does not automatically declare cells or variants false."
            )

        except Exception as error:

            st.error(
                "The analysis could not be completed."
            )

            st.warning(
                f"Reason: {error}"
            )

            st.info(
                "Your original uploaded dataset was not "
                "modified or overwritten."
            )


# ============================================================
# DNA VARIANT DEBUGGER
# ============================================================

else:

    st.header(
        "DNA Variant Debugger"
    )

    st.write(
        "Evaluate available sequencing evidence "
        "and identify variants requiring investigation."
    )

    uploaded_file = st.file_uploader(
        "Upload a VCF or VCF.GZ file",
        type=[
            "vcf",
            "gz",
        ],
    )

    if uploaded_file is not None:

        if st.button(
            "Run DNA Variant Debugger",
            type="primary",
            use_container_width=True,
        ):

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

                records = load_vcf(
                    temp_path
                )

                reports, summary = (
                    debug_variants(
                        records
                    )
                )

                st.success(
                    "DNA variant debugging complete."
                )

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
                        summary[
                            "samples"
                        ]
                    ),
                )

                st.subheader(
                    "Variant Evidence"
                )

                display_rows = []

                for report in reports:

                    display_rows.append(
                        {
                            "Variant":
                                report[
                                    "variant"
                                ],
                            "Sample":
                                report[
                                    "sample"
                                ],
                            "Genotype":
                                report[
                                    "genotype"
                                ],
                            "Depth":
                                report[
                                    "depth"
                                ],
                            "Alt reads":
                                report[
                                    "alt_reads"
                                ],
                            "Allele fraction":
                                report[
                                    "allele_fraction"
                                ],
                            "Mapping quality":
                                report[
                                    "mapping_quality"
                                ],
                            "Variant quality":
                                report[
                                    "quality"
                                ],
                            "Evidence":
                                report[
                                    "artifact_evidence"
                                ],
                            "Status":
                                (
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

                    st.dataframe(
                        pd.DataFrame(
                            display_rows
                        ),
                        use_container_width=True,
                        hide_index=True,
                    )

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
                            report[
                                "variant"
                            ]
                        ):

                            for reason in report[
                                "reasons"
                            ]:

                                st.write(
                                    f"- {reason}"
                                )

                            st.caption(
                                "Evidence requiring investigation "
                                "does not by itself prove that the "
                                "variant is false."
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