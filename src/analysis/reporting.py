import pandas as pd


def build_scrna_summary(adata, reports):
    """Build a concise scRNA-seq debugging summary."""

    summary = {
        "cells_analyzed": int(adata.n_obs),
        "genes_analyzed": int(adata.n_vars),
        "cells_requiring_investigation": len(reports),
        "clusters": (
            int(adata.obs["leiden"].nunique())
            if "leiden" in adata.obs
            else 0
        ),
    }

    if "low_quality_flag" in adata.obs:
        summary["low_quality_cells"] = int(
            adata.obs["low_quality_flag"].sum()
        )

    if "predicted_doublet" in adata.obs:
        summary["potential_doublets"] = int(
            adata.obs["predicted_doublet"].sum()
        )

    if "stress_flag" in adata.obs:
        summary["high_stress_cells"] = int(
            adata.obs["stress_flag"].sum()
        )

    if "ambient_evidence" in adata.obs:
        summary["ambient_marker_evidence"] = int(
            adata.obs["ambient_evidence"].sum()
        )

    return summary


def build_scrna_report(adata, reports):
    """
    Create a human-readable debugging report.

    This report describes evidence rather than claiming
    biological truth.
    """

    summary = build_scrna_summary(
        adata,
        reports,
    )

    return {
        "summary": summary,
        "investigations": reports,
        "conclusion": (
            "Cells listed as requiring investigation "
            "have multiple artifact-related signals. "
            "These signals do not prove that a cell "
            "is biologically false."
        ),
    }


def reports_to_dataframe(reports):
    """Convert investigation reports into a UI-friendly table."""

    if not reports:
        return pd.DataFrame(
            columns=[
                "Cell",
                "Cluster",
                "Evidence count",
                "Reasons",
                "Recommendation",
            ]
        )

    rows = []

    for report in reports:

        rows.append(
            {
                "Cell": report["cell"],
                "Cluster": report["cluster"],
                "Evidence count": report[
                    "artifact_evidence_count"
                ],
                "Reasons": "; ".join(
                    report["reasons"]
                ),
                "Recommendation": report[
                    "recommendation"
                ],
            }
        )

    return pd.DataFrame(rows)