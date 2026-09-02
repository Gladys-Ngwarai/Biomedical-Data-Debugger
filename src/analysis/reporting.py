import pandas as pd


def build_scrna_report(
    adata,
    before_cells,
    after_cells,
    before_clusters,
    after_clusters,
    flagged_cells,
    trust_score,
    interpretation,
):
    """
    Build a human-readable scRNA-seq debugging report.
    """

    rows = []

    suspicious = adata.obs[
        adata.obs["requires_investigation"]
    ]

    for cell_id, row in suspicious.iterrows():

        reasons = []

        if row["low_quality_flag"]:
            reasons.append("Low quality")

        if row["predicted_doublet"]:
            reasons.append("Potential doublet")

        if row["stress_flag"]:
            reasons.append("High stress response")

        if row["ambient_evidence"]:
            source = row.get(
                "ambient_source",
                ""
            )

            if source:
                reasons.append(
                    f"Unexpected {source} marker signal"
                )
            else:
                reasons.append(
                    "Lineage-marker evidence"
                )

        cluster = row.get(
            "leiden",
            "Unknown"
        )

        cell_type = adata.uns.get(
            "cluster_cell_types",
            {}
        ).get(
            str(cluster),
            "Unknown"
        )

        rows.append(
            {
                "Cell": cell_id,
                "Cluster": cluster,
                "Likely Cell Type": cell_type,
                "Evidence Count": int(
                    row[
                        "artifact_evidence_count"
                    ]
                ),
                "Reasons": "; ".join(reasons),
                "Recommendation": (
                    "Review before exclusion"
                ),
            }
        )

    report_df = pd.DataFrame(rows)

    summary = {
        "Cells Before": before_cells,
        "Cells After": after_cells,
        "Cells Removed": (
            before_cells - after_cells
        ),
        "Clusters Before": before_clusters,
        "Clusters After": after_clusters,
        "Suspicious Cells": flagged_cells,
        "Trust Score": trust_score,
        "Interpretation": interpretation,
    }

    return summary, report_df