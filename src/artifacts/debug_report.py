def generate_debug_report(adata):
    """
    Generate an interpretable report for suspicious cells.
    """

    reports = []

    for cell_id, row in adata.obs.iterrows():

        if not row["requires_investigation"]:
            continue

        reasons = []

        if row["low_quality_flag"]:
            reasons.append("low quality")

        if row["predicted_doublet"]:
            reasons.append("potential doublet")

        if row["stress_flag"]:
            reasons.append("high stress response")

        if row["ambient_evidence"]:
            reasons.append(
                f"unexpected {row['ambient_source']} marker signal"
            )

        reports.append({
            "cell": cell_id,
            "cluster": row["leiden"],
            "artifact_evidence_count": int(
                row["artifact_evidence_count"]
            ),
            "reasons": reasons,
            "recommendation": "Review before exclusion",
        })

    return reports