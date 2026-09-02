def aggregate_artifact_evidence(adata):
    """
    Combine artifact signals into cell-level evidence.

    Strong artifact signals:
    - low quality
    - predicted doublet
    - high stress

    Ambient/lineage-marker signal is recorded separately
    because marker expression alone cannot prove ambient RNA.

    Cells are flagged for investigation rather than automatically
    declared biologically false.
    """

    # -------------------------
    # STRESS
    # -------------------------

    stress_threshold = adata.obs["stress_score"].quantile(0.90)

    adata.obs["stress_flag"] = (
        adata.obs["stress_score"] >= stress_threshold
    )

    # -------------------------
    # AMBIENT / LINEAGE EVIDENCE
    # -------------------------

    adata.obs["ambient_evidence"] = False
    adata.obs["ambient_source"] = ""

    ambient_scores = adata.uns.get(
        "ambient_marker_scores",
        {}
    )

    if ambient_scores:

        for cell_id in adata.obs.index:

            strongest_lineage = None
            strongest_score = None

            for lineage, score_name in ambient_scores.items():

                if score_name not in adata.obs.columns:
                    continue

                score = float(
                    adata.obs.loc[cell_id, score_name]
                )

                if strongest_score is None or score > strongest_score:
                    strongest_score = score
                    strongest_lineage = lineage

            # Record strong marker signal for investigation.
            # It is NOT counted as an automatic removal criterion.
            if (
                strongest_lineage is not None
                and strongest_score is not None
                and strongest_score > 1.0
            ):
                adata.obs.loc[
                    cell_id,
                    "ambient_evidence"
                ] = True

                adata.obs.loc[
                    cell_id,
                    "ambient_source"
                ] = strongest_lineage

    # -------------------------
    # STRONG ARTIFACT EVIDENCE
    # -------------------------

    adata.obs["artifact_evidence_count"] = (
        adata.obs["low_quality_flag"].astype(int)
        + adata.obs["predicted_doublet"].astype(int)
        + adata.obs["stress_flag"].astype(int)
    )

    # Ambient evidence does NOT contribute to automatic removal.
    #
    # This prevents ordinary biological marker expression from
    # causing large numbers of cells to be removed.

    adata.obs["requires_investigation"] = (
        adata.obs["artifact_evidence_count"] >= 2
    )

    return adata