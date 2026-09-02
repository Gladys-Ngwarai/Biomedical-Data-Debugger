def calculate_trust_score(
    before_cells,
    after_cells,
    before_clusters,
    after_clusters,
    flagged_cells,
):
    """
    Calculate a simple interpretable trust score.

    This is a prototype score for the Biomedical Data Debugger.
    It is not a validated statistical confidence measure.
    """

    if before_cells == 0:
        return 0.0

    removal_fraction = flagged_cells / before_cells

    # Start with full trust.
    score = 100.0

    # Penalize excessive removal.
    score -= min(removal_fraction * 100, 20)

    # Penalize major structural changes.
    cluster_change = abs(
        before_clusters - after_clusters
    )

    score -= cluster_change * 5

    # Keep within 0–100.
    score = max(0.0, min(100.0, score))

    return round(score, 1)


def interpret_trust_score(score):
    """
    Convert the numerical score into an interpretable conclusion.
    """

    if score >= 85:
        return "Biological signal appears relatively stable after debugging."

    if score >= 70:
        return "Biological signal is mostly preserved, but should be reviewed."

    if score >= 50:
        return "Biological structure changed substantially and requires review."

    return "Biological conclusion may be strongly affected by technical artifacts."