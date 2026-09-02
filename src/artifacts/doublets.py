import numpy as np
import scrublet as scr


def detect_doublets(
    adata,
    expected_doublet_rate=0.05,
):
    """
    Detect potential doublets using Scrublet.

    Cells are flagged, not automatically removed.

    If Scrublet cannot run on a particular dataset,
    the module records that limitation instead of
    crashing the complete analysis.
    """

    try:

        scrub = scr.Scrublet(
            adata.X,
            expected_doublet_rate=(
                expected_doublet_rate
            ),
        )

        scores, predicted = (
            scrub.scrub_doublets(
                verbose=False
            )
        )

        adata.obs[
            "doublet_score"
        ] = scores

        adata.obs[
            "predicted_doublet"
        ] = predicted

        adata.uns[
            "doublet_detection_status"
        ] = "completed"

    except Exception as error:

        adata.obs[
            "doublet_score"
        ] = np.nan

        adata.obs[
            "predicted_doublet"
        ] = False

        adata.uns[
            "doublet_detection_status"
        ] = (
            "unavailable: "
            + str(error)
        )

    return adata