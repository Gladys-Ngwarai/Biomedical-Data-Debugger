import scrublet as scr


def detect_doublets(adata, expected_doublet_rate=0.05):
    """
    Detect potential doublets using Scrublet.

    Cells are flagged but NOT removed.

    Parameters
    ----------
    adata : anndata.AnnData
        Raw-count single-cell dataset.

    expected_doublet_rate : float
        Expected proportion of doublets in the dataset.

    Returns
    -------
    anndata.AnnData
        AnnData object containing doublet scores and predictions.
    """

    if adata.n_obs == 0:
        raise ValueError("The dataset contains no cells.")

    if not 0 < expected_doublet_rate < 1:
        raise ValueError(
            "expected_doublet_rate must be between 0 and 1."
        )

    scrub = scr.Scrublet(
        adata.X,
        expected_doublet_rate=expected_doublet_rate
    )

    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        verbose=False
    )

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    return adata