import scanpy as sc


def load_pbmc3k():
    """
    Load the PBMC 3k demonstration dataset.

    Returns
    -------
    AnnData
        PBMC 3k single-cell expression dataset.
    """
    adata = sc.datasets.pbmc3k()

    return adata