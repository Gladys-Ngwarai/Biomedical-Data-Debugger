import scanpy as sc
import pandas as pd


def calculate_basic_qc(adata):
    """
    Calculate basic quality-control metrics for an AnnData object.

    Parameters
    ----------
    adata : anndata.AnnData
        Single-cell expression dataset.

    Returns
    -------
    anndata.AnnData
        AnnData object containing calculated QC metrics.
    """

    # Calculate mitochondrial gene percentage
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        inplace=True
    )

    return adata