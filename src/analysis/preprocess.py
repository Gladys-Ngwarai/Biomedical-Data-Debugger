import numpy as np
import scanpy as sc


def preprocess_scrna(
    adata,
    n_top_genes=2000,
):
    """
    Prepare an scRNA-seq dataset for clustering.

    Handles sparse matrices, zero-count genes,
    invalid values, and datasets with fewer than
    n_top_genes available genes.
    """

    adata = adata.copy()

    # Remove cells with no counts.
    sc.pp.filter_cells(
        adata,
        min_counts=1,
    )

    # Remove genes with no counts.
    sc.pp.filter_genes(
        adata,
        min_cells=1,
    )

    if adata.n_obs < 10:
        raise ValueError(
            "Too few cells remain after removing "
            "empty cells."
        )

    if adata.n_vars < 20:
        raise ValueError(
            "Too few genes remain after removing "
            "unexpressed genes."
        )

    # Clean invalid numeric values.
    if hasattr(
        adata.X,
        "data",
    ):

        invalid = ~np.isfinite(
            adata.X.data
        )

        adata.X.data[
            invalid
        ] = 0

    else:

        adata.X = np.nan_to_num(
            adata.X,
            nan=0.0,
            posinf=0.0,
            neginf=0.0,
        )

    # Normalize.
    sc.pp.normalize_total(
        adata,
        target_sum=1e4,
    )

    # Log transform.
    sc.pp.log1p(
        adata
    )

    # Number of HVGs must be smaller than
    # the number of available genes.
    n_top_genes = min(
        n_top_genes,
        adata.n_vars,
    )

    if n_top_genes < 10:

        raise ValueError(
            "Too few genes are available for "
            "highly-variable gene selection."
        )

    try:

        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor="seurat",
            subset=False,
        )

    except (
        ValueError,
        RuntimeError,
    ):

        # Robust fallback for datasets where
        # Seurat-style binning cannot be performed.
        adata.var[
            "highly_variable"
        ] = True

    selected = adata.var[
        "highly_variable"
    ]

    if selected.sum() < 10:

        adata.var[
            "highly_variable"
        ] = True

    adata = adata[
        :,
        adata.var[
            "highly_variable"
        ],
    ].copy()

    sc.pp.scale(
        adata,
        max_value=10,
    )

    return adata