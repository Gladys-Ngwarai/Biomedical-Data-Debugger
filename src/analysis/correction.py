import numpy as np
import scanpy as sc


def correct_suspicious_cells(adata):
    """
    Create a corrected copy by excluding cells with
    multiple independent artifact signals.

    The original dataset is preserved.
    """

    if (
        "requires_investigation"
        not in adata.obs
    ):

        raise ValueError(
            "Run artifact evidence aggregation first."
        )

    corrected = adata[
        ~adata.obs[
            "requires_investigation"
        ]
    ].copy()

    if corrected.n_obs < 10:

        raise ValueError(
            "Correction would leave fewer than "
            "10 cells. The original dataset has "
            "been preserved."
        )

    return corrected


def _clean_matrix(adata):
    """
    Remove empty cells/genes and invalid values.
    """

    adata = adata.copy()

    sc.pp.filter_cells(
        adata,
        min_counts=1,
    )

    sc.pp.filter_genes(
        adata,
        min_cells=1,
    )

    if adata.n_obs < 10:
        raise ValueError(
            "Too few valid cells remain."
        )

    if adata.n_vars < 20:
        raise ValueError(
            "Too few valid genes remain."
        )

    if hasattr(
        adata.X,
        "data",
    ):

        adata.X.data[
            ~np.isfinite(
                adata.X.data
            )
        ] = 0

    else:

        adata.X = np.nan_to_num(
            adata.X,
            nan=0.0,
            posinf=0.0,
            neginf=0.0,
        )

    return adata


def reanalyze(
    adata,
    n_top_genes=2000,
):
    """
    Re-run the biological analysis after correction.

    Steps:
        filtering
        normalization
        log transformation
        HVG selection
        scaling
        PCA
        neighbors
        UMAP
        Leiden clustering
    """

    adata = _clean_matrix(
        adata
    )

    sc.pp.normalize_total(
        adata,
        target_sum=1e4,
    )

    sc.pp.log1p(
        adata
    )

    n_top_genes = min(
        n_top_genes,
        adata.n_vars,
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

        adata.var[
            "highly_variable"
        ] = True

    if (
        "highly_variable"
        not in adata.var
        or
        adata.var[
            "highly_variable"
        ].sum() < 10
    ):

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

    n_pcs = min(
        30,
        adata.n_obs - 1,
        adata.n_vars - 1,
    )

    if n_pcs < 2:
        raise ValueError(
            "Not enough dimensions remain "
            "for PCA."
        )

    sc.tl.pca(
        adata,
        n_comps=n_pcs,
        svd_solver="arpack",
    )

    n_neighbors = min(
        10,
        adata.n_obs - 1,
    )

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
    )

    sc.tl.umap(
        adata
    )

    sc.tl.leiden(
        adata,
        resolution=0.5,
    )

    return adata