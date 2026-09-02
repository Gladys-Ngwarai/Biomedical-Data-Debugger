from pathlib import Path
import scanpy as sc


SUPPORTED_SINGLE_CELL_FORMATS = {
    ".h5ad",
    ".h5",
    ".mtx",
    ".gz",
}


def load_pbmc3k():
    """
    Load the PBMC3k demonstration dataset.
    """

    return sc.datasets.pbmc3k()


def validate_scrna(adata):
    """
    Validate a single-cell AnnData object before analysis.

    Returns
    -------
    dict
        Dataset validation information.
    """

    if adata is None:
        raise ValueError(
            "No single-cell dataset was loaded."
        )

    if adata.n_obs < 10:
        raise ValueError(
            "The dataset contains fewer than 10 cells."
        )

    if adata.n_vars < 20:
        raise ValueError(
            "The dataset contains fewer than 20 genes."
        )

    if adata.obs_names is None:
        raise ValueError(
            "Cell identifiers are missing."
        )

    if adata.var_names is None:
        raise ValueError(
            "Gene identifiers are missing."
        )

    if adata.obs_names.has_duplicates:
        adata.obs_names_make_unique()

    if adata.var_names.has_duplicates:
        adata.var_names_make_unique()

    return {
        "cells": adata.n_obs,
        "genes": adata.n_vars,
        "layers": list(adata.layers.keys()),
        "obsm": list(adata.obsm.keys()),
        "uns": list(adata.uns.keys()),
    }


def load_h5ad(path):
    """
    Load an AnnData .h5ad file.
    """

    return sc.read_h5ad(path)


def load_10x_h5(path):
    """
    Load a 10x Genomics HDF5 file.
    """

    return sc.read_10x_h5(
        path
    )


def load_10x_mtx(path):
    """
    Load a 10x Matrix Market directory.

    Expected files include:
        matrix.mtx / matrix.mtx.gz
        barcodes.tsv / barcodes.tsv.gz
        features.tsv / features.tsv.gz

    Older datasets may use genes.tsv.
    """

    path = Path(path)

    return sc.read_10x_mtx(
        path,
        var_names="gene_symbols",
        make_unique=True,
    )


def load_scrna_file(path):
    """
    Automatically load a supported single-cell file.
    """

    path = Path(path)

    suffix = path.suffix.lower()

    if suffix == ".h5ad":

        adata = load_h5ad(
            path
        )

    elif suffix == ".h5":

        adata = load_10x_h5(
            path
        )

    else:

        raise ValueError(
            "Unsupported single-cell file. "
            "Use .h5ad or 10x .h5."
        )

    validate_scrna(
        adata
    )

    return adata


def prepare_counts(adata):
    """
    Prepare the dataset for artifact detection.

    Uses raw counts when available.

    The original AnnData object is not modified.
    """

    adata = adata.copy()

    # Prefer raw counts if supplied.
    if adata.raw is not None:

        try:
            raw = adata.raw.to_adata()

            # Only use raw if dimensions are sensible.
            if (
                raw.n_obs == adata.n_obs
                and raw.n_vars >= 20
            ):
                adata = raw

        except Exception:
            pass

    return adata