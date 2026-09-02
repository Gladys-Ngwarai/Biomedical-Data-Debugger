import scanpy as sc


CELL_TYPE_MARKERS = {
    "T_cell": ["CD3D", "CD3E"],
    "NK_cell": ["NKG7", "GNLY"],
    "B_cell": ["MS4A1", "CD79A"],
    "Monocyte": ["LYZ", "S100A8", "S100A9"],
    "Platelet": ["PF4", "PPBP"],
}


def calculate_ambient_rna_evidence(adata):
    """
    Calculate marker-contamination evidence.

    This is an evidence signal, NOT proof of ambient RNA.
    Strong expression of another lineage's markers can also
    arise from doublets or genuine transitional biology.
    """

    present_markers = {
        cell_type: [
            gene for gene in markers
            if gene in adata.var_names
        ]
        for cell_type, markers in CELL_TYPE_MARKERS.items()
    }

    scores = {}

    for cell_type, genes in present_markers.items():
        if genes:
            score_name = f"ambient_{cell_type}"
            sc.tl.score_genes(
                adata,
                gene_list=genes,
                score_name=score_name
            )
            scores[cell_type] = score_name

    adata.uns["ambient_marker_scores"] = scores

    return adata