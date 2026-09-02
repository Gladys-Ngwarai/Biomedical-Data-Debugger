CELL_TYPE_MARKERS = {
    "T_cell": ["CD3D", "CD3E", "IL7R", "LTB"],
    "NK_cell": ["NKG7", "GNLY", "CCL5"],
    "B_cell": ["MS4A1", "CD79A", "CD74", "HLA-DRA"],
    "Monocyte": ["LYZ", "S100A8", "S100A9", "TYROBP"],
    "FCGR3A_Monocyte": ["FCGR3A", "MS4A7", "LST1"],
    "Dendritic": ["FCER1A", "CST3"],
    "Platelet": ["PF4", "PPBP"],
}


def annotate_clusters(adata):
    """
    Assign a likely cell type to each cluster using marker overlap.
    """

    cluster_types = {}

    for cluster in adata.obs["leiden"].unique():

        cluster_cells = adata[adata.obs["leiden"] == cluster]

        scores = {}

        for cell_type, markers in CELL_TYPE_MARKERS.items():
            present = [
                gene for gene in markers
                if gene in adata.var_names
            ]

            if not present:
                scores[cell_type] = 0
                continue

            expression = cluster_cells[:, present].X

            if hasattr(expression, "mean"):
                score = expression.mean()

            scores[cell_type] = float(score)

        best_type = max(scores, key=scores.get)
        cluster_types[cluster] = best_type

    adata.uns["cluster_cell_types"] = cluster_types

    return adata