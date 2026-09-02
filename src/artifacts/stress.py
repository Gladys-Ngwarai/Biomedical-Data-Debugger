import scanpy as sc


STRESS_GENES = [
    "FOS",
    "JUN",
    "HSPA1A",
    "HSPA1B",
    "DNAJB1",
    "DUSP1",
    "ATF3",
]


def calculate_stress_score(adata):
    """
    Calculate a stress-response score for each cell.
    """

    genes_present = [
        gene for gene in STRESS_GENES
        if gene in adata.var_names
    ]

    if not genes_present:
        adata.obs["stress_score"] = 0.0
        return adata

    sc.tl.score_genes(
        adata,
        gene_list=genes_present,
        score_name="stress_score"
    )

    return adata