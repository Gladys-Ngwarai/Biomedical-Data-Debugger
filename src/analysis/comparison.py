def compare_before_after(before, after):
    """
    Compare biological structure before and after correction.
    """

    before_clusters = before.obs["leiden"].value_counts().sort_index()
    after_clusters = after.obs["leiden"].value_counts().sort_index()

    comparison = {
        "before_cells": before.n_obs,
        "after_cells": after.n_obs,
        "cells_removed": before.n_obs - after.n_obs,
        "before_clusters": before.obs["leiden"].nunique(),
        "after_clusters": after.obs["leiden"].nunique(),
        "before_cluster_sizes": before_clusters.to_dict(),
        "after_cluster_sizes": after_clusters.to_dict(),
    }

    return comparison