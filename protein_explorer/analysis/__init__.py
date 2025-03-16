"""
Analysis modules for protein structures and networks.
"""

from protein_explorer.analysis.structure import (
    calculate_pca,
    compute_structural_similarity,
    calculate_distance_matrix,
    extract_coordinates
)

from protein_explorer.analysis.networks import (
    compute_eigenvector_centrality,
    perform_spectral_clustering,
    identify_key_proteins,
    calculate_network_metrics
)