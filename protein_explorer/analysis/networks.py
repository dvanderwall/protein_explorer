"""
Functions for analyzing protein interaction networks using linear algebra.
"""

import networkx as nx
import numpy as np
from scipy import sparse
from sklearn.cluster import KMeans
import logging
from typing import Dict, List, Tuple, Optional

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def compute_eigenvector_centrality(network: nx.Graph, 
                                 weight: Optional[str] = 'confidence') -> Dict[str, float]:
    """
    Calculate eigenvector centrality for proteins in the network.
    
    Args:
        network: NetworkX Graph representing the interaction network
        weight: Edge attribute to use as weight (None for unweighted)
        
    Returns:
        Dictionary mapping protein IDs to centrality scores
    """
    try:
        centrality = nx.eigenvector_centrality(network, weight=weight)
        return centrality
    except nx.PowerIterationFailedConvergence:
        logger.warning("Eigenvector centrality failed to converge, using approximate method")
        # Fall back to approximate eigenvector centrality
        return nx.eigenvector_centrality_numpy(network, weight=weight)

def perform_spectral_clustering(network: nx.Graph, 
                              n_clusters: int = 5,
                              weight: Optional[str] = 'confidence') -> Dict[str, int]:
    """
    Perform spectral clustering on the protein interaction network.
    
    Args:
        network: NetworkX Graph representing the interaction network
        n_clusters: Number of clusters to identify
        weight: Edge attribute to use as weight (None for unweighted)
        
    Returns:
        Dictionary mapping protein IDs to cluster assignments
    """
    if network.number_of_nodes() < n_clusters:
        logger.warning(f"Number of nodes ({network.number_of_nodes()}) less than n_clusters ({n_clusters})")
        n_clusters = max(2, network.number_of_nodes() // 2)
    
    # Get adjacency matrix with weights
    if weight:
        adjacency_matrix = nx.to_numpy_array(network, weight=weight)
    else:
        adjacency_matrix = nx.to_numpy_array(network)
    
    # Compute the normalized graph Laplacian
    # L = D^(-1/2) * (D - A) * D^(-1/2)
    degree_matrix = np.diag(np.sum(adjacency_matrix, axis=1))
    degree_inv_sqrt = np.linalg.inv(np.sqrt(degree_matrix))
    laplacian = np.eye(adjacency_matrix.shape[0]) - degree_inv_sqrt @ adjacency_matrix @ degree_inv_sqrt
    
    # Compute eigenvalues and eigenvectors of the Laplacian
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
    
    # Select the eigenvectors corresponding to the k smallest non-zero eigenvalues
    # (the first eigenvalue should be close to zero)
    indices = np.argsort(eigenvalues)[1:n_clusters+1]
    features = eigenvectors[:, indices]
    
    # Normalize the rows to have unit length
    row_norms = np.sqrt(np.sum(features**2, axis=1))
    features = features / row_norms[:, np.newaxis]
    
    # Apply k-means clustering on the selected eigenvectors
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    cluster_labels = kmeans.fit_predict(features)
    
    # Map node indices to protein IDs
    protein_ids = list(network.nodes())
    clustering = {protein_ids[i]: int(cluster_labels[i]) for i in range(len(protein_ids))}
    
    return clustering

def identify_key_proteins(network: nx.Graph, 
                        metrics: List[str] = ['eigenvector', 'betweenness', 'degree'],
                        top_n: int = 10) -> Dict[str, Dict[str, float]]:
    """
    Identify key proteins in the network using various centrality metrics.
    
    Args:
        network: NetworkX Graph representing the interaction network
        metrics: List of centrality metrics to compute
        top_n: Number of top proteins to return for each metric
        
    Returns:
        Dictionary mapping metrics to dictionaries of protein IDs and scores
    """
    results = {}
    
    # Compute each requested centrality metric
    for metric in metrics:
        if metric == 'eigenvector':
            centrality = compute_eigenvector_centrality(network)
        elif metric == 'betweenness':
            centrality = nx.betweenness_centrality(network)
        elif metric == 'degree':
            centrality = nx.degree_centrality(network)
        elif metric == 'closeness':
            centrality = nx.closeness_centrality(network)
        elif metric == 'pagerank':
            centrality = nx.pagerank(network)
        else:
            logger.warning(f"Unknown centrality metric: {metric}")
            continue
        
        # Sort proteins by centrality (highest first)
        sorted_proteins = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
        
        # Keep top N proteins
        top_proteins = dict(sorted_proteins[:top_n])
        results[metric] = top_proteins
    
    return results

def calculate_network_metrics(network: nx.Graph) -> Dict[str, float]:
    """
    Calculate various network-level metrics.
    
    Args:
        network: NetworkX Graph representing the interaction network
        
    Returns:
        Dictionary of network metrics
    """
    metrics = {}
    
    # Basic network properties
    metrics['num_nodes'] = network.number_of_nodes()
    metrics['num_edges'] = network.number_of_edges()
    metrics['density'] = nx.density(network)
    
    # Check if the network is connected
    if nx.is_connected(network):
        # Connected network metrics
        metrics['average_shortest_path_length'] = nx.average_shortest_path_length(network)
        metrics['diameter'] = nx.diameter(network)
    else:
        # For disconnected networks, compute metrics on largest component
        components = list(nx.connected_components(network))
        largest_component = max(components, key=len)
        subgraph = network.subgraph(largest_component).copy()
        
        metrics['largest_component_size'] = subgraph.number_of_nodes()
        metrics['largest_component_fraction'] = subgraph.number_of_nodes() / network.number_of_nodes()
        metrics['num_components'] = len(components)
        
        if subgraph.number_of_nodes() > 1:
            metrics['average_shortest_path_length_largest_component'] = nx.average_shortest_path_length(subgraph)
            metrics['diameter_largest_component'] = nx.diameter(subgraph)
    
    # Clustering coefficient
    metrics['average_clustering'] = nx.average_clustering(network)
    
    # Degree distribution statistics
    degrees = [d for _, d in network.degree()]
    metrics['min_degree'] = min(degrees)
    metrics['max_degree'] = max(degrees)
    metrics['average_degree'] = sum(degrees) / len(degrees)
    
    return metrics

def find_modules(network: nx.Graph, resolution: float = 1.0) -> Dict[str, int]:
    """
    Find protein modules/communities using the Louvain algorithm.
    
    Args:
        network: NetworkX Graph representing the interaction network
        resolution: Resolution parameter (higher values give smaller communities)
        
    Returns:
        Dictionary mapping protein IDs to community assignments
    """
    try:
        import community as community_louvain
        
        # Apply Louvain community detection
        partition = community_louvain.best_partition(network, resolution=resolution)
        
        # Count communities
        num_communities = len(set(partition.values()))
        logger.info(f"Found {num_communities} communities")
        
        return partition
    except ImportError:
        logger.warning("python-louvain package not found, falling back to spectral clustering")
        return perform_spectral_clustering(network)