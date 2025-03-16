"""
Functions for building and traversing protein-protein interaction networks.
"""

import networkx as nx
from collections import deque
import logging
from typing import Dict, List, Set, Tuple, Optional

from protein_explorer.data.scaffold import get_protein_interactions

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def build_interaction_network(seed_proteins: List[str], 
                             max_depth: int = 2,
                             confidence_threshold: float = 0.7) -> nx.Graph:
    """
    Build a protein-protein interaction network starting from seed proteins.
    
    Args:
        seed_proteins: List of UniProt IDs to use as starting points
        max_depth: Maximum traversal depth from seed proteins
        confidence_threshold: Minimum confidence score for interactions
        
    Returns:
        NetworkX Graph object representing the interaction network
    """
    # Create empty graph
    network = nx.Graph()
    
    # Add seed proteins as nodes
    for protein in seed_proteins:
        network.add_node(protein, depth=0, seed=True)
    
    # Track visited proteins to avoid cycles
    visited = set(seed_proteins)
    
    # Initialize queue with seed proteins and their depth
    queue = deque([(protein, 0) for protein in seed_proteins])
    
    while queue:
        current_protein, current_depth = queue.popleft()
        
        # Stop if reached maximum depth
        if current_depth >= max_depth:
            continue
        
        # Get interactions for current protein
        try:
            interactions = get_protein_interactions(
                current_protein, 
                confidence_score=confidence_threshold
            )
        except Exception as e:
            logger.error(f"Error getting interactions for {current_protein}: {e}")
            continue
        
        # Add interactions to network
        for target, score in interactions.items():
            # Add edge with confidence score
            network.add_edge(current_protein, target, confidence=score)
            
            # Add target node
            if target not in visited:
                network.add_node(target, depth=current_depth + 1, seed=False)
                visited.add(target)
                queue.append((target, current_depth + 1))
    
    logger.info(f"Built network with {network.number_of_nodes()} nodes and {network.number_of_edges()} edges")
    return network

def bfs_traverse(network: nx.Graph, 
                start_node: str, 
                max_depth: int = None) -> Dict[str, int]:
    """
    Perform breadth-first search on a protein interaction network.
    
    Args:
        network: NetworkX Graph representing the interaction network
        start_node: UniProt ID of the starting protein
        max_depth: Maximum traversal depth (None for no limit)
        
    Returns:
        Dictionary mapping node IDs to their depth from start_node
    """
    if start_node not in network:
        raise ValueError(f"Start node {start_node} not in network")
    
    # Initialize BFS
    queue = deque([(start_node, 0)])  # (node, depth)
    visited = {start_node: 0}  # node: depth
    
    while queue:
        current_node, current_depth = queue.popleft()
        
        # Stop if reached maximum depth
        if max_depth is not None and current_depth >= max_depth:
            continue
        
        # Process neighbors
        for neighbor in network.neighbors(current_node):
            if neighbor not in visited:
                new_depth = current_depth + 1
                visited[neighbor] = new_depth
                queue.append((neighbor, new_depth))
    
    return visited

def get_subnetwork(network: nx.Graph, 
                  proteins: List[str], 
                  include_neighbors: bool = False) -> nx.Graph:
    """
    Extract a subnetwork containing specified proteins.
    
    Args:
        network: NetworkX Graph representing the interaction network
        proteins: List of UniProt IDs to include
        include_neighbors: Whether to include direct neighbors of specified proteins
        
    Returns:
        NetworkX Graph representing the subnetwork
    """
    # Find nodes to include
    nodes_to_include = set(proteins)
    
    # Add neighbors if requested
    if include_neighbors:
        for protein in proteins:
            if protein in network:
                nodes_to_include.update(network.neighbors(protein))
    
    # Extract subgraph
    subnetwork = network.subgraph(nodes_to_include).copy()
    
    logger.info(f"Created subnetwork with {subnetwork.number_of_nodes()} nodes and {subnetwork.number_of_edges()} edges")
    return subnetwork