"""
Functions for finding paths and relationships in protein interaction networks.
"""

import networkx as nx
from collections import deque
import logging
from typing import Dict, List, Set, Tuple, Optional

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def find_path(network: nx.Graph, 
             protein_a: str, 
             protein_b: str,
             weight: Optional[str] = None) -> List[str]:
    """
    Find the shortest path between two proteins in the interaction network.
    
    Args:
        network: NetworkX Graph representing the interaction network
        protein_a: UniProt ID of the first protein
        protein_b: UniProt ID of the second protein
        weight: Edge attribute to use as weight (None for unweighted)
        
    Returns:
        List of proteins in the path, or empty list if no path exists
    """
    if protein_a not in network:
        logger.error(f"Protein {protein_a} not found in network")
        return []
    
    if protein_b not in network:
        logger.error(f"Protein {protein_b} not found in network")
        return []
    
    try:
        if weight:
            # Use weighted shortest path
            path = nx.shortest_path(network, protein_a, protein_b, weight=weight)
        else:
            # Use unweighted shortest path (BFS)
            path = nx.shortest_path(network, protein_a, protein_b)
            
        return path
    except nx.NetworkXNoPath:
        logger.info(f"No path found between {protein_a} and {protein_b}")
        return []

def find_common_interactors(network: nx.Graph, protein_list: List[str]) -> Dict[str, Set[str]]:
    """
    Identify proteins that interact with multiple proteins from the input list.
    
    Args:
        network: NetworkX Graph representing the interaction network
        protein_list: List of UniProt IDs to analyze
        
    Returns:
        Dictionary mapping interactor IDs to the set of proteins they interact with
    """
    # Filter protein list to those present in the network
    valid_proteins = [p for p in protein_list if p in network]
    
    if not valid_proteins:
        logger.warning("None of the provided proteins found in network")
        return {}
    
    # Find all interactors for each protein
    interactors = {}
    
    for protein in valid_proteins:
        neighbors = set(network.neighbors(protein))
        
        for neighbor in neighbors:
            if neighbor in interactors:
                interactors[neighbor].add(protein)
            else:
                interactors[neighbor] = {protein}
    
    # Filter for interactors that connect to multiple input proteins
    common_interactors = {i: proteins for i, proteins in interactors.items() 
                         if len(proteins) > 1 and i not in protein_list}
    
    return common_interactors

def find_bridges(network: nx.Graph, protein_list: List[str]) -> List[Tuple[str, str]]:
    """
    Find bridges (proteins that mediate connections) between input proteins.
    
    Args:
        network: NetworkX Graph representing the interaction network
        protein_list: List of UniProt IDs to analyze
        
    Returns:
        List of bridges (pairs of proteins connected through a third protein)
    """
    # Check if proteins exist in network
    valid_proteins = [p for p in protein_list if p in network]
    
    if len(valid_proteins) < 2:
        logger.warning("Need at least two valid proteins to find bridges")
        return []
    
    bridges = []
    
    # Check each pair of proteins
    for i, protein_a in enumerate(valid_proteins[:-1]):
        for protein_b in valid_proteins[i+1:]:
            # Skip if directly connected
            if network.has_edge(protein_a, protein_b):
                continue
                
            # Find paths of length 2 (with exactly one intermediate node)
            paths = []
            
            for neighbor in network.neighbors(protein_a):
                if network.has_edge(neighbor, protein_b) and neighbor not in valid_proteins:
                    paths.append((protein_a, neighbor, protein_b))
            
            bridges.extend(paths)
    
    return bridges

def find_cliques(network: nx.Graph, min_size: int = 3) -> List[List[str]]:
    """
    Find protein cliques (fully connected subgraphs) in the network.
    
    Args:
        network: NetworkX Graph representing the interaction network
        min_size: Minimum size of cliques to find
        
    Returns:
        List of cliques (each a list of protein IDs)
    """
    # Find all maximal cliques
    cliques = list(nx.find_cliques(network))
    
    # Filter by minimum size
    large_cliques = [c for c in cliques if len(c) >= min_size]
    
    # Sort by size (largest first)
    large_cliques.sort(key=len, reverse=True)
    
    return large_cliques