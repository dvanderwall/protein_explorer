"""
Functions for analyzing protein structures using linear algebra techniques.
"""

import numpy as np
from scipy import linalg
import logging
from typing import Dict, List, Tuple, Optional
from Bio.PDB import PDBParser, Superimposer
import io

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def extract_coordinates(pdb_data: str, 
                      atom_type: str = "CA", 
                      chain_id: Optional[str] = None) -> np.ndarray:
    """
    Extract atom coordinates from PDB structure.
    
    Args:
        pdb_data: PDB format data as string
        atom_type: Atom type to extract (default: "CA" for alpha carbons)
        chain_id: Chain ID to extract (None for all chains)
        
    Returns:
        Numpy array of 3D coordinates
    """
    # Create a PDBParser object
    parser = PDBParser(QUIET=True)
    
    # Parse the PDB data
    structure = parser.get_structure("protein", io.StringIO(pdb_data))
    
    # Extract coordinates
    coordinates = []
    
    for model in structure:
        for chain in model:
            # Skip if not the requested chain
            if chain_id and chain.get_id() != chain_id:
                continue
                
            for residue in chain:
                # Skip non-standard residues and hetero-residues
                if residue.get_id()[0] != " ":
                    continue
                    
                # Check if the requested atom exists in this residue
                if atom_type in residue:
                    atom = residue[atom_type]
                    coordinates.append(atom.get_coord())
    
    # Convert to numpy array
    return np.array(coordinates)

def calculate_pca(coordinates: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Perform Principal Component Analysis on protein structure coordinates.
    
    Args:
        coordinates: Nx3 numpy array of atom coordinates
        
    Returns:
        Tuple of (eigenvalues, eigenvectors, projected_coordinates)
    """
    if len(coordinates) < 3:
        raise ValueError("Need at least 3 coordinates for PCA")
    
    # Center the coordinates
    mean_coord = np.mean(coordinates, axis=0)
    centered_coords = coordinates - mean_coord
    
    # Calculate covariance matrix
    cov_matrix = np.cov(centered_coords, rowvar=False)
    
    # Calculate eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
    
    # Sort by eigenvalues in descending order
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    # Project data onto principal components
    projected_coords = np.dot(centered_coords, eigenvectors)
    
    # Calculate explained variance
    total_var = np.sum(eigenvalues)
    explained_var = eigenvalues / total_var
    
    logger.info(f"Explained variance ratio: {explained_var}")
    
    return eigenvalues, eigenvectors, projected_coords

def calculate_distance_matrix(coordinates: np.ndarray) -> np.ndarray:
    """
    Calculate pairwise distance matrix between atoms.
    
    Args:
        coordinates: Nx3 numpy array of atom coordinates
        
    Returns:
        NxN numpy array of pairwise distances
    """
    # Number of atoms
    n_atoms = coordinates.shape[0]
    
    # Initialize distance matrix
    distance_matrix = np.zeros((n_atoms, n_atoms))
    
    # Calculate pairwise distances
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            # Euclidean distance
            dist = np.linalg.norm(coordinates[i] - coordinates[j])
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist
    
    return distance_matrix

def compute_structural_similarity(coords_a: np.ndarray, 
                                coords_b: np.ndarray) -> Dict:
    """
    Calculate structural similarity between two protein structures.
    
    Args:
        coords_a: Nx3 numpy array of coordinates for structure A
        coords_b: Mx3 numpy array of coordinates for structure B
        
    Returns:
        Dictionary with similarity metrics
    """
    # Check if structures have the same number of atoms
    if coords_a.shape[0] != coords_b.shape[0]:
        logger.warning("Structures have different number of atoms, using minimum length")
        min_length = min(coords_a.shape[0], coords_b.shape[0])
        coords_a = coords_a[:min_length]
        coords_b = coords_b[:min_length]
    
    # Create Superimposer object
    sup = Superimposer()
    
    # Set the two coordinate sets
    sup.set_atoms(coords_a, coords_b)
    
    # Calculate RMSD
    rmsd = sup.rms
    
    # Apply rotation and translation to coords_b
    rotated_coords = coords_b.copy()
    sup.apply(rotated_coords)
    
    # Calculate per-residue RMSD
    per_residue_rmsd = np.sqrt(np.sum((coords_a - rotated_coords)**2, axis=1))
    
    # Calculate TM-score (simplified)
    # TM-score ranges from 0 to 1, with 1 indicating perfect match
    d0 = 1.24 * (coords_a.shape[0] - 15)**(1/3) - 1.8
    tm_score = np.mean(1 / (1 + (per_residue_rmsd / d0)**2))
    
    return {
        "rmsd": rmsd,
        "tm_score": tm_score,
        "per_residue_rmsd": per_residue_rmsd.tolist()
    }

def identify_domains(coordinates: np.ndarray, 
                   distance_threshold: float = 8.0) -> List[List[int]]:
    """
    Identify potential protein domains using hierarchical clustering.
    
    Args:
        coordinates: Nx3 numpy array of atom coordinates
        distance_threshold: Maximum distance for atoms to be in the same domain
        
    Returns:
        List of lists, where each sublist contains residue indices in a domain
    """
    from scipy.cluster.hierarchy import linkage, fcluster
    
    # Calculate distance matrix
    distance_matrix = calculate_distance_matrix(coordinates)
    
    # Perform hierarchical clustering
    Z = linkage(distance_matrix, method='ward')
    
    # Cluster based on distance threshold
    clusters = fcluster(Z, distance_threshold, criterion='distance')
    
    # Group residues by cluster
    domains = {}
    for i, cluster_id in enumerate(clusters):
        if cluster_id not in domains:
            domains[cluster_id] = []
        domains[cluster_id].append(i)
    
    # Convert to list of lists
    domain_list = list(domains.values())
    
    # Sort by domain size (largest first)
    domain_list.sort(key=len, reverse=True)
    
    return domain_list

def calculate_contact_map(coordinates: np.ndarray, 
                        distance_cutoff: float = 8.0) -> np.ndarray:
    """
    Calculate a contact map for the protein structure.
    
    Args:
        coordinates: Nx3 numpy array of atom coordinates
        distance_cutoff: Maximum distance for residues to be in contact
        
    Returns:
        NxN numpy array where 1 indicates residues in contact, 0 otherwise
    """
    # Calculate distance matrix
    distance_matrix = calculate_distance_matrix(coordinates)
    
    # Create contact map
    contact_map = (distance_matrix <= distance_cutoff).astype(int)
    
    # Remove self-contacts
    np.fill_diagonal(contact_map, 0)
    
    return contact_map