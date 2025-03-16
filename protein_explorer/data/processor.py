"""
Functions for processing protein data.
"""

import numpy as np
import logging
from typing import Dict, List, Tuple, Optional
from Bio.PDB import PDBParser, Select
import io

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def parse_pdb_structure(pdb_data: str) -> Dict:
    """
    Parse PDB structure data.
    
    Args:
        pdb_data: PDB format data as string
        
    Returns:
        Dictionary with parsed structure information
    """
    # Create a PDBParser object
    parser = PDBParser(QUIET=True)
    
    # Parse the PDB data from a string
    structure = parser.get_structure("protein", io.StringIO(pdb_data))
    
    # Extract basic information
    chains = list(structure.get_chains())
    residues = list(structure.get_residues())
    atoms = list(structure.get_atoms())
    
    # Create result dictionary
    result = {
        "num_chains": len(chains),
        "num_residues": len(residues),
        "num_atoms": len(atoms),
        "chains": {},
    }
    
    # Extract chain information
    for chain in chains:
        chain_id = chain.get_id()
        chain_residues = list(chain.get_residues())
        
        result["chains"][chain_id] = {
            "num_residues": len(chain_residues),
            "sequence": "".join(get_residue_sequence(chain_residues))
        }
    
    return result

def get_residue_sequence(residues: List) -> List[str]:
    """
    Convert residue objects to one-letter amino acid codes.
    
    Args:
        residues: List of Biopython residue objects
        
    Returns:
        List of one-letter amino acid codes
    """
    # Dictionary mapping three-letter codes to one-letter codes
    three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    
    sequence = []
    for residue in residues:
        res_name = residue.get_resname()
        if res_name in three_to_one:
            sequence.append(three_to_one[res_name])
        else:
            sequence.append('X')  # Unknown amino acid
            
    return sequence

def extract_coordinates(pdb_data: str, 
                      chain_id: Optional[str] = None, 
                      atom_type: str = "CA") -> np.ndarray:
    """
    Extract atom coordinates from PDB structure.
    
    Args:
        pdb_data: PDB format data as string
        chain_id: Chain ID to extract (None for all chains)
        atom_type: Atom type to extract (default: "CA" for alpha carbons)
        
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
                # Check if the requested atom exists in this residue
                if atom_type in residue:
                    atom = residue[atom_type]
                    coordinates.append(atom.get_coord())
    
    # Convert to numpy array
    return np.array(coordinates)

def parse_interaction_data(interaction_data: Dict) -> Dict:
    """
    Process protein interaction data.
    
    Args:
        interaction_data: Dictionary of protein interactions
        
    Returns:
        Processed interaction data
    """
    # Sort interactions by score (highest first)
    sorted_interactions = sorted(
        interaction_data.items(), 
        key=lambda x: x[1], 
        reverse=True
    )
    
    # Create result dictionary
    result = {
        "total_interactions": len(sorted_interactions),
        "high_confidence": len([i for i in sorted_interactions if i[1] >= 0.9]),
        "medium_confidence": len([i for i in sorted_interactions if 0.7 <= i[1] < 0.9]),
        "low_confidence": len([i for i in sorted_interactions if i[1] < 0.7]),
        "interactions": dict(sorted_interactions)
    }
    
    return result