"""
Phosphosite Analyzer module - Handles phosphosite structural analysis for protein structures.

This module provides functions to analyze phosphorylation sites in proteins
and find structural similarities with other known sites.
"""

import os
import pandas as pd
import requests
from typing import Dict, List, Optional, Union, Tuple
import logging
from protein_explorer.analysis.phospho import analyze_phosphosites
from protein_explorer.data.scaffold import get_protein_by_id, get_alphafold_structure

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Global variables to store loaded data
STRUCTURAL_SIMILARITY_DF = None
PHOSPHOSITE_SUPP_DATA = None
PHOSPHOSITE_SUPP_DICT = None

def preload_structural_data(file_path: str = None) -> None:
    """
    Preload structural similarity data at application startup.
    
    Args:
        file_path: Path to the data file (feather preferred, parquet as fallback)
    """
    global STRUCTURAL_SIMILARITY_DF
    
    if STRUCTURAL_SIMILARITY_DF is not None:
        logger.info("Structural data already loaded")
        return
    
    # Find the data file
    if file_path is None:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(current_dir)
        # Check for feather file first
        feather_file = os.path.join(parent_dir, 'Combined_Kinome_10A_Master_Filtered_2.feather')
        if os.path.exists(feather_file):
            file_path = feather_file
            is_feather = True
        else:
            # Fall back to parquet
            parquet_file = os.path.join(parent_dir, 'Combined_Kinome_10A_Master_Filtered_2.parquet')
            if os.path.exists(parquet_file):
                file_path = parquet_file
                is_feather = False
            else:
                logger.warning("No structural data file found, will load on first request")
                return
    else:
        # Determine file type from extension
        is_feather = file_path.endswith('.feather')
    
    try:
        # Load the data
        logger.info(f"Preloading structural data from: {file_path}")
        import pandas as pd
        if is_feather:
            STRUCTURAL_SIMILARITY_DF = pd.read_feather(file_path)
        else:
            STRUCTURAL_SIMILARITY_DF = pd.read_parquet(file_path)
        
        # Create index for faster querying
        logger.info("Creating query index")
        STRUCTURAL_SIMILARITY_DF.set_index('Query', drop=False, inplace=True)
        
        logger.info(f"Successfully preloaded {len(STRUCTURAL_SIMILARITY_DF)} structural similarity records")
    except Exception as e:
        logger.error(f"Error preloading structural data: {e}")
        logger.warning("Will attempt to load data on first request")


def load_phosphosite_supp_data(file_path: str = None) -> pd.DataFrame:
    """
    Load the phosphosite supplementary data.
    
    Args:
        file_path: Path to the supplementary data file
        
    Returns:
        Pandas DataFrame with supplementary data
    """
    global PHOSPHOSITE_SUPP_DATA, PHOSPHOSITE_SUPP_DICT
    
    if PHOSPHOSITE_SUPP_DATA is not None:
        logger.info("Using cached phosphosite supplementary data")
        return PHOSPHOSITE_SUPP_DATA
    
    # Find the data file if not provided
    if file_path is None:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(current_dir)
        file_path = os.path.join(parent_dir, 'PhosphositeSuppData.feather')
    
    if not os.path.exists(file_path):
        logger.warning(f"Supplementary data file not found: {file_path}")
        return None
    
    try:
        # Load the data
        logger.info(f"Loading phosphosite supplementary data from: {file_path}")
        PHOSPHOSITE_SUPP_DATA = pd.read_feather(file_path)
        
        # Create dictionary for faster lookups
        logger.info("Creating lookup dictionary")
        PHOSPHOSITE_SUPP_DICT = {row['site_id']: row.to_dict() for _, row in PHOSPHOSITE_SUPP_DATA.iterrows()}
        
        logger.info(f"Successfully loaded {len(PHOSPHOSITE_SUPP_DATA)} phosphosite supplementary records")
        return PHOSPHOSITE_SUPP_DATA
    except Exception as e:
        logger.error(f"Error loading supplementary data: {e}")
        return None

        


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_structural_similarity_data(parquet_file: str = None) -> pd.DataFrame:
    """
    Load structural similarity data from feather file.
    
    Args:
        parquet_file: Path to the feather file with structural similarity data
                    (kept the parameter name for backward compatibility)
        
    Returns:
        Pandas DataFrame with structural similarity data
    """
    global STRUCTURAL_SIMILARITY_DF
    
    if STRUCTURAL_SIMILARITY_DF is not None:
        logger.info("Using cached structural similarity data")
        return STRUCTURAL_SIMILARITY_DF
    
    # Check if feather file exists, use default if not provided
    if parquet_file is None:
        # Try to find the feather file in the parent directory
        current_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(current_dir)
        feather_file = os.path.join(parent_dir, 'Combined_Kinome_10A_Master_Filtered_2.feather')
    else:
        # Replace .parquet with .feather
        feather_file = parquet_file.replace('.parquet', '.feather')
    
    if not os.path.exists(feather_file):
        logger.error(f"Feather file not found: {feather_file}")
        raise FileNotFoundError(f"Structural similarity data file not found: {feather_file}")
    
    try:
        # Read feather file
        logger.info(f"Reading feather file: {feather_file}")
        import pandas as pd
        STRUCTURAL_SIMILARITY_DF = pd.read_feather(feather_file)
        
        # Create indexes for faster querying
        logger.info("Creating query index")
        STRUCTURAL_SIMILARITY_DF.set_index('Query', drop=False, inplace=True)
        
        logger.info(f"Loaded {len(STRUCTURAL_SIMILARITY_DF)} structural similarity records")
        return STRUCTURAL_SIMILARITY_DF
    except Exception as e:
        logger.error(f"Error reading feather file: {e}")
        raise ValueError(f"Error reading structural similarity data: {e}")


def get_protein_data(identifier: str, id_type: str = 'uniprot') -> Dict:
    """
    Retrieve protein data by UniProt ID or gene symbol.
    
    Args:
        identifier: UniProt ID or gene symbol
        id_type: 'uniprot' or 'gene'
        
    Returns:
        Dictionary with protein data
    """
    logger.info(f"Getting protein data for {identifier} (type: {id_type})")
    
    try:
        if id_type.lower() == 'uniprot':
            protein_data = get_protein_by_id(uniprot_id=identifier)
        else:
            protein_data = get_protein_by_id(gene_symbol=identifier)
            
        # Extract protein info
        uniprot_id = protein_data.get('uniprot_id')
        gene_symbol = protein_data.get('gene_symbol', 'Unknown')
        name = protein_data.get('metadata', {}).get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown Protein')
        
        protein_info = {
            'uniprot_id': uniprot_id,
            'gene_symbol': gene_symbol,
            'name': name,
            'full_data': protein_data
        }
        
        return protein_info
    except Exception as e:
        logger.error(f"Error retrieving protein data: {e}")
        raise ValueError(f"Error retrieving protein data: {e}")

    
def get_phosphosite_data(site_id: str) -> Optional[Dict]:
    """
    Get supplementary data for a specific phosphosite.
    
    Args:
        site_id: The site ID in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary with supplementary data or None if not found
    """
    global PHOSPHOSITE_SUPP_DATA, PHOSPHOSITE_SUPP_DICT
    
    # Load data if not already loaded
    if PHOSPHOSITE_SUPP_DATA is None:
        load_phosphosite_supp_data()
    
    # If still None, return None
    if PHOSPHOSITE_SUPP_DATA is None:
        return None
    
    # Use dictionary for faster lookup
    if PHOSPHOSITE_SUPP_DICT is not None:
        return PHOSPHOSITE_SUPP_DICT.get(site_id)
    
    # Fall back to DataFrame lookup
    try:
        site_data = PHOSPHOSITE_SUPP_DATA[PHOSPHOSITE_SUPP_DATA['site_id'] == site_id]
        if not site_data.empty:
            return site_data.iloc[0].to_dict()
    except Exception as e:
        logger.error(f"Error getting phosphosite data: {e}")
    
    return None

def enhance_phosphosite(phosphosite: Dict, uniprot_id: str) -> Dict:
    """
    Enhance a phosphosite dictionary with supplementary data.
    
    Args:
        phosphosite: Dictionary with phosphosite info
        uniprot_id: UniProt ID
        
    Returns:
        Enhanced phosphosite dictionary
    """
    if 'resno' not in phosphosite:
        return phosphosite
    
    # Get site ID
    site_id = f"{uniprot_id}_{phosphosite['resno']}"
    
    # Get supplementary data
    supp_data = get_phosphosite_data(site_id)
    if not supp_data:
        return phosphosite
    
    # Create a new dictionary to avoid modifying the original
    enhanced_site = phosphosite.copy()
    
    # Enhance with supplementary data
    if 'motif_plddt' in supp_data and supp_data['motif_plddt'] is not None:
        enhanced_site['mean_plddt'] = f"{supp_data['motif_plddt']:.1f}"
    
    if 'nearby_count' in supp_data and supp_data['nearby_count'] is not None:
        enhanced_site['nearby_count'] = supp_data['nearby_count']
    
    if 'motif' in supp_data and supp_data['motif'] is not None:
        enhanced_site['motif'] = supp_data['motif']
    
    # Add any other supplementary fields
    for key in ['site_plddt', 'surface_accessibility', 'secondary_structure']:
        if key in supp_data and supp_data[key] is not None:
            enhanced_site[key] = supp_data[key]
    
    return enhanced_site

def enhance_structural_matches(matches: List[Dict], site: str) -> List[Dict]:
    """
    Enhance structural matches with supplementary data.
    
    Args:
        matches: List of structural match dictionaries
        site: Query site string for logging
        
    Returns:
        Enhanced list of matches
    """
    if not matches:
        return matches
    
    logger.info(f"Enhancing {len(matches)} structural matches for {site}")
    
    # Ensure supplementary data is loaded
    if PHOSPHOSITE_SUPP_DATA is None:
        load_phosphosite_supp_data()
    
    enhanced_matches = []
    for match in matches:
        # Skip self-matches (RMSD ≈ 0)
        if match.get('rmsd', 0) < 0.01:
            continue
        
        # Get target info
        target_uniprot = match.get('target_uniprot')
        target_site = match.get('target_site')
        
        # Parse site number from target_site
        import re
        site_match = re.match(r'(\d+)', target_site)
        if site_match:
            resno = int(site_match.group(1))
            target_id = f"{target_uniprot}_{resno}"
            
            # Get supplementary data
            target_supp = get_phosphosite_data(target_id)
            if target_supp:
                # Create enhanced match
                enhanced_match = match.copy()
                
                # Add supplementary data
                if 'motif_plddt' in target_supp and target_supp['motif_plddt'] is not None:
                    enhanced_match['plddt'] = f"{target_supp['motif_plddt']:.1f}"
                
                if 'nearby_count' in target_supp and target_supp['nearby_count'] is not None:
                    enhanced_match['nearby_count'] = target_supp['nearby_count']
                
                if 'motif' in target_supp and target_supp['motif'] is not None:
                    enhanced_match['motif'] = target_supp['motif']
                
                # Add additional fields
                for key in ['site_plddt', 'surface_accessibility', 'secondary_structure']:
                    if key in target_supp and target_supp[key] is not None:
                        enhanced_match[key] = target_supp[key]
                
                enhanced_matches.append(enhanced_match)
                continue
        
        # If no supplementary data, just add the original match
        enhanced_matches.append(match)
    
    return enhanced_matches

def get_phosphosites(uniprot_id: str) -> List[Dict]:
    """
    Analyze potential phosphorylation sites for a protein, with supplementary data.
    
    Args:
        uniprot_id: UniProt ID of the protein
        
    Returns:
        List of dictionaries with phosphosite information
    """
    logger.info(f"Analyzing phosphosites for {uniprot_id}")
    
    try:
        # Get protein data
        protein_data = get_protein_by_id(uniprot_id=uniprot_id)
        
        # Get sequence
        sequence = protein_data.get('metadata', {}).get('sequence', {}).get('value')
        if not sequence:
            logger.warning(f"Protein sequence not found for {uniprot_id}")
            raise ValueError(f"Protein sequence not found for {uniprot_id}")
            
        # Get structure
        structure = get_alphafold_structure(uniprot_id)
        if not structure:
            logger.warning(f"Protein structure not found for {uniprot_id}. Checking alternative sources...")
            
            # Try a mock structure for testing purposes when no real structure is available
            # This could be replaced with other structure sources (PDB, etc.) in a production environment
            mock_structure = generate_mock_structure(sequence)
            if mock_structure:
                logger.info(f"Using mock structure for {uniprot_id}")
                structure = mock_structure
            else:
                raise ValueError(f"Protein structure not found for {uniprot_id}")
        
        # Analyze phosphosites
        phosphosites = analyze_phosphosites(sequence, structure, uniprot_id)
        
        # Enhance with supplementary data
        enhanced_sites = []
        for site in phosphosites:
            enhanced_site = enhance_phosphosite(site, uniprot_id)
            enhanced_sites.append(enhanced_site)
        
        return enhanced_sites
    except Exception as e:
        logger.error(f"Error analyzing phosphosites: {e}")
        raise ValueError(f"Error analyzing phosphosites: {e}")

    

def generate_mock_structure(sequence: str) -> Optional[str]:
    """
    Generate a mock PDB structure for cases where the AlphaFold structure is not available.
    This is for demonstration purposes only and should be replaced with real structures in production.
    
    Args:
        sequence: Protein sequence
        
    Returns:
        PDB format structure as string, or None if generation fails
    """
    try:
        # Create a very basic linear structure
        # This is extremely simplified and not biologically accurate
        pdb_lines = []
        
        # PDB header
        pdb_lines.append("HEADER    MOCK STRUCTURE")
        pdb_lines.append("TITLE     MOCK STRUCTURE FOR SEQUENCE")
        
        # Add atoms - just alpha carbons in a straight line
        atom_num = 1
        for i, aa in enumerate(sequence):
            x = i * 3.8  # ~3.8Å is typical CA-CA distance
            y = 0
            z = 0
            
            # B-factor (PLDDT) set to 70 (medium confidence) for all residues
            b_factor = 70.0
            
            line = f"ATOM  {atom_num:5d}  CA  {aa}   A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{b_factor:6.2f}           C  "
            pdb_lines.append(line)
            atom_num += 1
        
        # End of file
        pdb_lines.append("END")
        
        return "\n".join(pdb_lines)
    except Exception as e:
        logger.error(f"Error generating mock structure: {e}")
        return None


def find_structural_matches(uniprot_id: str, phosphosites: List[Dict], 
                           parquet_file: str = None, top_n: int = None) -> Dict[str, List[Dict]]:
    """
    Find structural matches for phosphosites in the kinome dataset.
    
    Args:
        uniprot_id: UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries from analyze_phosphosites
        parquet_file: Path to the data file (parameter kept for backward compatibility)
        top_n: Number of top matches to return per site
        
    Returns:
        Dictionary mapping site IDs to lists of match dictionaries
    """
    global STRUCTURAL_SIMILARITY_DF
    logger.info(f"Finding structural matches for {uniprot_id}")
    
    try:
        # Use preloaded data if available, otherwise load it now
        if STRUCTURAL_SIMILARITY_DF is None:
            preload_structural_data(parquet_file)
        
        # If still None after trying to load, raise error
        if STRUCTURAL_SIMILARITY_DF is None:
            logger.error("Structural similarity data not available")
            raise ValueError("Structural similarity data not available")
        
        df = STRUCTURAL_SIMILARITY_DF
        
        # Create site IDs in the format UniprotID_ResNo
        site_ids = [f"{uniprot_id}_{site['resno']}" for site in phosphosites]
        
        # Find matches and sort by RMSD
        matches = []
        for site_id in site_ids:
            # Use efficient lookup with index
            if site_id in df.index:
                if top_n is not None:
                    # Only take top N matches
                    site_matches = df.loc[[site_id]].sort_values('RMSD').head(top_n)
                else:
                    # Take all matches, still sorted by RMSD
                    site_matches = df.loc[[site_id]].sort_values('RMSD')
                
                if not site_matches.empty:
                    for _, row in site_matches.iterrows():
                        query_parts = row['Query'].split('_')
                        target_parts = row['Target'].split('_')
                        
                        # Only add if we can parse the site numbers
                        if len(query_parts) > 1 and len(target_parts) > 1:
                            try:
                                query_site = int(query_parts[-1])
                                target_uniprot = target_parts[0]
                                target_site = int(target_parts[-1])
                                
                                # Find the corresponding site data
                                site_data = next((s for s in phosphosites if s['resno'] == query_site), None)
                                site_type = site_data['site'][0] if site_data else '?'
                                
                                matches.append({
                                    'query_uniprot': uniprot_id,
                                    'query_site': f"{site_type}{query_site}",
                                    'target_uniprot': target_uniprot,
                                    'target_site': target_parts[-1],
                                    'rmsd': row['RMSD']
                                })
                            except (ValueError, IndexError) as e:
                                logger.warning(f"Error parsing site ID: {e}")
        
        # Group by query site
        structural_matches = {}
        for match in matches:
            query_site = match['query_site']
            if query_site not in structural_matches:
                structural_matches[query_site] = []
            structural_matches[query_site].append(match)
        
        return structural_matches
    except Exception as e:
        logger.error(f"Error finding structural matches: {e}")
        raise ValueError(f"Error finding structural matches: {e}")

def analyze_protein(identifier: str, id_type: str = 'uniprot', 
                   parquet_file: str = None) -> Dict:
    """
    Complete analysis of phosphosites and structural matches for a protein.
    
    Args:
        identifier: UniProt ID or gene symbol
        id_type: 'uniprot' or 'gene'
        parquet_file: Path to the parquet file with structural similarity data
        
    Returns:
        Dictionary with protein info, phosphosites, and structural matches
    """
    try:
        # Get protein data
        protein_info = get_protein_data(identifier, id_type)
        uniprot_id = protein_info['uniprot_id']
        
        # Try to get phosphosites
        phosphosites = []
        structural_matches = None
        error_message = None
        
        try:
            # Get phosphosites with supplementary data
            phosphosites = get_phosphosites(uniprot_id)
            
            # Find structural matches
            try:
                structural_matches = find_structural_matches(uniprot_id, phosphosites, parquet_file)
                
                # Enhance with supplementary data
                for site, matches in structural_matches.items():
                    structural_matches[site] = enhance_structural_matches(matches, site)
                    
            except FileNotFoundError:
                error_message = "Structural similarity data file not found"
                logger.warning(error_message)
            except Exception as e:
                error_message = f"Error analyzing structural matches: {str(e)}"
                logger.error(error_message)
        except Exception as e:
            error_message = f"Error analyzing phosphosites: {str(e)}"
            logger.error(f"Error analyzing phosphosites: {e}")
        
        # Compile results
        results = {
            'protein_info': protein_info,
            'phosphosites': phosphosites,
            'structural_matches': structural_matches,
            'error': error_message
        }
        
        return results
    except Exception as e:
        logger.error(f"Error in complete analysis: {e}")
        raise ValueError(f"Error in complete analysis: {e}")
    

def analyze_phosphosite_context(structure_data, site_number, site_type):
    """
    Analyze structural context around a phosphorylation site.
    
    Args:
        structure_data: PDB format data as string
        site_number: The residue number
        site_type: The residue type (S, T, or Y)
        
    Returns:
        Dictionary with structural context information
    """
    from Bio.PDB import PDBParser, NeighborSearch, Selection, Vector
    import io
    import numpy as np
    
    # Parse the PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", io.StringIO(structure_data))
    
    # Get all atoms
    all_atoms = list(structure.get_atoms())
    ns = NeighborSearch(all_atoms)
    
    # Find the target residue
    target_residue = None
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[1] == site_number:
                    target_residue = residue
                    break
    
    if not target_residue:
        return {"error": f"Site {site_type}{site_number} not found in structure"}
    
    # Get center atom for the residue (CA or first atom)
    if 'CA' in target_residue:
        center_atom = target_residue['CA']
    else:
        # Use the first atom if CA not available
        center_atom = next(target_residue.get_atoms())
    
    # Get the residue coordinates
    center_coords = center_atom.get_coord()
    
    # Find nearby residues (8Å radius)
    nearby_atoms = ns.search(center_coords, 8.0)
    
    # Group nearby atoms by residue
    nearby_residues = {}
    for atom in nearby_atoms:
        residue = atom.get_parent()
        resno = residue.get_id()[1]
        
        # Skip the target residue itself
        if resno == site_number:
            continue
            
        residue_name = residue.get_resname()
        
        # Add to nearby residues dictionary
        if resno not in nearby_residues:
            nearby_residues[resno] = {
                "resname": residue_name,
                "atoms": [],
                "min_distance": float('inf')
            }
            
        # Store the atom and its distance
        dist = np.linalg.norm(atom.get_coord() - center_coords)
        nearby_residues[resno]["atoms"].append({
            "atom_name": atom.get_name(),
            "distance": dist
        })
        
        # Update minimum distance
        if dist < nearby_residues[resno]["min_distance"]:
            nearby_residues[resno]["min_distance"] = dist
            
    # Sort nearby residues by distance
    sorted_nearby = sorted(
        nearby_residues.items(), 
        key=lambda x: x[1]["min_distance"]
    )
    
    # Prepare results
    nearby_info = [
        {
            "resno": resno,
            "resname": data["resname"],
            "min_distance": round(data["min_distance"], 2),
            "atoms": len(data["atoms"])
        }
        for resno, data in sorted_nearby
    ]
    
    # Count amino acid types within contact distance (5Å)
    amino_acid_groups = {
        "polar": ["SER", "THR", "TYR", "CYS", "ASN", "GLN"],
        "nonpolar": ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "GLY"],
        "acidic": ["ASP", "GLU"],
        "basic": ["LYS", "ARG", "HIS"]
    }
    
    contact_counts = {group: 0 for group in amino_acid_groups}
    
    for resno, data in sorted_nearby:
        if data["min_distance"] <= 5.0:  # Only count residues within 5Å
            for group, residues in amino_acid_groups.items():
                if data["resname"] in residues:
                    contact_counts[group] += 1
                    break
                    
    # Calculate secondary structure
    try:
        from Bio.PDB.DSSP import DSSP
        model = structure[0]  # Use first model
        dssp = DSSP(model, io.StringIO(structure_data))
        
        # Get DSSP data for the target residue
        chain_id = list(model.child_dict.keys())[0]  # Get first chain ID
        dssp_key = (chain_id, site_number)
        
        if dssp_key in dssp:
            # DSSP assigns: H (alpha helix), B (beta bridge), E (strand),
            # G (3-10 helix), I (pi helix), T (turn), S (bend), or - (other)
            ss_code = dssp[dssp_key][1]
            
            # Simplify to three main categories
            if ss_code in ['H', 'G', 'I']:
                secondary_structure = 'Helix'
            elif ss_code in ['E', 'B']:
                secondary_structure = 'Sheet'
            else:
                secondary_structure = 'Loop'
        else:
            secondary_structure = 'Unknown'
    except:
        # If DSSP fails, leave as unknown
        secondary_structure = 'Unknown'
    
    # Calculate solvent accessibility
    # We'll use a simple proxy based on neighbor count
    # Fewer neighbors = more exposed
    max_neighbors = 30  # Approximate maximum reasonable number of neighbors in 8Å
    nearby_count = len(nearby_residues)
    solvent_accessibility = max(0, min(100, (max_neighbors - nearby_count) / max_neighbors * 100))
    
    # Extract B-factor (pLDDT in AlphaFold) for the target residue
    b_factors = [atom.get_bfactor() for atom in target_residue]
    mean_plddt = sum(b_factors) / len(b_factors) if b_factors else None
    
    return {
        "site": f"{site_type}{site_number}",
        "nearby_residues": nearby_info[:10],  # Show top 10
        "nearby_count": len(nearby_info),
        "contact_distribution": contact_counts,
        "secondary_structure": secondary_structure,
        "solvent_accessibility": round(solvent_accessibility, 1),
        "plddt": round(mean_plddt, 1) if mean_plddt else None
    }


def enhance_site_visualization(uniprot_id, site, supplementary_data=None):
    """
    Create an enhanced visualization of a phosphorylation site.
    Integrates supplementary structural data and highlights structural features.
    
    Args:
        uniprot_id: UniProt ID of the protein
        site: Site identifier (e.g., "S15")
        supplementary_data: Supplementary data for the site
        
    Returns:
        HTML/JavaScript code for the visualization
    """
    from protein_explorer.data.scaffold import get_alphafold_structure
    import re
    import base64
    
    # Parse site to get residue number and type
    site_match = re.match(r'([A-Z])(\d+)', site)
    if not site_match:
        return f"<div class='alert alert-danger'>Invalid site format: {site}</div>"
    
    site_type = site_match.group(1)
    site_number = int(site_match.group(2))
    
    # Get structure
    structure_data = get_alphafold_structure(uniprot_id)
    if not structure_data:
        return f"<div class='alert alert-danger'>Could not retrieve structure for {uniprot_id}</div>"
    
    # Base64 encode the structure
    pdb_base64 = base64.b64encode(structure_data.encode()).decode()
    
    # Analyze structural context if not provided
    if not supplementary_data:
        from protein_explorer.analysis.phospho_analyzer import get_phosphosite_data
        site_id = f"{uniprot_id}_{site_number}"
        supplementary_data = get_phosphosite_data(site_id)
    
    # Default values if data not available
    site_plddt = supplementary_data.get('site_plddt', 'N/A') if supplementary_data else 'N/A'
    surface_access = supplementary_data.get('surface_accessibility', 'N/A') if supplementary_data else 'N/A'
    nearby_count = supplementary_data.get('nearby_count', 'N/A') if supplementary_data else 'N/A'
    secondary_structure = supplementary_data.get('secondary_structure', 'Unknown') if supplementary_data else 'Unknown'
    
    # Create visualization code
    js_code = f"""
    <div class="card mb-4">
        <div class="card-header">
            <h5 class="mb-0">3D Site Visualization</h5>
        </div>
        <div class="card-body">
            <div class="row">
                <div class="col-md-8">
                    <div id="site-viewer" style="width: 100%; height: 450px; border: 1px solid #ddd; border-radius: 5px;"></div>
                </div>
                <div class="col-md-4">
                    <div class="card">
                        <div class="card-header">
                            <h6 class="mb-0">Site Information</h6>
                        </div>
                        <div class="card-body">
                            <p><strong>Site:</strong> {site}</p>
                            <p><strong>pLDDT:</strong> {site_plddt}</p>
                            <p><strong>Surface Accessibility:</strong> {surface_access}%</p>
                            <p><strong>Nearby Residues:</strong> {nearby_count}</p>
                            <p><strong>Secondary Structure:</strong> {secondary_structure}</p>
                        </div>
                    </div>
                    <div class="mt-3">
                        <button id="reset-view" class="btn btn-sm btn-outline-primary">Reset View</button>
                        <button id="toggle-view" class="btn btn-sm btn-outline-secondary">Full Protein</button>
                        <button id="toggle-color" class="btn btn-sm btn-outline-success">Color by Type</button>
                    </div>
                    <div class="mt-3">
                        <div class="d-flex align-items-center mb-1">
                            <div style="width:12px; height:12px; background-color:#FF4500; border-radius:50%; margin-right:5px;"></div>
                            <small>Target Site</small>
                        </div>
                        <div class="d-flex align-items-center mb-1">
                            <div style="width:12px; height:12px; background-color:#87CEFA; border-radius:50%; margin-right:5px;"></div>
                            <small>Polar Residues</small>
                        </div>
                        <div class="d-flex align-items-center mb-1">
                            <div style="width:12px; height:12px; background-color:#FFD700; border-radius:50%; margin-right:5px;"></div>
                            <small>Non-polar Residues</small>
                        </div>
                        <div class="d-flex align-items-center mb-1">
                            <div style="width:12px; height:12px; background-color:#FF6347; border-radius:50%; margin-right:5px;"></div>
                            <small>Acidic Residues</small>
                        </div>
                        <div class="d-flex align-items-center mb-1">
                            <div style="width:12px; height:12px; background-color:#98FB98; border-radius:50%; margin-right:5px;"></div>
                            <small>Basic Residues</small>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
    
    <script>
    document.addEventListener('DOMContentLoaded', function() {{
        // Initialize NGL viewer
        const viewer = new NGL.Stage('site-viewer', {{backgroundColor: "white"}});
        
        // Handle window resizing
        window.addEventListener('resize', function() {{
            viewer.handleResize();
        }});
        
        // Define residue type groupings
        const aminoAcidGroups = {{
            polar: ["SER", "THR", "TYR", "CYS", "ASN", "GLN"],
            nonpolar: ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "GLY"],
            acidic: ["ASP", "GLU"],
            basic: ["LYS", "ARG", "HIS"]
        }};
        
        // Define colors for each group
        const groupColors = {{
            polar: [135/255, 206/255, 250/255],     // Light blue
            nonpolar: [255/255, 215/255, 0/255],    // Gold
            acidic: [255/255, 99/255, 71/255],      // Tomato
            basic: [152/255, 251/255, 152/255]      // Pale green
        }};
        
        // Function to determine AA group
        function getAminoAcidGroup(resname) {{
            for (const [group, residues] of Object.entries(aminoAcidGroups)) {{
                if (residues.includes(resname)) {{
                    return group;
                }}
            }}
            return "other";
        }}
        
        // Color function based on amino acid type
        function colorByType(atom) {{
            // Special color for the target site
            if (atom.resno === {site_number}) {{
                return [1.0, 0.27, 0.0];  // #FF4500 orange-red
            }}
            
            // Color by amino acid type
            const group = getAminoAcidGroup(atom.resname);
            if (group in groupColors) {{
                return groupColors[group];
            }}
            
            // Default grey for others
            return [0.5, 0.5, 0.5];
        }}
        
        // Load structure
        const pdbBlob = new Blob([atob('{pdb_base64}')], {{type: 'text/plain'}});
        
        viewer.loadFile(pdbBlob, {{ext: 'pdb'}}).then(function(component) {{
            // Get target selection
            const siteSelection = "{site_number} and .{site_type}";
            const environmentSelection = siteSelection + " or (" + siteSelection + " around 5)";
            
            // State variables
            let isFullView = false;
            let colorMode = "element";  // "element" or "type"
            
            // Button handlers
            document.getElementById('reset-view').addEventListener('click', function() {{
                updateRepresentations();
            }});
            
            document.getElementById('toggle-view').addEventListener('click', function() {{
                isFullView = !isFullView;
                this.textContent = isFullView ? 'Site Focus' : 'Full Protein';
                updateRepresentations();
            }});
            
            document.getElementById('toggle-color').addEventListener('click', function() {{
                colorMode = colorMode === "element" ? "type" : "element";
                this.textContent = colorMode === "element" ? 'Color by Type' : 'Color by Element';
                updateRepresentations();
            }});
            
            // Update all representations based on current state
            function updateRepresentations() {{
                // Remove all existing representations
                component.removeAllRepresentations();
                
                // Add cartoon representation for entire protein
                component.addRepresentation("cartoon", {{
                    color: colorMode === "type" ? colorByType : "chainid",
                    opacity: 0.7,
                    smoothSheet: true
                }});
                
                // Add ball and stick for target residue
                component.addRepresentation("ball+stick", {{
                    sele: siteSelection,
                    color: colorMode === "type" ? colorByType : "element",
                    aspectRatio: 1.5,
                    scale: 1.2
                }});
                
                // Add licorice for environment (if not full view)
                if (!isFullView) {{
                    component.addRepresentation("licorice", {{
                        sele: environmentSelection + " and not " + siteSelection,
                        color: colorMode === "type" ? colorByType : "element",
                        opacity: 0.8,
                        scale: 0.8
                    }});
                    
                    // Add labels
                    component.addRepresentation("label", {{
                        sele: environmentSelection,
                        color: "#333333",
                        labelType: "format",
                        labelFormat: "{{resname}}{{resno}}",
                        labelGrouping: "residue",
                        attachment: "middle-center",
                        showBackground: true,
                        backgroundColor: "white",
                        backgroundOpacity: 0.5
                    }});
                }}
                
                // Set view
                if (isFullView) {{
                    component.autoView();
                }} else {{
                    component.autoView(environmentSelection, 2000);
                }}
            }}
            
            // Initial setup
            updateRepresentations();
        }}).catch(function(error) {{
            console.error("Error loading structure:", error);
            document.getElementById('site-viewer').innerHTML = 
                '<div class="alert alert-danger mt-3">Error loading structure: ' + error.message + '</div>';
        }});
    }});
    </script>
    """
    
    return js_code

def create_comparative_motif_visualization(primary_site, matches):
    """
    Create a comparative visualization of sequence motifs for the primary site
    and its structural matches, showing -5 to +5 range around the phosphosite.
    
    Args:
        primary_site: Dictionary with primary site information
        matches: List of dictionaries with match information
        
    Returns:
        HTML code for the visualization
    """
    if not primary_site or 'motif' not in primary_site:
        return "<div class='alert alert-warning'>Motif data not available for primary site</div>"

    # Get primary site motif
    primary_motif = primary_site.get('motif', '')
    primary_site_name = primary_site.get('site', 'Unknown')

    # IMPROVED: More aggressive approach to get the UniProt ID
    # First try the direct uniprot_id field
    primary_uniprot = primary_site.get('uniprot_id', '')

    # If that's empty, try other common field names
    if not primary_uniprot:
        for field in ['query_uniprot', 'protein_id', 'uniprotid', 'uniprot']:
            if field in primary_site and primary_site[field]:
                primary_uniprot = primary_site[field]
                break

    # If still empty, try to get from protein dictionary if it exists
    if not primary_uniprot and isinstance(primary_site.get('protein'), dict):
        primary_uniprot = primary_site['protein'].get('uniprot_id', '')

    # If still empty, check if it's in any matches (as query_uniprot)
    if not primary_uniprot and matches:
        for match in matches:
            if 'query_uniprot' in match and match['query_uniprot']:
                primary_uniprot = match['query_uniprot']
                break

    # If still empty, try to infer from site_id if present
    if not primary_uniprot and 'site_id' in primary_site:
        site_id = primary_site['site_id']
        if '_' in site_id:
            primary_uniprot = site_id.split('_')[0]

    # EXPLICIT DEBUG: Print what we found
    print(f"DEBUG - Primary UniProt: {primary_uniprot}")
    if not primary_uniprot:
        print("DEBUG - Failed to find UniProt ID in:", primary_site.keys())


    # Get position information for proper X padding
    primary_site_pos = None
    if primary_site_name:
        import re
        site_match = re.match(r'([STY])(\d+)', primary_site_name)
        if site_match:
            primary_site_pos = int(site_match.group(2))
    
    # Filter matches that have motif data
    valid_matches = [m for m in matches if 'motif' in m and m['motif']]
    
    if not valid_matches:
        return "<div class='alert alert-warning'>No motif data available for matches</div>"
    
    # Sort by RMSD (closest matches first)
    sorted_matches = sorted(valid_matches, key=lambda x: x.get('rmsd', float('inf')))
    
    # Take top N matches
    top_matches = sorted_matches[:10]  # Limit to 10 for visualization
    
    # Create HTML
    html = """
    <style>
        .motif-comparison {
            font-family: monospace;
            margin-bottom: 20px;
        }
        .motif-row {
            display: flex;
            align-items: center;
            margin-bottom: 5px;
        }
        .motif-label {
            width: 130px;
            font-weight: bold;
            text-align: right;
            padding-right: 10px;
            font-size: 0.9rem;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }
        .motif-sequence {
            display: flex;
        }
        .motif-aa {
            width: 24px;
            height: 24px;
            display: flex;
            align-items: center;
            justify-content: center;
            margin: 0 1px;
            border-radius: 3px;
        }
        .motif-aa.highlighted {
            background-color: #ff5722;
            color: white;
            font-weight: bold;
        }
        .motif-aa.sty {
            background-color: #bbdefb;
        }
        .motif-aa.nq {
            background-color: #b39ddb;
        }
        .motif-aa.cys {
            background-color: #ffcc80;
        }
        .motif-aa.proline {
            background-color: #81c784;
        }
        .motif-aa.nonpolar {
            background-color: #ffecb3;
        }
        .motif-aa.acidic {
            background-color: #ffcdd2;
        }
        .motif-aa.basic {
            background-color: #c8e6c9;
        }
        .motif-aa.special {
            background-color: #e1bee7;
        }
        .motif-aa.aa-x {
            background-color: #e0e0e0;
            color: #9e9e9e;
        }
        .match-info {
            margin-left: 10px;
            font-size: 12px;
            color: #333;
        }
        .motif-position {
            display: flex;
            padding-left: 130px;
            margin-bottom: 10px;
        }
        .motif-position span {
            width: 24px;
            text-align: center;
            font-size: 10px;
            color: #666;
        }
    </style>
    
    <div class="motif-comparison">
        <h5 class="mb-3">Motif Comparison</h5>
        
        <!-- Position markers - CHANGED to -5 to +5 -->
        <div class="motif-position">
    """
    
    # CHANGED: Add position markers from -5 to +5 instead of -7 to +7
    for i in range(-5, 6):
        html += f'<span>{i}</span>'
    
    html += """
        </div>
    """
    
    # Helper function to get amino acid class for coloring
    def get_aa_class(aa):
        if aa == 'X':
            return "aa-x"
        elif aa in 'STY':
            return "sty"
        elif aa in 'NQ':
            return "nq"
        elif aa == 'C':
            return "cys"
        elif aa == 'P':
            return "proline"
        elif aa in 'AVILMFWG':
            return "nonpolar"
        elif aa in 'DE':
            return "acidic"
        elif aa in 'KRH':
            return "basic"
        else:
            return "special"
    
    # Full standardization function for -7:+7 motifs
    def standardize_motif(motif, site_position=None):
        """
        Standardize a phosphosite motif to have exactly 7 positions before and after
        the phosphosite, with proper X padding.
        
        Args:
            motif (str): The motif sequence
            site_position (int, optional): The position of the site in the protein sequence
            
        Returns:
            str: The standardized motif
        """
        # Find the center position (phosphosite)
        center_pos = len(motif) // 2
        
        # Get the phosphosite and parts before/after
        site_char = motif[center_pos]
        before_site = motif[:center_pos]
        after_site = motif[center_pos + 1:]
        
        # If we have the absolute site position
        if site_position is not None:
            # Calculate padding needed at beginning based on site position
            aas_before = site_position - 1  # e.g., for S6, this is 5
            padding_needed = max(0, 7 - aas_before)
            
            # Create the before part: add X padding at beginning if needed
            if len(before_site) <= 7:
                # If we have 7 or fewer residues, use all of them with padding
                padded_before = "X" * padding_needed + before_site
            else:
                # If we have more than 7, take the last 7
                padded_before = before_site[-7:]
            
            # Create the after part: take exactly 7 chars, NO padding at end
            padded_after = after_site[:7]
        else:
            # Default behavior when we don't know the site position
            # Ensure we have exactly 7 characters before
            if len(before_site) < 7:
                padded_before = "X" * (7 - len(before_site)) + before_site
            else:
                padded_before = before_site[-7:]
            
            # Ensure we have exactly 7 characters after, no padding unless needed
            padded_after = after_site[:7]
        
        return padded_before + site_char + padded_after
    
    # NEW: Function to trim to just -5:+5 range
    def trim_to_central_range(motif_str):
        """Trim a standardized 15-char motif to just the central 11 positions (-5:+5)"""
        # Assuming the motif is standardized to 15 chars with the phosphosite at position 7 (0-indexed)
        # We want to keep positions 2-12 (0-indexed), which are -5 to +5 around the phosphosite
        return motif_str[2:13]
    
    # Modified helper function to create HTML for a motif
    def create_motif_html(motif, site_pos=None):
        # First standardize motif to full 15 chars (7+1+7)
        std_motif = standardize_motif(motif, site_pos)
        
        # Then trim to just -5:+5 range
        trimmed_motif = trim_to_central_range(std_motif)
        
        # Create HTML for each amino acid in the trimmed motif
        html = '<div class="motif-sequence" style="display: flex; flex-wrap: nowrap;">'
        for i, aa in enumerate(trimmed_motif):
            aa_class = get_aa_class(aa)
            # Center position (phosphosite) is now at position 5 (0-indexed) in the trimmed motif
            highlight_class = "highlighted" if i == 5 else aa_class
            html += f'<div class="motif-aa {highlight_class}" style="width: 24px; height: 24px; display: flex; align-items: center; justify-content: center; margin: 0 1px; border-radius: 3px;">{aa}</div>'
        html += '</div>'
        return html
    
    # Add primary site motif with UniProt ID
    html += f"""
        <div class="motif-row">
            <div class="motif-label">{primary_uniprot}_{primary_site_name}:</div>
            {create_motif_html(primary_motif, primary_site_pos)}
        </div>
    """

    # Add match motifs
    for match in top_matches:
        motif = match.get('motif', '')
        target_site = match.get('target_site', 'Unknown')
        target_uniprot = match.get('target_uniprot', 'Unknown')
        rmsd = match.get('rmsd', 0.0)
        
        # Extract site position from target_site if possible
        target_site_pos = None
        import re
        site_match = re.match(r'([STY])(\d+)', target_site)
        if site_match:
            target_site_pos = int(site_match.group(2))
        else:
            # Try another pattern like just digits
            site_match = re.match(r'(\d+)', target_site)
            if site_match:
                target_site_pos = int(site_match.group(1))
        
        html += f"""
        <div class="motif-row">
            <div class="motif-label">{target_uniprot}_{target_site}:</div>
            {create_motif_html(motif, target_site_pos)}
            <div class="match-info">
                RMSD: {rmsd:.2f}Å | <a href="/site/{target_uniprot}/{target_site}" class="text-decoration-none">View site</a>
            </div>
        </div>
        """
    
    html += """
    </div>
    """
    
    return html


def get_aa_bg_color(aa):
    """Get the background color for an amino acid based on its type."""
    if aa == 'X':
        return "#e0e0e0"  # Light gray for placeholder X
    elif aa in 'STY':
        return "#bbdefb"  # Light blue for STY
    elif aa in 'NQ':
        return "#b39ddb"  # Light purple for NQ
    elif aa == 'C':
        return "#ffcc80"  # Light orange for Cysteine
    elif aa == 'P':
        return "#81c784"  # Light green for Proline
    elif aa in 'AVILMFWG':
        return "#ffecb3"  # Light yellow for other nonpolar
    elif aa in 'DE':
        return "#ffcdd2"  # Light red for acidic
    elif aa in 'KRH':
        return "#c8e6c9"  # Pale green for basic
    else:
        return "#e1bee7"  # Light pink for special cases




def analyze_residue_distributions(structural_matches):
    """
    Analyze the distribution of residues across structural matches to identify
    potential conservation patterns.
    
    Args:
        structural_matches: List of match dictionaries
        
    Returns:
        Dictionary with analysis results
    """
    if not structural_matches:
        return None
        
    # Get motifs from matches
    motifs = []
    for match in structural_matches:
        if 'motif' in match and match['motif']:
            # We assume the phosphosite is in the middle of the motif
            motifs.append(match['motif'])
    
    if not motifs:
        return None
        
    # Determine the motif length (use longest motif)
    motif_length = max(len(m) for m in motifs)
    
    # Calculate the center position (where the phosphosite is)
    center_pos = motif_length // 2
    
    # Count amino acids at each position
    position_counts = []
    for i in range(motif_length):
        counts = {}
        for motif in motifs:
            if i < len(motif):
                aa = motif[i]
                counts[aa] = counts.get(aa, 0) + 1
        position_counts.append(counts)
    
    # Calculate frequencies and identify consensus
    frequencies = []
    consensus = []
    
    for i, counts in enumerate(position_counts):
        total = sum(counts.values())
        freq = {aa: count/total for aa, count in counts.items()}
        frequencies.append(freq)
        
        # Find most common AA
        if counts:
            max_aa = max(counts.items(), key=lambda x: x[1])
            consensus.append(max_aa[0])
        else:
            consensus.append('-')
    
    # Generate relative position labels
    positions = [i - center_pos for i in range(motif_length)]
    
    # Identify conserved positions (>50% same AA)
    conserved = []
    for i, counts in enumerate(position_counts):
        total = sum(counts.values())
        max_count = max(counts.values()) if counts else 0
        
        if max_count / total >= 0.5:
            conserved.append({
                'position': positions[i],
                'amino_acid': consensus[i],
                'frequency': max_count / total * 100
            })
    
    # Group amino acids by type
    aa_groups = {
        'polar': 'STYCNQ',
        'nonpolar': 'AVILMFWPG',
        'acidic': 'DE',
        'basic': 'KRH'
    }
    
    # Count by group at each position
    group_counts = []
    for i in range(motif_length):
        counts = {group: 0 for group in aa_groups}
        for motif in motifs:
            if i < len(motif):
                aa = motif[i]
                for group, aas in aa_groups.items():
                    if aa in aas:
                        counts[group] += 1
                        break
        group_counts.append(counts)
    
    # Generate consensus motif
    consensus_str = ''.join(consensus)
    
    return {
        'motif_count': len(motifs),
        'positions': positions,
        'consensus': consensus_str,
        'frequencies': frequencies,
        'position_counts': position_counts,
        'group_counts': group_counts,
        'conserved': conserved
    }
    
load_phosphosite_supp_data()



if __name__ == "__main__":
    # Example usage when run as a script
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python phospho_analyzer.py <uniprot_id or gene_symbol> [uniprot|gene]")
        sys.exit(1)
    
    identifier = sys.argv[1]
    id_type = sys.argv[2] if len(sys.argv) > 2 else 'uniprot'
    
    try:
        results = analyze_protein(identifier, id_type)
        
        # Print basic info
        print(f"\nProtein: {results['protein_info']['gene_symbol']} ({results['protein_info']['uniprot_id']})")
        print(f"Name: {results['protein_info']['name']}")
        print(f"Phosphosites found: {len(results['phosphosites'])}")
        
        if results['structural_matches']:
            match_count = sum(len(matches) for matches in results['structural_matches'].values())
            print(f"Structural matches found: {match_count}")
        else:
            print("No structural matches found")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)