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
# Global variable to store the loaded structural similarity data
STRUCTURAL_SIMILARITY_DF = None

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

def get_phosphosites(uniprot_id: str) -> List[Dict]:
    """
    Analyze potential phosphorylation sites for a protein.
    
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
        
        return phosphosites
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
            # Get phosphosites
            phosphosites = get_phosphosites(uniprot_id)
            
            # Find structural matches
            try:
                structural_matches = find_structural_matches(uniprot_id, phosphosites, parquet_file)
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