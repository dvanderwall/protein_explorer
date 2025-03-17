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
                           parquet_file: str = None, top_n: int = 5) -> Dict[str, List[Dict]]:
    """
    Find structural matches for phosphosites in the kinome dataset.
    
    Args:
        uniprot_id: UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries from analyze_phosphosites
        parquet_file: Path to the parquet file with structural similarity data
        top_n: Number of top matches to return per site
        
    Returns:
        Dictionary mapping site IDs to lists of match dictionaries
    """
    logger.info(f"Finding structural matches for {uniprot_id}")
    
    # Check if parquet file exists, use default if not provided
    if parquet_file is None:
        # Try to find the parquet file in the parent directory
        current_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(current_dir)
        parquet_file = os.path.join(parent_dir, 'Combined_Kinome_10A_Master_Filtered_2.parquet')
    
    if not os.path.exists(parquet_file):
        logger.error(f"Parquet file not found: {parquet_file}")
        raise FileNotFoundError(f"Structural similarity data file not found: {parquet_file}")
    
    try:
        # Read parquet file
        logger.info(f"Reading parquet file: {parquet_file}")
        df = pd.read_parquet(parquet_file)
        
        # Create site IDs in the format UniprotID_ResNo
        site_ids = [f"{uniprot_id}_{site['resno']}" for site in phosphosites]
        
        # Find matches and sort by RMSD
        matches = []
        for site_id in site_ids:
            site_matches = df[df['Query'] == site_id].sort_values('RMSD').head(top_n)
            
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