"""
This script creates a supplementary data file for phosphosites with structural information.
The script processes data by UniProt ID to maximize efficiency, analyzing all sites
for a given protein at once after retrieving its structure and sequence.

Usage:
python precompute_structural_data.py [input_file] [output_file]

If no arguments are provided, it will look for Combined_Kinome_10A_Master_Filtered_2.feather
in the current directory and save the supplementary data as PhosphositeSuppData.feather
"""

import os
import sys
import time
import pandas as pd
import numpy as np
import logging
import requests
import re
from tqdm import tqdm
from Bio.PDB import PDBParser, NeighborSearch, Selection
import io
from collections import defaultdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("phosphosite_data.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Constants
CACHE_DIR = os.path.expanduser("~/.protein_explorer/cache")
BATCH_SIZE = 50  # Process proteins in batches

def ensure_cache_directory():
    """Ensure the cache directory exists."""
    global CACHE_DIR
    
    if not os.path.exists(CACHE_DIR):
        try:
            os.makedirs(CACHE_DIR, exist_ok=True)
            logger.info(f"Created cache directory: {CACHE_DIR}")
        except Exception as e:
            logger.error(f"Failed to create cache directory: {e}")
            # Fall back to temporary directory
            import tempfile
            CACHE_DIR = os.path.join(tempfile.gettempdir(), "protein_explorer_cache")
            os.makedirs(CACHE_DIR, exist_ok=True)
            logger.info(f"Using alternative cache directory: {CACHE_DIR}")

def get_alphafold_structure(uniprot_id):
    """Get AlphaFold structure for a given UniProt ID."""
    # Check cache first
    cache_filename = f"alphafold_{uniprot_id}.pdb"
    cache_file = os.path.join(CACHE_DIR, cache_filename)
    
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r') as f:
                return f.read()
        except Exception as e:
            logger.warning(f"Error reading from cache for {uniprot_id}: {e}")
    
    # Try different AlphaFold DB URLs
    urls_to_try = [
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb",
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v3.pdb",
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v2.pdb",
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v1.pdb"
    ]
    
    for url in urls_to_try:
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                structure = response.text
                # Cache the result
                try:
                    with open(cache_file, 'w') as f:
                        f.write(structure)
                except Exception as e:
                    logger.warning(f"Error caching structure for {uniprot_id}: {e}")
                
                return structure
        except requests.exceptions.RequestException as e:
            continue
    
    logger.warning(f"Failed to retrieve structure for {uniprot_id}")
    return None

def get_uniprot_sequence(uniprot_id):
    """Get protein sequence from UniProt."""
    # Check cache first
    cache_filename = f"uniprot_sequence_{uniprot_id}.txt"
    cache_file = os.path.join(CACHE_DIR, cache_filename)
    
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r') as f:
                return f.read()
        except Exception as e:
            logger.warning(f"Error reading sequence from cache for {uniprot_id}: {e}")
    
    # Try to get from UniProt REST API
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            lines = response.text.strip().split('\n')
            if len(lines) > 1:
                # Remove header line and join sequence lines
                sequence = ''.join(lines[1:])
                
                # Cache the result
                try:
                    with open(cache_file, 'w') as f:
                        f.write(sequence)
                except Exception as e:
                    logger.warning(f"Error caching sequence for {uniprot_id}: {e}")
                
                return sequence
    except requests.exceptions.RequestException as e:
        logger.warning(f"Error getting sequence for {uniprot_id}: {e}")
    
    logger.warning(f"Failed to retrieve sequence for {uniprot_id}")
    return None

def process_protein_sites(uniprot_id, site_numbers):
    """
    Process all phosphorylation sites for a single protein.
    This is much more efficient than processing each site individually.
    
    Args:
        uniprot_id: UniProt ID of the protein
        site_numbers: List of residue numbers to process
        
    Returns:
        List of dictionaries containing data for each site
    """
    try:
        # Get protein structure and sequence (only once)
        structure_data = get_alphafold_structure(uniprot_id)
        sequence = get_uniprot_sequence(uniprot_id)
        
        if not structure_data or not sequence:
            logger.warning(f"Missing structure or sequence for {uniprot_id}")
            return []
        
        # Parse the structure once (this is an expensive operation)
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", io.StringIO(structure_data))
        
        # Create NeighborSearch object once (also expensive)
        all_atoms = list(structure.get_atoms())
        ns = NeighborSearch(all_atoms)
        
        # Create a map of residue IDs to residue objects for quick lookup
        residue_map = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Store only standard residues (not hetero-residues)
                    if residue.get_id()[0] == " ":
                        residue_map[residue.get_id()[1]] = residue
        
        # Process each site
        results = []
        for resno in site_numbers:
            # Skip if residue number is invalid
            if resno <= 0 or resno > len(sequence):
                logger.warning(f"Residue number {resno} out of range for {uniprot_id}")
                continue
            
            site_type = sequence[resno-1]
            site_key = f"{uniprot_id}_{resno}"
            
            # Initialize data for this site
            site_data = {
                'site_id': site_key,
                'uniprot_id': uniprot_id,
                'residue_number': resno,
                'residue_type': site_type,
                'is_known_phosphosite': True  # All sites in Combined_Kinome are known phosphosites
            }
            
            # Only proceed if the residue exists in the structure
            if resno in residue_map:
                target_residue = residue_map[resno]
                
                # Get motif sequence
                motif = extract_motif(sequence, resno)
                site_data['motif'] = motif
                
                # Calculate AA frequencies if motif is available
                if motif:
                    aa_freq = calculate_aa_frequency(motif)
                    if aa_freq:
                        for group, percent in aa_freq.items():
                            site_data[f'{group}_aa_percent'] = percent
                
                # Calculate site pLDDT
                site_data['site_plddt'] = calculate_site_plddt(target_residue)
                
                # Calculate motif pLDDT
                site_data['motif_plddt'] = calculate_motif_plddt(
                    structure, sequence, resno, residue_map)
                
                # Get center atom for the residue
                if 'CA' in target_residue:
                    center_atom = target_residue['CA']
                else:
                    # Use the first atom if CA not available
                    center_atom = next(target_residue.get_atoms())
                
                # Calculate nearby residues
                nearby_atoms = ns.search(center_atom.get_coord(), 10)  # 10Å radius
                
                # Count unique residues (excluding the target residue)
                nearby_residues = set()
                for atom in nearby_atoms:
                    parent = atom.get_parent()
                    if parent != target_residue:
                        nearby_residues.add(parent)
                
                site_data['nearby_count'] = len(nearby_residues)
                
                # Calculate surface accessibility
                site_data['surface_accessibility'] = calculate_surface_accessibility(
                    len(nearby_residues))
            
            results.append(site_data)
        
        return results
    except Exception as e:
        logger.error(f"Error processing protein {uniprot_id}: {e}")
        return []

def extract_motif(sequence, resno):
    """Extract the -7 to +7 amino acid motif around the specified residue."""
    if not sequence:
        return None
    
    try:
        # Convert to 0-based index
        idx = resno - 1
        if idx < 0 or idx >= len(sequence):
            return None
        
        # Extract surrounding residues
        start = max(0, idx - 7)
        end = min(len(sequence), idx + 8)  # +8 because end index is exclusive
        
        return sequence[start:end]
    except Exception as e:
        logger.warning(f"Error extracting motif: {e}")
        return None

def calculate_site_plddt(residue):
    """Calculate the pLDDT score for a specific residue."""
    try:
        # Calculate average B-factor for this residue
        b_factors = [atom.get_bfactor() for atom in residue if not atom.get_name().startswith('H')]
        if b_factors:
            return float(np.mean(b_factors))
        
        return None
    except Exception as e:
        logger.warning(f"Error calculating site pLDDT: {e}")
        return None

def calculate_motif_plddt(structure, sequence, resno, residue_map):
    """Calculate the mean pLDDT for the motif around the specified residue."""
    try:
        # Extract motif sequence indices (7 residues before and after)
        sequence_length = len(sequence)
        motif_start = max(0, resno - 1 - 7)  # Convert to 0-based index and go back 7
        motif_end = min(sequence_length, resno - 1 + 8)  # +8 because range is exclusive
        
        # Map to PDB residue numbers (1-based)
        pdb_indices = list(range(motif_start + 1, motif_end + 1))
        
        # Extract B-factors (pLDDT scores in AlphaFold)
        plddt_values = []
        
        for idx in pdb_indices:
            if idx in residue_map:
                residue = residue_map[idx]
                b_factors = [atom.get_bfactor() for atom in residue if not atom.get_name().startswith('H')]
                if b_factors:
                    plddt_values.append(np.mean(b_factors))
        
        if plddt_values:
            return float(np.mean(plddt_values))
        
        return None
    except Exception as e:
        logger.warning(f"Error calculating motif pLDDT: {e}")
        return None

def calculate_surface_accessibility(nearby_count):
    """Calculate the relative surface accessibility based on neighbor count."""
    try:
        # Invert and normalize - fewer neighbors means more exposed
        max_neighbors = 50  # Approximate maximum reasonable number of neighbors
        return max(0, min(100, (max_neighbors - nearby_count) / max_neighbors * 100))
    except Exception as e:
        logger.warning(f"Error calculating surface accessibility: {e}")
        return None

def calculate_aa_frequency(motif):
    """Calculate amino acid frequencies in the motif."""
    if not motif:
        return None
    
    try:
        # Count amino acids by type
        aa_groups = {
            'polar': 'STYCNQ',
            'nonpolar': 'AVILMFWPG',
            'acidic': 'DE',
            'basic': 'KRH'
        }
        
        counts = {group: 0 for group in aa_groups}
        for aa in motif:
            for group, aas in aa_groups.items():
                if aa in aas:
                    counts[group] += 1
                    break
        
        # Convert to percentages
        motif_length = len(motif)
        percentages = {group: count/motif_length*100 for group, count in counts.items()}
        
        return percentages
    except Exception as e:
        logger.warning(f"Error calculating AA frequencies: {e}")
        return None

def group_sites_by_uniprot(site_list):
    """
    Group site IDs by UniProt ID for more efficient processing.
    
    Args:
        site_list: List of site IDs in the format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary mapping UniProt IDs to lists of residue numbers
    """
    grouped_sites = defaultdict(list)
    
    for site in site_list:
        match = re.match(r'([A-Z0-9_]+)_(\d+)', site)
        if match:
            uniprot_id, resno_str = match.groups()
            try:
                resno = int(resno_str)
                grouped_sites[uniprot_id].append(resno)
            except ValueError:
                logger.warning(f"Invalid residue number in site {site}")
    
    return grouped_sites

def main():
    # Parse command line arguments
    input_file = sys.argv[1] if len(sys.argv) > 1 else "C:/Users/mz30/protein_explorer/Combined_Kinome_10A_Master_Filtered_2.feather"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "C:/Users/mz30/protein_explorer/PhosphositeSuppData.feather"
    
    logger.info(f"Starting creation of phosphosite supplementary data from {input_file}")
    logger.info(f"Output will be saved to {output_file}")
    
    # Ensure cache directory exists
    ensure_cache_directory()
    
    # Load the data
    try:
        if input_file.endswith('.feather'):
            df = pd.read_feather(input_file)
        elif input_file.endswith('.parquet'):
            df = pd.read_parquet(input_file)
        else:
            raise ValueError("Input file must be .feather or .parquet format")
    except Exception as e:
        logger.error(f"Error loading input file: {e}")
        return
    
    logger.info(f"Loaded {len(df)} rows from {input_file}")
    
    # Extract unique site IDs (from both Query and Target columns)
    unique_sites = set(df['Query'].unique()).union(set(df['Target'].unique()))
    logger.info(f"Found {len(unique_sites)} unique phosphosites to process")
    
    # Group sites by UniProt ID
    grouped_sites = group_sites_by_uniprot(unique_sites)
    logger.info(f"Grouped into {len(grouped_sites)} unique proteins")
    
    # Create temp directory for intermediate results
    temp_dir = "phosphosite_data_temp"
    os.makedirs(temp_dir, exist_ok=True)
    
    # Process proteins in batches
    uniprot_ids = list(grouped_sites.keys())
    total_batches = (len(uniprot_ids) + BATCH_SIZE - 1) // BATCH_SIZE
    
    all_results = []
    
    # Process batches of proteins
    for i in range(0, len(uniprot_ids), BATCH_SIZE):
        batch_proteins = uniprot_ids[i:i+BATCH_SIZE]
        batch_number = i // BATCH_SIZE + 1
        
        logger.info(f"Processing batch {batch_number}/{total_batches} ({len(batch_proteins)} proteins)")
        
        # Process each protein in the batch
        batch_results = []
        with tqdm(total=len(batch_proteins), desc=f"Batch {batch_number}/{total_batches}") as progress:
            for uniprot_id in batch_proteins:
                site_numbers = grouped_sites[uniprot_id]
                protein_results = process_protein_sites(uniprot_id, site_numbers)
                
                # Add all site results for this protein
                batch_results.extend(protein_results)
                progress.update(1)
                
                # Brief pause to prevent rate limiting
                time.sleep(0.01)
        
        # Save intermediate results
        if batch_results:
            batch_df = pd.DataFrame(batch_results)
            temp_file = os.path.join(temp_dir, f"batch_{i}.feather")
            batch_df.to_feather(temp_file)
            logger.info(f"Saved intermediate batch {batch_number} with {len(batch_results)} processed sites to {temp_file}")
            
            all_results.extend(batch_results)
    
    # Combine all results
    logger.info("Combining all results...")
    
    if all_results:
        final_df = pd.DataFrame(all_results)
    else:
        # Read from temp files if memory was an issue
        result_dfs = []
        for i in range(0, len(uniprot_ids), BATCH_SIZE):
            temp_file = os.path.join(temp_dir, f"batch_{i}.feather")
            if os.path.exists(temp_file):
                batch_df = pd.read_feather(temp_file)
                result_dfs.append(batch_df)
        
        if result_dfs:
            final_df = pd.concat(result_dfs, ignore_index=True)
        else:
            logger.error("No results found!")
            return
    
    # Save the final result
    logger.info(f"Saving final result to {output_file}")
    final_df.to_feather(output_file)
    
    # Clean up temp files
    logger.info("Cleaning up temporary files")
    for f in os.listdir(temp_dir):
        os.remove(os.path.join(temp_dir, f))
    os.rmdir(temp_dir)
    
    logger.info(f"Phosphosite supplementary data creation complete. Processed {len(final_df)} sites.")

if __name__ == "__main__":
    start_time = time.time()
    main()
    elapsed_time = time.time() - start_time
    logger.info(f"Total execution time: {elapsed_time:.2f} seconds")