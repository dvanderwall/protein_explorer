"""
Functions for analyzing sequence similarity between phosphorylation sites.

This module handles loading sequence similarity data from parquet files,
performing sequence-based analysis, and providing visualization-ready outputs.
It integrates with supplementary data to include motif information.
"""

import os
import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Tuple, Optional, Union
from collections import Counter, defaultdict
import re

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Global variables to store loaded data
SEQUENCE_SIMILARITY_DF = None
SEQ_SIMILARITY_CACHE = {}

def preload_sequence_data(file_path: str = None) -> None:
    """
    Preload sequence similarity data at application startup.
    
    Args:
        file_path: Path to the sequence similarity data file (parquet format)
    """
    global SEQUENCE_SIMILARITY_DF
    
    if SEQUENCE_SIMILARITY_DF is not None:
        logger.info("Sequence similarity data already loaded")
        return
    
    # Find the data file
    if file_path is None:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(os.path.dirname(current_dir))
        file_path = os.path.join(parent_dir, 'Sequence_Similarity_Edges.parquet')
    
    try:
        # Load the data
        logger.info(f"Preloading sequence similarity data from: {file_path}")
        SEQUENCE_SIMILARITY_DF = pd.read_parquet(file_path)
        
        # Create indices for faster querying
        logger.info("Creating query indices for sequence data")
        if 'ID1' in SEQUENCE_SIMILARITY_DF.columns:
            SEQUENCE_SIMILARITY_DF.set_index('ID1', drop=False, inplace=True)
        
        logger.info(f"Successfully preloaded {len(SEQUENCE_SIMILARITY_DF)} sequence similarity records")
    except Exception as e:
        logger.error(f"Error preloading sequence data: {e}")
        logger.warning("Will attempt to load data on first request")

# Run preloading at module import time if not already run
if __name__ == "__main__":
    # Example usage when run as a script
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python sequence_analyzer.py <site_id> [min_similarity]")
        sys.exit(1)
    
    site_id = sys.argv[1]
    min_similarity = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4
    
    try:
        # Ensure data is loaded
        if SEQUENCE_SIMILARITY_DF is None:
            load_sequence_similarity_data()
        
        # Find matches
        matches = find_sequence_matches(site_id, min_similarity=min_similarity)
        print(f"Found {len(matches)} sequence matches for {site_id}")
        
        # Print first few matches
        for i, match in enumerate(matches[:5]):
            print(f"Match {i+1}: {match['target_id']} (Similarity: {match['similarity']:.2f})")
            if 'motif' in match and match['motif']:
                print(f"  Motif: {match['motif']}")
            print()
            
        # Simple conservation analysis
        try:
            from protein_explorer.analysis.phospho_analyzer import get_phosphosite_data
            # Get motif for query site
            query_motif = None
            site_data = get_phosphosite_data(site_id)
            if site_data and 'motif' in site_data:
                query_motif = site_data['motif']
                print(f"Query motif: {query_motif}")
            
            # Analyze conservation
            conservation = analyze_motif_conservation(matches, query_motif)
            if conservation['motif_count'] > 0:
                print(f"Conservation analysis for {conservation['motif_count']} motifs:")
                print(f"Consensus: {conservation['consensus_motif']}")
                if conservation['conserved_positions']:
                    print("Conserved positions:")
                    for pos in conservation['conserved_positions']:
                        print(f"  {pos['position']}: {pos['amino_acid']} ({pos['frequency']:.1f}%)")
        except Exception as e:
            print(f"Error in conservation analysis: {e}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

def load_sequence_similarity_data(file_path: str = None) -> pd.DataFrame:
    """
    Load sequence similarity data from parquet file.
    
    Args:
        file_path: Path to the parquet file with sequence similarity data
        
    Returns:
        Pandas DataFrame with sequence similarity data
    """
    global SEQUENCE_SIMILARITY_DF
    
    if SEQUENCE_SIMILARITY_DF is not None:
        logger.info("Using cached sequence similarity data")
        return SEQUENCE_SIMILARITY_DF
    
    # Check if parquet file exists, use default if not provided
    if file_path is None:
        # Try to find the parquet file in the parent directory
        current_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(os.path.dirname(current_dir))
        file_path = os.path.join(parent_dir, 'Sequence_Similarity_Edges.parquet')
    
    if not os.path.exists(file_path):
        logger.error(f"Sequence similarity file not found: {file_path}")
        raise FileNotFoundError(f"Sequence similarity data file not found: {file_path}")
    
    try:
        # Read parquet file
        logger.info(f"Reading sequence similarity file: {file_path}")
        SEQUENCE_SIMILARITY_DF = pd.read_parquet(file_path)
        
        # Create indexes for faster querying
        logger.info("Creating query index for sequence data")
        if 'ID1' in SEQUENCE_SIMILARITY_DF.columns:
            SEQUENCE_SIMILARITY_DF.set_index('ID1', drop=False, inplace=True)
        
        logger.info(f"Loaded {len(SEQUENCE_SIMILARITY_DF)} sequence similarity records")
        return SEQUENCE_SIMILARITY_DF
    except Exception as e:
        logger.error(f"Error reading sequence similarity file: {e}")
        raise ValueError(f"Error reading sequence similarity data: {e}")

def find_sequence_matches(site_id: str, 
                         top_n: int = 200, 
                         min_similarity: float = 0.4) -> List[Dict]:
    """
    Find sequence similarity matches for a site.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        top_n: Maximum number of results to return
        min_similarity: Minimum similarity score to include (0-1)
        
    Returns:
        List of dictionaries with match information
    """
    global SEQUENCE_SIMILARITY_DF, SEQ_SIMILARITY_CACHE
    
    # Check cache first
    cache_key = f"{site_id}_{top_n}_{min_similarity}"
    if cache_key in SEQ_SIMILARITY_CACHE:
        logger.info(f"Using cached sequence matches for {site_id}")
        return SEQ_SIMILARITY_CACHE[cache_key]
    
    try:
        # Load data if not already loaded
        if SEQUENCE_SIMILARITY_DF is None:
            load_sequence_similarity_data()
        
        # Check if data was loaded successfully
        if SEQUENCE_SIMILARITY_DF is None:
            logger.error("Sequence similarity data not available")
            return []
        
        # Query for matches where the site is in ID1
        df = SEQUENCE_SIMILARITY_DF
        
        # Handle the case whether we have an index or not
        if isinstance(df.index, pd.MultiIndex) or df.index.name == 'ID1':
            # Use efficient index lookup
            if site_id in df.index:
                matches_df = df.loc[[site_id]]
            else:
                logger.warning(f"Site {site_id} not found in sequence similarity index")
                matches_df = pd.DataFrame()
        else:
            # Fallback to regular query
            matches_df = df[df['ID1'] == site_id]
        
        # If empty, check for matches where the site is in ID2
        if matches_df.empty:
            matches_df = df[df['ID2'] == site_id]
            
            # If still empty, return empty list
            if matches_df.empty:
                logger.warning(f"No sequence matches found for {site_id}")
                return []
            
            # Swap columns for consistency
            matches_df = matches_df.rename(columns={'ID1': 'ID2', 'ID2': 'ID1'})
            matches_df['ID1'] = site_id
        
        # Filter by minimum similarity
        matches_df = matches_df[matches_df['Similarity'] >= min_similarity]
        
        # Sort by similarity (highest first)
        matches_df = matches_df.sort_values('Similarity', ascending=False)
        
        # Take top N matches
        if top_n is not None:
            matches_df = matches_df.head(top_n)
        
        # Convert to list of dictionaries
        matches = []
        for _, row in matches_df.iterrows():
            target_id = row['ID2']
            
            # Parse target_id to extract UniProt ID and site number
            target_parts = target_id.split('_')
            if len(target_parts) >= 2:
                target_uniprot = target_parts[0]
                target_site = target_parts[1]
                
                # Extract site type (S/T/Y) if the information is available
                site_type = None
                site_match = re.match(r'([STY])(\d+)', target_site)
                if site_match:
                    site_type = site_match.group(1)
                    target_site = site_match.group(0)  # Get the full match
                else:
                    # Try another pattern like just digits
                    site_match = re.match(r'(\d+)', target_site)
                    if site_match:
                        site_type = None  # We don't know the type
                        target_site = site_match.group(0)
                    
                matches.append({
                    'query_id': site_id,
                    'target_id': target_id,
                    'target_uniprot': target_uniprot,
                    'target_site': target_site,
                    'site_type': site_type,
                    'similarity': float(row['Similarity'])
                })
        
        # Enhance matches with motif data from supplementary data
        try:
            # Import here to avoid circular imports
            from protein_explorer.analysis.phospho_analyzer import get_phosphosite_data, enhance_phosphosite
            
            enhanced_matches = []
            for match in matches:
                target_id = match['target_id']
                
                # Get supplementary data for this target site
                try:
                    target_supp = get_phosphosite_data(target_id)
                    
                    # Create enhanced match
                    enhanced_match = match.copy()
                    
                    # Add motif if available in supplementary data
                    if target_supp and 'motif' in target_supp and target_supp['motif'] is not None:
                        enhanced_match['motif'] = target_supp['motif']
                        logger.debug(f"Added motif for {target_id}: {target_supp['motif']}")
                    
                    # Add other useful supplementary data
                    for key in ['site_plddt', 'surface_accessibility', 'secondary_structure', 'nearby_count']:
                        if target_supp and key in target_supp and target_supp[key] is not None:
                            enhanced_match[key] = target_supp[key]
                    
                    enhanced_matches.append(enhanced_match)
                except Exception as e:
                    logger.error(f"Error enhancing match for {target_id}: {e}")
                    # Still include the original match if enhancement fails
                    enhanced_matches.append(match)
            
            # Use enhanced matches with supplementary data
            matches = enhanced_matches
            logger.info(f"Enhanced {len(enhanced_matches)} sequence matches with supplementary data")
        except Exception as e:
            logger.error(f"Error integrating with supplementary data: {e}")
            logger.warning("Using matches without supplementary data")
        
        # Cache results
        SEQ_SIMILARITY_CACHE[cache_key] = matches
        
        return matches
    except Exception as e:
        logger.error(f"Error finding sequence matches: {e}")
        return []

def extract_motif_from_site_id(site_id: str, motif_db: Dict = None) -> Optional[str]:
    """
    Extract motif sequence for a given site ID using a motif database.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        motif_db: Dictionary mapping site IDs to motifs
        
    Returns:
        Motif sequence or None if not found
    """
    # Check if we have a motif database to query
    if motif_db and site_id in motif_db:
        return motif_db[site_id]
    
    # If no motif database, try to get from supplementary data
    try:
        from protein_explorer.analysis.phospho_analyzer import get_phosphosite_data
        supp_data = get_phosphosite_data(site_id)
        if supp_data and 'motif' in supp_data and supp_data['motif'] is not None:
            return supp_data['motif']
    except Exception as e:
        logger.error(f"Error getting motif from supplementary data: {e}")
    
    # If no motif database or ID not found, return None
    return None

def enhance_sequence_matches(matches: List[Dict], 
                           motif_db: Dict = None,
                           enhanced_db: Dict = None) -> List[Dict]:
    """
    Enhance sequence matches with additional information.
    
    Args:
        matches: List of match dictionaries from find_sequence_matches
        motif_db: Dictionary mapping site IDs to motifs
        enhanced_db: Dictionary mapping site IDs to additional metadata
        
    Returns:
        Enhanced list of match dictionaries
    """
    # If no matches, return empty list
    if not matches:
        return []
    
    enhanced_matches = []
    
    for match in matches:
        # Create enhanced match from original
        enhanced = match.copy()
        
        # Add motif if available
        target_id = match['target_id']
        
        # First check if motif is already in the match
        if 'motif' not in enhanced or not enhanced['motif']:
            # If not, try to get from motif_db
            if motif_db and target_id in motif_db:
                enhanced['motif'] = motif_db[target_id]
            # If still not found, try supplementary data
            elif 'motif' not in enhanced or not enhanced['motif']:
                try:
                    motif = extract_motif_from_site_id(target_id)
                    if motif:
                        enhanced['motif'] = motif
                except Exception as e:
                    logger.error(f"Error extracting motif for {target_id}: {e}")
        
        # Add additional metadata if available
        if enhanced_db and target_id in enhanced_db:
            for key, value in enhanced_db[target_id].items():
                # Avoid overwriting existing keys
                if key not in enhanced:
                    enhanced[key] = value
        
        enhanced_matches.append(enhanced)
    
    return enhanced_matches

def analyze_motif_conservation(matches: List[Dict], 
                             query_motif: str = None,
                             motif_db: Dict = None) -> Dict:
    """
    Analyze conservation patterns in motifs of sequence-similar sites.
    
    Args:
        matches: List of match dictionaries from find_sequence_matches
        query_motif: Motif of the query site
        motif_db: Dictionary mapping site IDs to motifs
        
    Returns:
        Dictionary with conservation analysis results
    """
    logger.info(f"Analyzing conservation with {len(matches)} matches")
    # If no matches or no motifs available, return empty results
    if not matches or (not query_motif and not motif_db):
        return {
            'motif_count': 0,
            'conserved_positions': [],
            'n_term_analysis': {},
            'c_term_analysis': {},
            'position_frequencies': {}
        }
    
    # Collect motifs from matches
    motifs = []
    for match in matches:
        motif = None
        if 'motif' in match and match['motif']:
            motif = match['motif']
            motifs.append(motif)
            logger.info(f"Found motif in match: {motif}")
        elif motif_db and match['target_id'] in motif_db:
            motif = motif_db[match['target_id']]
            motifs.append(motif)
            logger.info(f"Found motif in db: {motif}")
        else:
            # Try to get motif from supplementary data
            try:
                target_id = match['target_id']
                motif = extract_motif_from_site_id(target_id)
                if motif:
                    motifs.append(motif)
                    logger.info(f"Found motif in supplementary data: {motif}")
            except Exception as e:
                logger.error(f"Error getting motif for {match.get('target_id', 'unknown')}: {e}")
    
    # Add query motif if provided
    if query_motif:
        motifs = [query_motif] + motifs
        logger.info(f"Added query motif, total motifs: {len(motifs)}")
    
    # If no motifs found, return empty results
    if not motifs:
        return {
            'motif_count': 0,
            'conserved_positions': [],
            'n_term_analysis': {},
            'c_term_analysis': {},
            'position_frequencies': {}
        }
    
    # Standardize motif length by assuming phosphosite is in the middle
    # and padding with X if necessary
    std_motifs = []
    for motif in motifs:
        center_pos = len(motif) // 2
        site_char = motif[center_pos]
        
        before_site = motif[:center_pos]
        after_site = motif[center_pos+1:]
        
        # Ensure we have 7 positions before and after
        if len(before_site) < 7:
            before_site = 'X' * (7 - len(before_site)) + before_site
        else:
            before_site = before_site[-7:]
            
        if len(after_site) < 7:
            after_site = after_site + 'X' * (7 - len(after_site))
        else:
            after_site = after_site[:7]
            
        std_motifs.append(before_site + site_char + after_site)
    
    # Analyze conservation at each position
    position_counts = []
    for i in range(15):  # -7 to +7 positions
        counts = Counter()
        for motif in std_motifs:
            if i < len(motif):
                counts[motif[i]] += 1
        position_counts.append(counts)
    
    # Calculate frequency of each amino acid at each position
    position_frequencies = {}
    motif_count = len(std_motifs)
    
    for i, counts in enumerate(position_counts):
        position = i - 7  # Convert to -7 to +7 positions
        position_frequencies[position] = {}
        
        for aa, count in counts.items():
            position_frequencies[position][aa] = count / motif_count
    
    # Identify positions with strong conservation (>50% same AA)
    conserved_positions = []
    for pos, freqs in position_frequencies.items():
        if pos == 0:  # Skip the phosphosite itself
            continue
            
        most_common = max(freqs.items(), key=lambda x: x[1], default=(None, 0))
        if most_common[1] >= 0.5:  # 50% or more conservation
            conserved_positions.append({
                'position': pos,
                'amino_acid': most_common[0],
                'frequency': most_common[1] * 100
            })
    
    # Analyze N-terminal and C-terminal regions separately
    n_term_motifs = [m[:7] for m in std_motifs]  # -7 to -1 positions
    c_term_motifs = [m[8:] for m in std_motifs]  # +1 to +7 positions
    
    # Amino acid group classification
    aa_groups = {
        'polar': 'STYCNQ',
        'nonpolar': 'AVILMFWPG',
        'acidic': 'DE',
        'basic': 'KRH',
        'other': 'X'
    }
    
    # Function to analyze region
    def analyze_region(motifs):
        aa_composition = defaultdict(int)
        aa_group_composition = defaultdict(int)
        total_aa = len(motifs) * 7  # 7 positions per motif
        
        for motif in motifs:
            for aa in motif:
                aa_composition[aa] += 1
                
                # Classify by group
                for group, members in aa_groups.items():
                    if aa in members:
                        aa_group_composition[group] += 1
                        break
        
        # Convert to percentages
        aa_percentages = {aa: count/total_aa*100 for aa, count in aa_composition.items()}
        group_percentages = {group: count/total_aa*100 for group, count in aa_group_composition.items()}
        
        return {
            'aa_composition': dict(sorted(aa_percentages.items(), key=lambda x: x[1], reverse=True)),
            'group_composition': dict(sorted(group_percentages.items(), key=lambda x: x[1], reverse=True))
        }
    
    n_term_analysis = analyze_region(n_term_motifs)
    c_term_analysis = analyze_region(c_term_motifs)
    
    # Calculate a consensus motif
    consensus_motif = ""
    for i in range(15):  # -7 to +7 positions
        if i == 7:  # Phosphosite position
            consensus_motif += std_motifs[0][7] if std_motifs else "X"
            continue
            
        counts = position_counts[i]
        if counts:
            # Get the most common AA, but exclude X unless it's the only one
            filtered_counts = {aa: count for aa, count in counts.items() if aa != 'X'}
            if filtered_counts:
                most_common = max(filtered_counts.items(), key=lambda x: x[1])[0]
            else:
                most_common = 'X'
            consensus_motif += most_common
        else:
            consensus_motif += "X"
    
    return {
        'motif_count': len(motifs),
        'consensus_motif': consensus_motif,
        'conserved_positions': conserved_positions,
        'n_term_analysis': n_term_analysis,
        'c_term_analysis': c_term_analysis,
        'position_frequencies': position_frequencies
    }

def create_sequence_network_data(query_site_id: str, 
                               matches: List[Dict],
                               query_motif: str = None) -> Dict:
    """
    Create data for sequence similarity network visualization.
    
    Args:
        query_site_id: Site ID of the query
        matches: List of match dictionaries from find_sequence_matches
        query_motif: Motif of the query site
        
    Returns:
        Dictionary with nodes and links for network visualization
    """
    if not matches:
        return {'nodes': [], 'links': []}
    
    # Extract parts of query_site_id
    query_parts = query_site_id.split('_')
    query_uniprot = query_parts[0] if len(query_parts) > 0 else ""
    query_site = query_parts[1] if len(query_parts) > 1 else ""
    
    # Extract site type if possible
    query_site_type = None
    site_match = re.match(r'([STY])(\d+)', query_site)
    if site_match:
        query_site_type = site_match.group(1)
        query_site_number = site_match.group(2)
        query_site = f"{query_site_type}{query_site_number}"
    
    # Create nodes list starting with query node
    nodes = [{
        'id': query_site_id,
        'name': query_site,
        'display_name': query_site_id,
        'uniprot': query_uniprot,
        'type': 'query',
        'site_type': query_site_type,
        'size': 12,
        'motif': query_motif
    }]
    
    # Create links list
    links = []
    
    # Process each match
    seen_nodes = {query_site_id}  # Track nodes we've already added
    
    for match in matches:
        target_id = match['target_id']
        
        # Skip self-matches
        if target_id == query_site_id:
            continue
            
        # Skip duplicates
        if target_id in seen_nodes:
            continue
            
        seen_nodes.add(target_id)
        
        # Create node for this match
        nodes.append({
            'id': target_id,
            'name': match['target_site'],
            'display_name': target_id,
            'uniprot': match['target_uniprot'],
            'type': 'target',
            'site_type': match.get('site_type'),
            'similarity': match['similarity'],
            'size': 8,
            'motif': match.get('motif')
        })
        
        # Create link between query and target
        links.append({
            'source': query_site_id,
            'target': target_id,
            'similarity': match['similarity']
        })
    
    return {
        'nodes': nodes,
        'links': links
    }

def get_motif_enrichment(matches: List[Dict], 
                       query_motif: str = None,
                       background_frequencies: Dict = None) -> Dict:
    """
    Calculate amino acid enrichment in motifs compared to background frequencies.
    
    Args:
        matches: List of match dictionaries
        query_motif: Motif of the query site
        background_frequencies: Background AA frequencies (defaults to UniProt averages)
        
    Returns:
        Dictionary with enrichment analysis results
    """
    # Default background frequencies from UniProt
    if background_frequencies is None:
        background_frequencies = {
            'A': 0.0825, 'R': 0.0553, 'N': 0.0406, 'D': 0.0545,
            'C': 0.0137, 'Q': 0.0393, 'E': 0.0675, 'G': 0.0707,
            'H': 0.0227, 'I': 0.0595, 'L': 0.0965, 'K': 0.0584,
            'M': 0.0241, 'F': 0.0386, 'P': 0.0470, 'S': 0.0656,
            'T': 0.0534, 'W': 0.0108, 'Y': 0.0292, 'V': 0.0687
        }
    
    # Collect motifs
    motifs = []
    for match in matches:
        if 'motif' in match and match['motif']:
            motifs.append(match['motif'])
    
    # Add query motif if provided
    if query_motif:
        motifs = [query_motif] + motifs
    
    # If no motifs, return empty results
    if not motifs:
        return {
            'position_enrichment': {},
            'overall_enrichment': {}
        }
    
    # Standardize motifs
    std_motifs = []
    for motif in motifs:
        center_pos = len(motif) // 2
        site_char = motif[center_pos]
        
        before_site = motif[:center_pos]
        after_site = motif[center_pos+1:]
        
        # Ensure we have 7 positions before and after
        if len(before_site) < 7:
            before_site = 'X' * (7 - len(before_site)) + before_site
        else:
            before_site = before_site[-7:]
            
        if len(after_site) < 7:
            after_site = after_site + 'X' * (7 - len(after_site))
        else:
            after_site = after_site[:7]
            
        std_motifs.append(before_site + site_char + after_site)
    
    # Count AAs at each position
    position_counts = []
    for i in range(15):  # -7 to +7 positions
        counts = Counter()
        for motif in std_motifs:
            if i < len(motif) and motif[i] != 'X':  # Skip placeholder X
                counts[motif[i]] += 1
        position_counts.append(counts)
    
    # Count overall AA frequencies (excluding phosphosite and X)
    overall_counts = Counter()
    for motif in std_motifs:
        for i, aa in enumerate(motif):
            if i != 7 and aa != 'X':  # Skip phosphosite and placeholder X
                overall_counts[aa] += 1
    
    # Calculate position-specific enrichment
    position_enrichment = {}
    for i, counts in enumerate(position_counts):
        position = i - 7  # Convert to -7 to +7 positions
        
        if position == 0:  # Skip phosphosite position
            continue
            
        position_enrichment[position] = {}
        total_aas = sum(counts.values())
        
        if total_aas == 0:
            continue
            
        for aa, count in counts.items():
            if aa == 'X' or aa not in background_frequencies:
                continue
                
            observed_freq = count / total_aas
            expected_freq = background_frequencies[aa]
            enrichment = observed_freq / expected_freq if expected_freq > 0 else 0
            
            position_enrichment[position][aa] = {
                'observed_freq': observed_freq,
                'expected_freq': expected_freq,
                'enrichment': enrichment,
                'count': count
            }
    
    # Calculate overall enrichment
    overall_enrichment = {}
    total_aas = sum(overall_counts.values())
    
    if total_aas > 0:
        for aa, count in overall_counts.items():
            if aa not in background_frequencies:
                continue
                
            observed_freq = count / total_aas
            expected_freq = background_frequencies[aa]
            enrichment = observed_freq / expected_freq if expected_freq > 0 else 0
            
            overall_enrichment[aa] = {
                'observed_freq': observed_freq,
                'expected_freq': expected_freq,
                'enrichment': enrichment,
                'count': count
            }
    
    return {
        'position_enrichment': position_enrichment,
        'overall_enrichment': overall_enrichment
    }

def create_sequence_motif_visualization(query_site_id: str, 
                                         query_motif: str,
                                         matches: List[Dict],
                                         max_matches: int = 10) -> str:
    """
    Create HTML for a comparative visualization of sequence motifs.
    
    Args:
        query_site_id: ID of the query site
        query_motif: Motif of the query site
        matches: List of match dictionaries
        max_matches: Maximum number of matches to display
        
    Returns:
        HTML code for the visualization
    """
    logger.info(f"Creating motif comparison with {len(matches)} matches")
    logger.info(f"Query motif: {query_motif}")

    # Extract uniprot and site from query_site_id
    query_parts = query_site_id.split('_')
    query_uniprot = query_parts[0] if len(query_parts) > 0 else ""
    query_site = query_parts[1] if len(query_parts) > 1 else ""
    
    site_match = re.match(r'([STY])(\d+)', query_site)
    if site_match:
        query_site = site_match.group(0)
    
    # If the query motif is empty or None, try to get it from supplementary data
    if not query_motif:
        try:
            from protein_explorer.analysis.phospho_analyzer import get_phosphosite_data
            supp_data = get_phosphosite_data(query_site_id)
            if supp_data and 'motif' in supp_data and supp_data['motif'] is not None:
                query_motif = supp_data['motif']
                logger.info(f"Retrieved query motif from supplementary data: {query_motif}")
        except Exception as e:
            logger.error(f"Error retrieving query motif from supplementary data: {e}")
    
    # Filter matches with motifs
    valid_matches = []
    for match in matches:
        # Check if match already has a motif
        if 'motif' in match and match['motif']:
            valid_matches.append(match)
            continue
            
        # If not, try to retrieve from supplementary data
        try:
            from protein_explorer.analysis.phospho_analyzer import get_phosphosite_data
            target_id = match['target_id']
            supp_data = get_phosphosite_data(target_id)
            if supp_data and 'motif' in supp_data and supp_data['motif'] is not None:
                # Create a new match with the motif added
                enhanced_match = match.copy()
                enhanced_match['motif'] = supp_data['motif']
                valid_matches.append(enhanced_match)
                logger.debug(f"Added motif for {target_id} from supplementary data")
            else:
                logger.debug(f"No motif found for {target_id} in supplementary data")
        except Exception as e:
            logger.error(f"Error retrieving motif for {match.get('target_id', 'unknown')}: {e}")
    
    logger.info(f"Found {len(valid_matches)} matches with motifs after enhancement")

    # If no valid matches, return simple message
    if not valid_matches:
        return f"""
        <div class="alert alert-info">
            No sequence motif data available for comparison with {query_site_id}.
        </div>
        """
    
    # Sort by similarity (highest first)
    sorted_matches = sorted(valid_matches, key=lambda x: x.get('similarity', 0), reverse=True)
    
    # Take top N matches
    top_matches = sorted_matches[:max_matches]
    
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
    
    # Standardization function for motifs
    def standardize_motif(motif):
        # Find the center position (phosphosite)
        center_pos = len(motif) // 2
        
        # Get the phosphosite and parts before/after
        site_char = motif[center_pos]
        before_site = motif[:center_pos]
        after_site = motif[center_pos + 1:]
        
        # Ensure we have exactly 7 characters before and after
        if len(before_site) < 7:
            before_site = "X" * (7 - len(before_site)) + before_site
        else:
            before_site = before_site[-7:]
            
        if len(after_site) < 7:
            after_site = after_site + "X" * (7 - len(after_site))
        else:
            after_site = after_site[:7]
            
        return before_site + site_char + after_site
    
    # Function to trim to just -5 to +5 range for display
    def trim_to_central_range(motif_str):
        # We want to keep positions 2-12 (0-indexed) from a 15-char motif
        # which represent positions -5 to +5 around the phosphosite
        return motif_str[2:13]
    
    # Modified helper function to create HTML for a motif
    def create_motif_html(motif):
        # First standardize to full 15 chars
        std_motif = standardize_motif(motif)
        
        # Then trim to -5 to +5 range
        trimmed_motif = trim_to_central_range(std_motif)
        
        # Create HTML
        html = '<div class="motif-sequence" style="display: flex; flex-wrap: nowrap;">'
        for i, aa in enumerate(trimmed_motif):
            aa_class = get_aa_class(aa)
            highlight_class = "highlighted" if i == 5 else aa_class  # Position 5 is the phosphosite in trimmed motif
            html += f'<div class="motif-aa {highlight_class}" style="width: 24px; height: 24px; display: flex; align-items: center; justify-content: center; margin: 0 1px; border-radius: 3px;">{aa}</div>'
        html += '</div>'
        return html
    
    # Create HTML for the visualization
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
        <h5 class="mb-3">Sequence Similarity Motif Comparison</h5>
        
        <!-- Position markers -->
        <div class="motif-position">
    """
    
    # Add position markers from -5 to +5
    for i in range(-5, 6):
        html += f'<span>{i}</span>'
    
    html += """
        </div>
    """
    
    # Add query motif row
    html += f"""
        <div class="motif-row">
            <div class="motif-label">{query_uniprot}_{query_site}:</div>
            {create_motif_html(query_motif)}
            <div class="match-info">
                Query
            </div>
        </div>
    """
    
    # Add match motifs
    for match in top_matches:
        motif = match.get('motif', '')
        target_site = match.get('target_site', 'Unknown')
        target_uniprot = match.get('target_uniprot', 'Unknown')
        similarity = match.get('similarity', 0.0)
        
        html += f"""
        <div class="motif-row">
            <div class="motif-label">{target_uniprot}_{target_site}:</div>
            {create_motif_html(motif)}
            <div class="match-info">
                Similarity: {similarity:.2f} | <a href="/site/{target_uniprot}/{target_site}" class="text-decoration-none">View site</a>
            </div>
        </div>
        """
    
    html += """
    </div>
    """
    
    return html