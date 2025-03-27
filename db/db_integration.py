"""
Integration module for KinoPlex application with Cloud SQL database.
Drop-in replacement for existing data access functions in KinoPlex.
"""

import os
import logging
import pandas as pd
import numpy as np
import re
from typing import Dict, List, Optional, Tuple, Union

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Flag to control database usage - can be set via environment variable
USE_DATABASE = os.environ.get('USE_DATABASE', 'True').lower() in ('true', '1', 't', 'yes')

# Original file-based functions - import if available
try:
    from protein_explorer.analysis.phospho_analyzer import (
        load_structural_similarity_data as file_load_structural_similarity_data,
        find_structural_matches as file_find_structural_matches,
        get_phosphosite_data as file_get_phosphosite_data,
        enhance_phosphosite as file_enhance_phosphosite,
        enhance_structural_matches as file_enhance_structural_matches,
        analyze_protein as file_analyze_protein,
        get_phosphosites as file_get_phosphosites
    )
    
    from protein_explorer.analysis.sequence_analyzer import (
        find_sequence_matches as file_find_sequence_matches,
        analyze_motif_conservation as file_analyze_motif_conservation,
        create_sequence_network_data as file_create_sequence_network_data,
        create_sequence_motif_visualization as file_create_sequence_motif_visualization
    )
    
    from protein_explorer.analysis.kinase_predictor import (
        load_kinase_scores as file_load_kinase_scores,
        get_site_kinase_scores as file_get_site_kinase_scores,
        predict_kinases as file_predict_kinases,
        get_heatmap_data as file_get_heatmap_data,
        get_kinase_comparison_data as file_get_kinase_comparison_data,
        get_known_kinase_info as file_get_known_kinase_info,
        categorize_kinases_by_family as file_categorize_kinases_by_family
    )
    
    from protein_explorer.analysis.network_kinase_predictor import (
        get_similar_sites as file_get_similar_sites,
        get_network_kinase_scores as file_get_network_kinase_scores,
        compute_aggregated_kinase_scores as file_compute_aggregated_kinase_scores,
        predict_kinases_network as file_predict_kinases_network,
        get_network_heatmap_data as file_get_network_heatmap_data,
        get_network_kinase_comparison as file_get_network_kinase_comparison,
        get_kinase_family_distribution_network as file_get_kinase_family_distribution_network
    )
    
    FILE_FUNCTIONS_AVAILABLE = True
except ImportError as e:
    import traceback
    print(f"Error importing file-based functions: {e}")
    print(traceback.format_exc())
    logger.warning("Original file-based functions not found, database functions will be used exclusively")
    FILE_FUNCTIONS_AVAILABLE = False

# Import database functions - this will fail if the module is not installed yet
try:
    # Use relative import to avoid circular dependency issues
    from .cloud_sql_connector import (
        get_db_engine,
        execute_query,
        cache_result,
        clear_cache,
        get_phosphosite_data as db_get_phosphosite_data,
        get_kinase_scores as db_get_kinase_scores,
        get_all_phosphosites as db_get_all_phosphosites,
        get_similar_sites as db_get_similar_sites,
        get_network_kinase_scores as db_get_network_kinase_scores,
        get_heatmap_data as db_get_heatmap_data,
        health_check
    )
    
    # These functions might not exist in cloud_sql_connector.py,
    # so import them conditionally
    try:
        from .cloud_sql_connector import get_structural_matches as db_get_structural_matches
    except ImportError:
        db_get_structural_matches = None
        
    try:
        from .cloud_sql_connector import get_sequence_matches as db_get_sequence_matches
    except ImportError:
        db_get_sequence_matches = None
        
    DB_FUNCTIONS_AVAILABLE = True
    print("Successfully imported database functions")
except ImportError as e:
    import traceback
    print(f"Failed to import cloud_sql_connector: {e}")
    print(traceback.format_exc())
    logger.warning(f"Database connector module not found: {e}")
    logger.warning("Falling back to file-based functions")
    DB_FUNCTIONS_AVAILABLE = False
    # Force USE_DATABASE to False if database functions are not available
    USE_DATABASE = False

# Import other required modules
try:
    from protein_explorer.data.scaffold import get_protein_by_id, get_alphafold_structure
    from protein_explorer.analysis.phospho import analyze_phosphosites
except ImportError as e:
    import traceback
    print(f"Error importing essential modules: {e}")
    print(traceback.format_exc())
    logger.warning("Could not import essential KinoPlex modules")


# Integration functions

def find_structural_matches(
    uniprot_id: str, 
    phosphosites: List[Dict], 
    parquet_file: str = None, 
    top_n: Optional[int] = None
) -> Dict[str, List[Dict]]:
    """
    Find structural matches for phosphosites in the database or feather file.
    Compatible with the original function signature for drop-in replacement.
    
    Args:
        uniprot_id: UniProt ID of the protein
        phosphosites: List of phosphosite dictionaries from analyze_phosphosites
        parquet_file: Path to the data file (ignored when using database)
        top_n: Number of top matches to return per site
        
    Returns:
        Dictionary mapping site IDs to lists of match dictionaries
    """
    if USE_DATABASE and DB_FUNCTIONS_AVAILABLE:
        try:
            # Database version
            logger.info(f"Using database to find structural matches for {uniprot_id}")
            
            structural_matches = {}
            
            for site_data in phosphosites:
                if 'site' not in site_data:
                    continue
                    
                # Get the site string (e.g., "S123")
                site = site_data['site']
                
                # Query database for matches
                matches_df = db_get_structural_matches(uniprot_id, site, top_n)
                
                if matches_df.empty:
                    structural_matches[site] = []
                    continue
                
                # Convert dataframe rows to dictionaries
                matches = []
                for _, row in matches_df.iterrows():
                    match_dict = {
                        'query_uniprot': row.get('query_uniprot', uniprot_id),
                        'query_site': row.get('query_site', site),
                        'target_uniprot': row.get('target_uniprot'),
                        'target_site': row.get('target_site'),
                        'rmsd': row.get('RMSD', row.get('rmsd')),
                    }
                    
                    # Include any other columns that might be present
                    for col in row.index:
                        if col not in match_dict and pd.notna(row[col]):
                            match_dict[col] = row[col]
                    
                    matches.append(match_dict)
                
                # Add to results
                structural_matches[site] = matches
            
            return structural_matches
            
        except Exception as e:
            logger.error(f"Error finding structural matches using database: {e}")
            
            if FILE_FUNCTIONS_AVAILABLE:
                logger.info("Falling back to file-based function")
                return file_find_structural_matches(uniprot_id, phosphosites, parquet_file, top_n)
            else:
                # Return empty dictionary if no fallback available
                return {}
    else:
        # Use original file-based function
        if FILE_FUNCTIONS_AVAILABLE:
            return file_find_structural_matches(uniprot_id, phosphosites, parquet_file, top_n)
        else:
            logger.error("File-based functions not available")
            return {}

def enhance_structural_matches(matches: List[Dict], site: str) -> List[Dict]:
    """
    Enhance structural matches with supplementary data.
    Compatible with the original function signature for drop-in replacement.
    
    Args:
        matches: List of structural match dictionaries
        site: Query site string for logging
        
    Returns:
        Enhanced list of matches
    """
    if not matches:
        return matches
    
    enhanced_matches = []
    for match in matches:
        # Skip self-matches (RMSD ≈ 0)
        if match.get('rmsd', 0) < 0.01:
            continue
        
        # Get target info
        target_uniprot = match.get('target_uniprot')
        target_site = match.get('target_site')
        
        if not target_uniprot or not target_site:
            continue
        
        # Parse site number from target_site
        import re
        site_match = re.match(r'([STY])?(\d+)', str(target_site))
        if site_match:
            site_type = site_match.group(1) if site_match.group(1) else ''
            resno = int(site_match.group(2))
            target_id = f"{target_uniprot}_{resno}"
            
            # Get supplementary data
            target_supp = get_phosphosite_data(target_id)
            if target_supp:
                # Create enhanced match
                enhanced_match = match.copy()
                
                # Add supplementary data
                if 'motif_plddt' in target_supp and target_supp['motif_plddt'] is not None:
                    enhanced_match['plddt'] = float(target_supp['motif_plddt'])
                
                if 'nearby_count' in target_supp and target_supp['nearby_count'] is not None:
                    enhanced_match['nearby_count'] = target_supp['nearby_count']
                
                if 'SITE_+/-7_AA' in target_supp and target_supp['SITE_+/-7_AA'] is not None:
                    enhanced_match['motif'] = target_supp['SITE_+/-7_AA']
                
                # Add additional fields
                for key in ['site_plddt', 'surface_accessibility']:
                    if key in target_supp and target_supp[key] is not None:
                        enhanced_match[key] = target_supp[key]
                
                enhanced_matches.append(enhanced_match)
                continue
        
        # If no supplementary data, just add the original match
        enhanced_matches.append(match)
    
    return enhanced_matches

def get_phosphosite_data(site_id: str) -> Optional[Dict]:
    """
    Get phosphosite data from the database or feather file.
    Compatible with the original function signature for drop-in replacement.
    
    Args:
        site_id: The site ID in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary with phosphosite data or None if not found
    """
    if USE_DATABASE and DB_FUNCTIONS_AVAILABLE:
        try:
            # Database version
            logger.debug(f"Using database to get phosphosite data for {site_id}")
            return db_get_phosphosite_data(site_id)
        except Exception as e:
            logger.error(f"Error getting phosphosite data using database: {e}")
            
            if FILE_FUNCTIONS_AVAILABLE:
                logger.info("Falling back to file-based function")
                return file_get_phosphosite_data(site_id)
            else:
                return None
    else:
        # Use original file-based function
        if FILE_FUNCTIONS_AVAILABLE:
            return file_get_phosphosite_data(site_id)
        else:
            logger.error("File-based functions not available")
            return None

def enhance_phosphosite(phosphosite: Dict, uniprot_id: str) -> Dict:
    """
    Enhance a phosphosite dictionary with supplementary data.
    Compatible with the original function signature for drop-in replacement.
    
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
    
    if 'SITE_+/-7_AA' in supp_data and supp_data['SITE_+/-7_AA'] is not None:
        enhanced_site['motif'] = supp_data['SITE_+/-7_AA']
    
    # Add any other supplementary fields
    for key in ['site_plddt', 'surface_accessibility']:
        if key in supp_data and supp_data[key] is not None:
            enhanced_site[key] = supp_data[key]
    
    return enhanced_site

def find_sequence_matches(
    site_id: str, 
    top_n: int = 200, 
    min_similarity: float = 0.4
) -> List[Dict]:
    """
    Find sequence similarity matches from the database or feather file.
    Compatible with the original function signature for drop-in replacement.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        top_n: Maximum number of results to return
        min_similarity: Minimum similarity score to include (0-1)
        
    Returns:
        List of dictionaries with match information
    """
    if USE_DATABASE and DB_FUNCTIONS_AVAILABLE:
        try:
            # Database version
            logger.debug(f"Using database to find sequence matches for {site_id}")
            
            # Get matches from database
            matches_df = db_get_sequence_matches(site_id, min_similarity, top_n)
            
            if matches_df.empty:
                return []
            
            # Convert dataframe rows to the expected dictionary format
            matches = []
            for _, row in matches_df.iterrows():
                match_dict = {
                    'query_id': row.get('query_id', site_id),
                    'target_id': row.get('target_id'),
                    'target_uniprot': row.get('target_uniprot'),
                    'target_site': row.get('target_site'),
                    'similarity': float(row.get('Similarity', row.get('similarity', 0.0))),
                }
                
                # Add motif if available
                if 'target_motif' in row and pd.notna(row['target_motif']):
                    match_dict['motif'] = row['target_motif']
                else:
                    # Try to get motif from supplementary data
                    target_id = row.get('target_id')
                    if target_id:
                        supp_data = get_phosphosite_data(target_id)
                        if supp_data and 'SITE_+/-7_AA' in supp_data and supp_data['SITE_+/-7_AA']:
                            match_dict['motif'] = supp_data['SITE_+/-7_AA']
                
                # Include any other columns that might be present
                for col in row.index:
                    if col not in match_dict and pd.notna(row[col]):
                        match_dict[col] = row[col]
                
                matches.append(match_dict)
            
            return matches
        
        except Exception as e:
            logger.error(f"Error finding sequence matches using database: {e}")
            
            if FILE_FUNCTIONS_AVAILABLE:
                logger.info("Falling back to file-based function")
                return file_find_sequence_matches(site_id, top_n, min_similarity)
            else:
                return []
    else:
        # Use original file-based function
        if FILE_FUNCTIONS_AVAILABLE:
            return file_find_sequence_matches(site_id, top_n, min_similarity)
        else:
            logger.error("File-based functions not available")
            return []

def get_site_kinase_scores(site_id: str, score_type: str = 'structure') -> Dict:
    """
    Get kinase scores for a specific site from the database or feather file.
    Compatible with the original function signature for drop-in replacement.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with kinase scores
    """
    if USE_DATABASE and DB_FUNCTIONS_AVAILABLE:
        try:
            # Database version
            logger.debug(f"Using database to get {score_type} kinase scores for {site_id}")
            return db_get_kinase_scores(site_id, score_type)
        except Exception as e:
            logger.error(f"Error getting {score_type} kinase scores using database: {e}")
            
            if FILE_FUNCTIONS_AVAILABLE:
                logger.info("Falling back to file-based function")
                return file_get_site_kinase_scores(site_id, score_type)
            else:
                return {}
    else:
        # Use original file-based function
        if FILE_FUNCTIONS_AVAILABLE:
            return file_get_site_kinase_scores(site_id, score_type)
        else:
            logger.error("File-based functions not available")
            return {}

def predict_kinases(site_id: str, top_n: int = 5, score_type: str = 'structure') -> List[Dict]:
    """
    Get top N predicted kinases for a site.
    Compatible with the original function signature for drop-in replacement.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        top_n: Number of top kinases to return
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        List of dictionaries with kinase names and scores
    """
    # Get all kinase scores
    site_data = get_site_kinase_scores(site_id, score_type)
    
    if not site_data or 'scores' not in site_data:
        logger.warning(f"No {score_type} kinase scores available for {site_id}")
        return []
    
    # Get all kinase scores
    scores = site_data['scores']
    
    # Sort by score (descending)
    sorted_scores = sorted(scores.items(), key=lambda x: float(x[1]), reverse=True)
    
    # Return top N kinases
    top_kinases = []
    for kinase, score in sorted_scores[:top_n]:
        top_kinases.append({
            'kinase': kinase,
            'score': float(score)  # Ensure score is a float
        })
    
    return top_kinases

def get_similar_sites(
    site_id: str, 
    uniprot_id: Optional[str] = None, 
    site_data: Optional[Dict] = None,
    similarity_threshold: float = 0.6, 
    rmsd_threshold: float = 3.0
) -> List[str]:
    """
    Get a list of sites similar to the query site based on both
    sequence similarity and structural similarity.
    
    Args:
        site_id: The query site ID (format: UniProtID_ResidueNumber)
        uniprot_id: The UniProt ID of the protein (extracted from site_id if not provided)
        site_data: Site data dictionary (optional, will be retrieved if not provided)
        similarity_threshold: Minimum sequence similarity score to include
        rmsd_threshold: Maximum RMSD value to include for structural matches
        
    Returns:
        List of site IDs (including the query site) for aggregation
    """
    if USE_DATABASE and DB_FUNCTIONS_AVAILABLE:
        try:
            # Database version
            logger.debug(f"Using database to find similar sites for {site_id}")
            return db_get_similar_sites(site_id, similarity_threshold, rmsd_threshold)
        except Exception as e:
            logger.error(f"Error finding similar sites using database: {e}")
            
            if FILE_FUNCTIONS_AVAILABLE:
                logger.info("Falling back to file-based function")
                return file_get_similar_sites(site_id, uniprot_id, site_data, 
                                            similarity_threshold, rmsd_threshold)
            else:
                return [site_id]  # Return just the query site in case of error
    else:
        # Use original file-based function
        if FILE_FUNCTIONS_AVAILABLE:
            return file_get_similar_sites(site_id, uniprot_id, site_data, 
                                         similarity_threshold, rmsd_threshold)
        else:
            logger.error("File-based functions not available")
            return [site_id]

def get_network_kinase_scores(
    site_id: str, 
    score_type: str = 'structure',
    similarity_threshold: float = 0.6, 
    rmsd_threshold: float = 3.0
) -> Dict[str, Dict[str, float]]:
    """
    Get kinase scores for a network of similar sites.
    
    Args:
        site_id: The query site ID
        score_type: Type of scores - 'structure' or 'sequence'
        similarity_threshold: Minimum sequence similarity score
        rmsd_threshold: Maximum RMSD value for structural matches
        
    Returns:
        Dictionary mapping site IDs to score dictionaries
    """
    if USE_DATABASE and DB_FUNCTIONS_AVAILABLE:
        try:
            # Database version
            logger.debug(f"Using database to get network kinase scores for {site_id}")
            return db_get_network_kinase_scores(site_id, score_type, 
                                              similarity_threshold, rmsd_threshold)
        except Exception as e:
            logger.error(f"Error getting network kinase scores using database: {e}")
            
            if FILE_FUNCTIONS_AVAILABLE:
                logger.info("Falling back to file-based function")
                return file_get_network_kinase_scores(site_id, score_type,
                                                    similarity_threshold, rmsd_threshold)
            else:
                return {}
    else:
        # Use original file-based function
        if FILE_FUNCTIONS_AVAILABLE:
            return file_get_network_kinase_scores(site_id, score_type,
                                                similarity_threshold, rmsd_threshold)
        else:
            logger.error("File-based functions not available")
            return {}

def compute_aggregated_kinase_scores(
    site_id: str, 
    score_type: str = 'structure',
    similarity_threshold: float = 0.6,
    rmsd_threshold: float = 3.0
) -> List[Dict]:
    """
    Compute aggregated kinase scores across a network of similar sites.
    
    Args:
        site_id: The query site ID
        score_type: Type of scores - 'structure' or 'sequence'
        similarity_threshold: Minimum sequence similarity score
        rmsd_threshold: Maximum RMSD value for structural matches
        
    Returns:
        List of dictionaries with aggregated kinase scores and statistics
    """
    # Get scores for all similar sites
    network_scores = get_network_kinase_scores(site_id, score_type, 
                                             similarity_threshold, rmsd_threshold)
    
    if not network_scores:
        logger.warning(f"No network scores available for {site_id}")
        return []
    
    # Get all unique kinases
    all_kinases = set()
    for site_scores in network_scores.values():
        all_kinases.update(site_scores.keys())
    
    # Calculate statistics for each kinase
    aggregated_scores = []
    for kinase in all_kinases:
        # Collect scores for this kinase across all sites
        kinase_scores = [site_scores.get(kinase, 0) for site_scores in network_scores.values()]
        
        # Filter out zero scores (missing values)
        non_zero_scores = [score for score in kinase_scores if score > 0]
        
        # Only include kinases that have scores for at least 2 sites
        if len(non_zero_scores) < 2:
            continue
        
        # Calculate statistics
        mean_score = sum(non_zero_scores) / len(non_zero_scores) if non_zero_scores else 0
        
        # Sort scores for median calculation
        sorted_scores = sorted(non_zero_scores)
        if sorted_scores:
            median_score = sorted_scores[len(sorted_scores)//2] if len(sorted_scores) % 2 == 1 else (sorted_scores[len(sorted_scores)//2-1] + sorted_scores[len(sorted_scores)//2])/2
            max_score = max(sorted_scores)
            min_score = min(sorted_scores)
            variability = max_score - min_score
        else:
            median_score = max_score = min_score = variability = 0
        
        # Skip if all zeros
        if mean_score == 0 and median_score == 0:
            continue
            
        # Create aggregated score entry
        aggregated_scores.append({
            'kinase': kinase,
            'mean_score': mean_score,
            'median_score': median_score,
            'max_score': max_score,
            'min_score': min_score,
            'variability': variability,
            'sample_size': len(non_zero_scores),
            'scores': non_zero_scores  # Store all scores for detailed analysis
        })
    
    # Sort by mean score (descending)
    return sorted(aggregated_scores, key=lambda x: x['mean_score'], reverse=True)

def predict_kinases_network(
    site_id: str, 
    top_n: int = 5, 
    score_type: str = 'structure',
    similarity_threshold: float = 0.6, 
    rmsd_threshold: float = 3.0
) -> List[Dict]:
    """
    Get top N predicted kinases based on aggregated network scores.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        top_n: Number of top kinases to return
        score_type: Type of scores - 'structure' or 'sequence'
        similarity_threshold: Minimum sequence similarity score
        rmsd_threshold: Maximum RMSD value for structural matches
        
    Returns:
        List of dictionaries with kinase aggregated scores
    """
    # Get aggregated scores
    aggregated_scores = compute_aggregated_kinase_scores(
        site_id, score_type, similarity_threshold, rmsd_threshold
    )
    
    # Return top N kinases
    return aggregated_scores[:top_n]

def get_heatmap_data(site_ids: List[str], top_n: int = 10, score_type: str = 'structure') -> Dict:
    """
    Get data for heatmap visualization of kinase scores.
    
    Args:
        site_ids: List of site IDs
        top_n: Number of top kinases to include
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with heatmap data
    """
    if USE_DATABASE and DB_FUNCTIONS_AVAILABLE:
        try:
            # Database version
            logger.debug(f"Using database to get heatmap data for {len(site_ids)} sites")
            return db_get_heatmap_data(site_ids, top_n, score_type)
        except Exception as e:
            logger.error(f"Error getting heatmap data using database: {e}")
            
            if FILE_FUNCTIONS_AVAILABLE:
                logger.info("Falling back to file-based function")
                return file_get_heatmap_data(site_ids, top_n, score_type)
            else:
                return {'sites': [], 'kinases': [], 'scores': []}
    else:
        # Use original file-based function
        if FILE_FUNCTIONS_AVAILABLE:
            return file_get_heatmap_data(site_ids, top_n, score_type)
        else:
            logger.error("File-based functions not available")
            return {'sites': [], 'kinases': [], 'scores': []}

def get_network_heatmap_data(
    site_id: str, 
    top_n: int = 10, 
    score_type: str = 'structure',
    similarity_threshold: float = 0.6, 
    rmsd_threshold: float = 3.0
) -> Dict:
    """
    Get data for heatmap visualization of kinase scores across a network of similar sites.
    
    Args:
        site_id: The query site ID
        top_n: Number of top kinases to include
        score_type: Type of scores - 'structure' or 'sequence'
        similarity_threshold: Minimum sequence similarity score
        rmsd_threshold: Maximum RMSD value for structural matches
        
    Returns:
        Dictionary with heatmap data
    """
    # Get network scores
    network_scores = get_network_kinase_scores(site_id, score_type, 
                                             similarity_threshold, rmsd_threshold)
    
    if not network_scores:
        logger.warning(f"No network scores available for {site_id}")
        return {'sites': [], 'kinases': [], 'scores': []}
    
    # Get aggregated scores to identify top kinases
    aggregated_scores = compute_aggregated_kinase_scores(
        site_id, score_type, similarity_threshold, rmsd_threshold
    )
    
    # Get top N kinases by mean score
    top_kinases = [score['kinase'] for score in aggregated_scores[:top_n]]
    
    # Prepare heatmap data
    heatmap_data = {
        'sites': list(network_scores.keys()),
        'kinases': top_kinases,
        'scores': []
    }
    
    # Add scores for each site and kinase
    for site, site_scores in network_scores.items():
        for kinase in top_kinases:
            score = site_scores.get(kinase, 0)
            
            # Always add data point even if score is 0 to ensure complete heatmap
            heatmap_data['scores'].append({
                'site': site,
                'kinase': kinase,
                'score': float(score)
            })
    
    return heatmap_data

def get_kinase_comparison_data(
    site_id: str, 
    score_types: List[str] = ['structure', 'sequence'], 
    top_n: int = 5
) -> Dict:
    """
    Get comparison data between structure and sequence kinase scores.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_types: List of score types to compare
        top_n: Number of top kinases to include
        
    Returns:
        Dictionary with comparison data
    """
    # Get all unique kinases from both score types
    all_kinases = set()
    
    for score_type in score_types:
        top_kinases = predict_kinases(site_id, top_n, score_type)
        all_kinases.update([k['kinase'] for k in top_kinases])
    
    # Prepare comparison data
    comparison_data = {
        'kinases': list(all_kinases),
        'datasets': []
    }
    
    # Add scores for each score type
    for score_type in score_types:
        site_data = get_site_kinase_scores(site_id, score_type)
        
        if site_data and 'scores' in site_data:
            scores = site_data['scores']
            
            dataset = {
                'label': f"{score_type.capitalize()} Score",
                'data': [scores.get(kinase, 0) for kinase in all_kinases]
            }
            
            comparison_data['datasets'].append(dataset)
    
    return comparison_data

def get_network_kinase_comparison(
    site_id: str, 
    score_types: List[str] = ['structure', 'sequence'],
    top_n: int = 5,
    similarity_threshold: float = 0.6, 
    rmsd_threshold: float = 3.0
) -> Dict:
    """
    Get comparison data between structure and sequence network kinase predictions.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_types: List of score types to compare
        top_n: Number of top kinases to include 
        similarity_threshold: Minimum sequence similarity score
        rmsd_threshold: Maximum RMSD value for structural matches
        
    Returns:
        Dictionary with comparison data
    """
    # Get predictions for each score type
    all_predictions = {}
    for score_type in score_types:
        predictions = predict_kinases_network(
            site_id, top_n=top_n, score_type=score_type,
            similarity_threshold=similarity_threshold,
            rmsd_threshold=rmsd_threshold
        )
        all_predictions[score_type] = predictions
    
    # Get all unique kinases
    all_kinases = set()
    for score_type, predictions in all_predictions.items():
        all_kinases.update([p['kinase'] for p in predictions])
    
    # Prepare comparison data
    comparison_data = {
        'kinases': list(all_kinases),
        'datasets': []
    }
    
    # Add data for each score type
    for score_type, predictions in all_predictions.items():
        # Create a lookup dictionary
        pred_dict = {p['kinase']: p for p in predictions}
        
        # Add dataset
        dataset = {
            'label': f"{score_type.capitalize()} Network Score",
            'data': [pred_dict.get(kinase, {}).get('mean_score', 0) for kinase in all_kinases],
            'error_bars': [pred_dict.get(kinase, {}).get('variability', 0) for kinase in all_kinases]
        }
        
        comparison_data['datasets'].append(dataset)
    
    return comparison_data

def get_all_phosphosites(uniprot_id: str) -> List[Dict]:
    """
    Get all phosphosites for a protein from the database.
    
    Args:
        uniprot_id: UniProt ID of the protein
        
    Returns:
        List of dictionaries with phosphosite information
    """
    if USE_DATABASE and DB_FUNCTIONS_AVAILABLE:
        try:
            # Database version
            logger.info(f"Using database to get all phosphosites for {uniprot_id}")
            return db_get_all_phosphosites(uniprot_id)
        except Exception as e:
            logger.error(f"Error getting phosphosites using database: {e}")
            
            # Try to fall back to file-based analysis
            if FILE_FUNCTIONS_AVAILABLE:
                try:
                    logger.info("Falling back to file-based function")
                    if hasattr(file_get_phosphosites, "__call__"):
                        return file_get_phosphosites(uniprot_id)
                    else:
                        # Fall back to other methods
                        return _analyze_phosphosites_fallback(uniprot_id)
                except Exception as fallback_error:
                    logger.error(f"Fallback function also failed: {fallback_error}")
                    return []
            else:
                return []
    else:
        # Use original file-based function sequence
        if FILE_FUNCTIONS_AVAILABLE:
            try:
                if hasattr(file_get_phosphosites, "__call__"):
                    return file_get_phosphosites(uniprot_id)
                else:
                    # Fall back to other methods
                    return _analyze_phosphosites_fallback(uniprot_id)
            except Exception as e:
                logger.error(f"Error in file-based phosphosite analysis: {e}")
                return []
        else:
            logger.error("File-based functions not available")
            return []

def _analyze_phosphosites_fallback(uniprot_id: str) -> List[Dict]:
    """
    Fallback method to analyze phosphosites using direct structure and sequence analysis.
    
    Args:
        uniprot_id: UniProt ID of the protein
        
    Returns:
        List of dictionaries with phosphosite information
    """
    try:
        # Get protein data
        protein_data = get_protein_by_id(uniprot_id=uniprot_id)
        
        # Get sequence
        sequence = protein_data.get('metadata', {}).get('sequence', {}).get('value')
        if not sequence:
            logger.warning(f"Protein sequence not found for {uniprot_id}")
            return []
            
        # Get structure
        structure = get_alphafold_structure(uniprot_id)
        if not structure:
            logger.warning(f"Protein structure not found for {uniprot_id}")
            return []
        
        # Analyze phosphosites
        phosphosites = analyze_phosphosites(sequence, structure, uniprot_id)
        
        # Enhance with supplementary data if available
        for site in phosphosites:
            if 'resno' in site:
                site_id = f"{uniprot_id}_{site['resno']}"
                supp_data = get_phosphosite_data(site_id)
                
                if supp_data:
                    for key, value in supp_data.items():
                        if key not in site and value is not None:
                            site[key] = value
        
        return phosphosites
    except Exception as e:
        logger.error(f"Error in phosphosite fallback analysis: {e}")
        return []

def get_known_kinase_info(site_id: str, score_type: str = 'structure') -> Dict:
    """
    Get information about the known kinase for a site, if available.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_type: Type of scores to check
        
    Returns:
        Dictionary with known kinase information
    """
    # Get site data
    site_data = get_site_kinase_scores(site_id, score_type)
    
    if not site_data:
        return {'has_known_kinase': False}
    
    known_kinase = site_data.get('known_kinase', 'unlabeled')
    
    if known_kinase == 'unlabeled':
        return {'has_known_kinase': False}
    
    return {
        'has_known_kinase': True,
        'kinase': known_kinase
    }

def categorize_kinases_by_family(kinases: List[Dict]) -> Dict:
    """
    Categorize kinases by family.
    
    Args:
        kinases: List of dictionaries with kinase names and scores
        
    Returns:
        Dictionary with kinase families and their scores
    """
    # Define kinase families
    kinase_families = {
        'CDK': ['CDK1', 'CDK2', 'CDK4', 'CDK5', 'CDK6', 'CDK7', 'CDK8', 'CDK9'],
        'MAPK': ['ERK1', 'ERK2', 'p38', 'JNK1', 'JNK2', 'JNK3', 'P38A', 'P38B', 'P38D', 'P38G'],
        'GSK': ['GSK3', 'GSK3A', 'GSK3B'],
        'CK': ['CK1', 'CK2', 'CSNK1', 'CSNK2', 'CK1A', 'CK1D', 'CK1E', 'CK1G1', 'CK1G2', 'CK2A1', 'CK2A2', 'CK2B'],
        'PKC': ['PKC', 'PKCALPHA', 'PKCBETA', 'PKCDELTA', 'PKCEPSILON', 'PKCGAMMA', 'PKCZETA', 'PKCA', 'PKCB', 'PKCD', 'PKCE', 'PKCG', 'PKCZ'],
        'PKA': ['PKA', 'PKACA', 'PKACB', 'PKACG'],
        'AKT': ['AKT', 'AKT1', 'AKT2', 'AKT3'],
        'SRC': ['SRC', 'FYN', 'LCK', 'LYN', 'HCK', 'FGR', 'BLK', 'YES'],
        'CAMK': ['CAMK', 'CAMK1', 'CAMK2', 'CAMK4', 'CAMK1A', 'CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G'],
        'ATM/ATR': ['ATM', 'ATR', 'DNAPK'],
        'PLK': ['PLK1', 'PLK2', 'PLK3', 'PLK4'],
        'AURORA': ['AURKA', 'AURKB', 'AURKC', 'AurA', 'AurB', 'AurC'],
        'Other': []
    }
    
    # Categorize kinases
    family_scores = {}
    
    for kinase_data in kinases:
        kinase_name = kinase_data['kinase']
        kinase_score = kinase_data['score']
        
        # Find the family for this kinase
        assigned = False
        for family, members in kinase_families.items():
            if any(member.upper() in kinase_name.upper() for member in members):
                if family not in family_scores:
                    family_scores[family] = 0
                family_scores[family] += kinase_score
                assigned = True
                break
        
        # If not assigned to any family, put in 'Other'
        if not assigned:
            if 'Other' not in family_scores:
                family_scores['Other'] = 0
            family_scores['Other'] += kinase_score
    
    # Return sorted by score (descending)
    return dict(sorted(family_scores.items(), key=lambda x: x[1], reverse=True))

def get_kinase_family_distribution_network(
    site_id: str, 
    score_type: str = 'structure',
    similarity_threshold: float = 0.6, 
    rmsd_threshold: float = 3.0
) -> Dict:
    """
    Get kinase family distribution based on network aggregated scores.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_type: Type of scores - 'structure' or 'sequence'
        similarity_threshold: Minimum sequence similarity score
        rmsd_threshold: Maximum RMSD value for structural matches
        
    Returns:
        Dictionary with kinase families and their scores
    """
    # Get top kinases
    predictions = predict_kinases_network(
        site_id, top_n=20, score_type=score_type,
        similarity_threshold=similarity_threshold,
        rmsd_threshold=rmsd_threshold
    )
    
    # Convert to format expected by categorize_kinases_by_family
    kinases_for_categorization = []
    for pred in predictions:
        kinases_for_categorization.append({
            'kinase': pred['kinase'],
            'score': pred['mean_score']
        })
    
    # Use categorize_kinases_by_family
    return categorize_kinases_by_family(kinases_for_categorization)

def analyze_motif_conservation(
    matches: List[Dict], 
    query_motif: str = None
) -> Dict:
    """
    Analyze conservation patterns in motifs of sequence-similar sites.
    
    Args:
        matches: List of match dictionaries from find_sequence_matches
        query_motif: Motif of the query site
        
    Returns:
        Dictionary with conservation analysis results
    """
    if FILE_FUNCTIONS_AVAILABLE and hasattr(file_analyze_motif_conservation, "__call__"):
        return file_analyze_motif_conservation(matches, query_motif)
    
    # Fallback implementation
    logger.warning("Using fallback implementation for motif conservation analysis")
    
    # If no matches or no motifs available, return empty results
    if not matches:
        return {
            'motif_count': 0,
            'conserved_positions': [],
            'consensus_motif': ""
        }
    
    # Collect motifs from matches
    motifs = []
    for match in matches:
        if 'motif' in match and match['motif']:
            motifs.append(match['motif'])
    
    # Add query motif if provided
    if query_motif:
        motifs = [query_motif] + motifs
    
    if not motifs:
        return {
            'motif_count': 0,
            'conserved_positions': [],
            'consensus_motif': ""
        }
    
    # Standardize motif length (assume phosphosite is in the middle)
    std_motifs = []
    for motif in motifs:
        center_pos = len(motif) // 2
        site_char = motif[center_pos]
        
        before_site = motif[:center_pos]
        after_site = motif[center_pos+1:]
        
        # Ensure we have 7 positions before and after (or pad with X)
        if len(before_site) < 7:
            before_site = 'X' * (7 - len(before_site)) + before_site
        else:
            before_site = before_site[-7:]
            
        if len(after_site) < 7:
            after_site = after_site + 'X' * (7 - len(after_site))
        else:
            after_site = after_site[:7]
            
        std_motifs.append(before_site + site_char + after_site)
    
    # Calculate conservation at each position
    position_counts = []
    for i in range(15):  # 15 positions in standardized motif (-7 to +7)
        counts = {}
        for motif in std_motifs:
            if i < len(motif):
                aa = motif[i]
                if aa not in counts:
                    counts[aa] = 0
                counts[aa] += 1
        position_counts.append(counts)
    
    # Generate consensus motif
    consensus = ""
    for i in range(15):
        counts = position_counts[i]
        if not counts:
            consensus += "X"
        else:
            max_aa = max(counts.items(), key=lambda x: x[1])[0]
            consensus += max_aa
    
    # Find conserved positions
    conserved_positions = []
    for i in range(15):
        if i == 7:  # Skip the phosphosite itself
            continue
            
        position = i - 7  # Convert to -7 to +7 positions
        counts = position_counts[i]
        if not counts:
            continue
            
        total = sum(counts.values())
        max_count = max(counts.values())
        max_aa = max(counts.items(), key=lambda x: x[1])[0]
        
        if max_count / total >= 0.5:  # 50% or more conservation
            conserved_positions.append({
                'position': position,
                'amino_acid': max_aa,
                'frequency': max_count / total * 100
            })
    
    return {
        'motif_count': len(motifs),
        'consensus_motif': consensus,
        'conserved_positions': conserved_positions
    }

# Function to toggle database usage at runtime
def use_database(enabled: bool = True):
    """Toggle database usage at runtime."""
    global USE_DATABASE
    
    # Only change if database functions are available
    if enabled and not DB_FUNCTIONS_AVAILABLE:
        logger.warning("Database functions not available, cannot enable database usage")
        return False
    
    USE_DATABASE = enabled
    logger.info(f"Database usage {'enabled' if enabled else 'disabled'}")
    return True

# Database health check
def check_database_health() -> Dict:
    """Check database connection and return health status."""
    if not DB_FUNCTIONS_AVAILABLE:
        return {
            'status': 'unavailable',
            'message': 'Database functions not available',
            'timestamp': pd.Timestamp.now().isoformat()
        }
    
    try:
        health_status = health_check()
        return health_status
    except Exception as e:
        import traceback
        logger.error(f"Database health check failed with exception: {e}")
        logger.error(traceback.format_exc())
        return {
            'status': 'error',
            'error': str(e),
            'timestamp': pd.Timestamp.now().isoformat()
        }