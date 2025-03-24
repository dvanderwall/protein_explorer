"""
Functions for kinase prediction based on structural and sequence similarity.
"""

import os
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Union
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Global variables to store loaded data
STRUCTURE_KINASE_SCORES = None
SEQUENCE_KINASE_SCORES = None

def load_kinase_scores(file_path: str = None, score_type: str = 'structure') -> pd.DataFrame:
    """
    Load kinase scores from the specified file.
    
    Args:
        file_path: Path to the kinase scores file (feather format)
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Pandas DataFrame with kinase scores
    """
    global STRUCTURE_KINASE_SCORES, SEQUENCE_KINASE_SCORES
    
    # Use the correct global variable based on score_type
    if score_type.lower() == 'structure':
        if STRUCTURE_KINASE_SCORES is not None:
            logger.info("Using cached structure kinase scores")
            return STRUCTURE_KINASE_SCORES
    else:
        if SEQUENCE_KINASE_SCORES is not None:
            logger.info("Using cached sequence kinase scores")
            return SEQUENCE_KINASE_SCORES
    
    # Find the file if path not provided
    if file_path is None:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(os.path.dirname(current_dir))
        
        if score_type.lower() == 'structure':
            file_path = os.path.join(parent_dir, 'Structure_Kinase_Scores.feather')
        else:
            file_path = os.path.join(parent_dir, 'Sequence_Kinase_Scores.feather')
    
    try:
        # Load the data
        logger.info(f"Loading {score_type} kinase scores from: {file_path}")
        scores_df = pd.read_feather(file_path)
        
        # Ensure the site ID column exists
        if 'site_id' not in scores_df.columns and scores_df.columns[0] != 'site_id':
            # Assume first column is the site ID
            scores_df = scores_df.rename(columns={scores_df.columns[0]: 'site_id'})
        
        # Create an index for faster lookups
        scores_df.set_index('site_id', inplace=True)
        
        # Store in the appropriate global variable
        if score_type.lower() == 'structure':
            STRUCTURE_KINASE_SCORES = scores_df
            logger.info(f"Loaded {len(scores_df)} structure kinase scores")
        else:
            SEQUENCE_KINASE_SCORES = scores_df
            logger.info(f"Loaded {len(scores_df)} sequence kinase scores")
        
        return scores_df
    except Exception as e:
        logger.error(f"Error loading {score_type} kinase scores: {e}")
        return pd.DataFrame()

def get_site_kinase_scores(site_id: str, score_type: str = 'structure') -> Dict:
    """
    Get kinase scores for a specific site.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with kinase names as keys and scores as values
    """
    # Load the data if not already loaded
    scores_df = load_kinase_scores(score_type=score_type)
    
    if scores_df.empty:
        logger.warning(f"No {score_type} kinase scores available")
        return {}
    
    try:
        # Get scores for the specific site
        if site_id in scores_df.index:
            site_scores = scores_df.loc[site_id]
            
            # Convert to dictionary, excluding the 'known_kinase' column if it exists
            scores_dict = {}
            for col in site_scores.index:
                if col != 'known_kinase':
                    # Make sure we convert all values to float
                    try:
                        value = site_scores[col]
                        if pd.notna(value):  # Check if not NaN
                            scores_dict[col] = float(value)
                        else:
                            scores_dict[col] = 0.0
                    except (ValueError, TypeError):
                        # Skip values that can't be converted to float
                        logger.warning(f"Skipping non-numeric value for kinase {col}: {site_scores[col]}")
            
            # Check if site has a known kinase
            known_kinase = site_scores.get('known_kinase', 'unlabeled')
            
            return {
                'known_kinase': known_kinase,
                'scores': scores_dict
            }
        else:
            logger.warning(f"Site {site_id} not found in {score_type} kinase scores")
            return {}
    except Exception as e:
        logger.error(f"Error getting {score_type} kinase scores for {site_id}: {e}")
        return {}

def predict_kinases(site_id: str, top_n: int = 5, score_type: str = 'structure') -> List[Dict]:
    """
    Get top N predicted kinases for a site.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        top_n: Number of top kinases to return
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        List of dictionaries with kinase names and scores
    """
    # Get all kinase scores for the site
    site_data = get_site_kinase_scores(site_id, score_type)
    
    if not site_data or 'scores' not in site_data:
        logger.warning(f"No {score_type} kinase scores available for {site_id}")
        return []
    
    # Get all kinase scores
    scores = site_data['scores']
    
    # Convert any string scores to float
    for kinase, score in list(scores.items()):
        try:
            if isinstance(score, str):
                # Try to convert string to float
                scores[kinase] = float(score)
        except (ValueError, TypeError):
            # If conversion fails, remove this item
            logger.warning(f"Removing non-numeric score for kinase {kinase}: {score}")
            del scores[kinase]
    
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

def compare_kinase_scores(site_ids: List[str], top_n: int = 5, score_type: str = 'structure') -> Dict:
    """
    Compare kinase scores across multiple sites.
    
    Args:
        site_ids: List of site IDs
        top_n: Number of top kinases to consider for each site
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with comparison data
    """
    if not site_ids:
        return {}
    
    # Get top kinases for each site
    site_kinases = {}
    all_kinases = set()
    
    for site_id in site_ids:
        top_kinases = predict_kinases(site_id, top_n, score_type)
        site_kinases[site_id] = top_kinases
        
        # Collect all unique kinases
        all_kinases.update([k['kinase'] for k in top_kinases])
    
    # Prepare comparison data
    comparison = {
        'sites': site_ids,
        'kinases': list(all_kinases),
        'data': {}
    }
    
    # Add scores for each site and kinase
    for site_id in site_ids:
        comparison['data'][site_id] = {}
        site_data = get_site_kinase_scores(site_id, score_type)
        
        if site_data and 'scores' in site_data:
            scores = site_data['scores']
            for kinase in all_kinases:
                comparison['data'][site_id][kinase] = scores.get(kinase, 0)
    
    return comparison

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
    if not site_ids:
        return {}
    
    # Load the data if not already loaded
    scores_df = load_kinase_scores(score_type=score_type)
    
    if scores_df.empty:
        logger.warning(f"No {score_type} kinase scores available")
        return {}
    
    try:
        # Get data for the specified sites
        sites_data = scores_df.loc[scores_df.index.isin(site_ids)]
        
        if sites_data.empty:
            logger.warning(f"None of the specified sites found in {score_type} kinase scores")
            return {}
        
        # Exclude 'known_kinase' column if it exists
        if 'known_kinase' in sites_data.columns:
            kinase_columns = [col for col in sites_data.columns if col != 'known_kinase']
        else:
            kinase_columns = list(sites_data.columns)
        
        # Convert all score columns to numeric, coercing errors to NaN
        for col in kinase_columns:
            sites_data[col] = pd.to_numeric(sites_data[col], errors='coerce')
        
        # Replace NaN with 0
        sites_data = sites_data.fillna(0)
        
        # Calculate mean score for each kinase across all sites
        mean_scores = sites_data[kinase_columns].mean()
        
        # Get top N kinases by mean score
        top_kinases = mean_scores.sort_values(ascending=False).head(top_n).index.tolist()
        
        # Prepare heatmap data
        heatmap_data = {
            'sites': list(sites_data.index),
            'kinases': top_kinases,
            'scores': []
        }
        
        # Add scores for each site and kinase
        for site in sites_data.index:
            for kinase in top_kinases:
                heatmap_data['scores'].append({
                    'site': site,
                    'kinase': kinase,
                    'score': float(sites_data.loc[site, kinase])
                })
        
        return heatmap_data
    except Exception as e:
        logger.error(f"Error generating heatmap data: {e}")
        return {}

def get_kinase_radar_data(site_id: str, top_n: int = 5, score_type: str = 'structure') -> Dict:
    """
    Get data for radar chart visualization of kinase scores.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        top_n: Number of top kinases to include
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with radar chart data
    """
    # Get top kinases for the site
    top_kinases = predict_kinases(site_id, top_n, score_type)
    
    if not top_kinases:
        return {}
    
    # Prepare radar chart data
    radar_data = {
        'labels': [k['kinase'] for k in top_kinases],
        'datasets': [{
            'label': f"{score_type.capitalize()} Kinase Scores",
            'data': [k['score'] for k in top_kinases]
        }]
    }
    
    return radar_data

def get_kinase_comparison_data(site_id: str, score_types: List[str] = ['structure', 'sequence'], top_n: int = 5) -> Dict:
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
        'MAPK': ['ERK1', 'ERK2', 'p38', 'JNK1', 'JNK2', 'JNK3'],
        'GSK': ['GSK3', 'GSK3A', 'GSK3B'],
        'CK': ['CK1', 'CK2', 'CSNK1', 'CSNK2'],
        'PKC': ['PKC', 'PKCALPHA', 'PKCBETA', 'PKCDELTA', 'PKCEPSILON', 'PKCGAMMA', 'PKCZETA'],
        'PKA': ['PKA', 'PKACA', 'PKACB', 'PKACG'],
        'AKT': ['AKT', 'AKT1', 'AKT2', 'AKT3'],
        'SRC': ['SRC', 'FYN', 'LCK', 'LYN', 'HCK', 'FGR', 'BLK', 'YES'],
        'CAMK': ['CAMK', 'CAMK1', 'CAMK2', 'CAMK4'],
        'ATM/ATR': ['ATM', 'ATR', 'DNAPK'],
        'PLK': ['PLK1', 'PLK2', 'PLK3', 'PLK4'],
        'AURORA': ['AURKA', 'AURKB', 'AURKC'],
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
            if any(member in kinase_name.upper() for member in members):
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