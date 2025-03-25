"""
Enhanced Kinase Prediction Based on Phosphosite Networks.

This module extends the kinase_predictor.py functionality to incorporate
network-based aggregation of kinase predictions across similar sites.
"""

import os
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Union, Tuple
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Import existing functions from other modules
from protein_explorer.analysis.kinase_predictor import (
    load_kinase_scores, get_site_kinase_scores, predict_kinases,
    get_known_kinase_info, categorize_kinases_by_family
)

from protein_explorer.analysis.phospho_analyzer import (
    find_structural_matches, get_phosphosite_data
)

from protein_explorer.analysis.sequence_analyzer import (
    find_sequence_matches
)

def get_similar_sites(site_id: str, 
                      uniprot_id: Optional[str] = None, 
                      site_data: Optional[Dict] = None,
                      similarity_threshold: float = 0.6, 
                      rmsd_threshold: float = 3.0) -> List[str]:
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
    # Parse site_id to get uniprot_id and site_number if not provided
    if not uniprot_id:
        parts = site_id.split('_')
        if len(parts) >= 2:
            uniprot_id = parts[0]
        else:
            logger.warning(f"Could not extract UniProt ID from site_id: {site_id}")
            return [site_id]  # Return just the query site if we can't extract UniProt ID
    
    # Parse site information if not provided
    site_number = None
    site_type = None  # We'll determine this from site_data if possible
    
    if not site_data:
        parts = site_id.split('_')
        if len(parts) >= 2:
            # Sites are in format UniProtID_ResidueNumber
            try:
                site_number = int(parts[1])
                # Set a default site_type for structural matching
                site_type = 'S'  # Default to Serine as most common
            except ValueError:
                logger.warning(f"Could not parse site number from {parts[1]}")
                return [site_id]
    else:
        # Extract from site_data
        site_number = site_data.get('resno')
        # Try to get site_type from site_data
        if 'siteType' in site_data:
            site_type = site_data['siteType']
        elif 'site' in site_data and site_data['site'] and len(site_data['site']) > 0:
            site_type = site_data['site'][0]  # First character
        else:
            site_type = 'S'  # Default
    
    if not site_number:
        logger.warning(f"Could not extract site number from site_id: {site_id}")
        return [site_id]  # Return just the query site if we can't extract site number
    
    # Ensure site_type is one of S, T, Y for structural analysis
    if not site_type or site_type not in 'STY':
        site_type = 'S'  # Default to Serine if no valid type
    
    logger.info(f"Finding similar sites for {site_id} (UniProt: {uniprot_id}, Site: {site_type}{site_number})")
    
    # Get sequence-similar sites
    sequence_similar_sites = []
    try:
        # This function handles the site_id format correctly
        seq_matches = find_sequence_matches(site_id, min_similarity=similarity_threshold)
        sequence_similar_sites = [match['target_id'] for match in seq_matches 
                                if match.get('similarity', 0) >= similarity_threshold]
        logger.info(f"Found {len(sequence_similar_sites)} sequence-similar sites")
    except Exception as e:
        logger.error(f"Error finding sequence-similar sites: {e}")
    
    # Get structurally-similar sites
    structural_similar_sites = []
    try:
        # For structural matching, we need site_type + number format
        site_str = f"{site_type}{site_number}"
        
        # Create minimal site data if not provided
        if not site_data:
            site_data = {
                'site': site_str,
                'resno': site_number,
                'siteType': site_type
            }
        elif 'site' not in site_data or not site_data['site']:
            # Ensure site is in site_data
            site_data['site'] = site_str
        
        # Find structural matches
        struct_matches = find_structural_matches(uniprot_id, [site_data])
        
        if site_str in struct_matches:
            raw_matches = struct_matches[site_str]
            # Convert to site_ids in format UniProt_ResNo
            for match in raw_matches:
                if match.get('rmsd', 10) < rmsd_threshold:
                    # Parse target site - might be in format S123 or just 123
                    target_site = match.get('target_site', '')
                    target_uniprot = match.get('target_uniprot', '')
                    
                    if not target_uniprot:
                        continue
                        
                    # Extract just the number from target_site
                    import re
                    site_match = re.match(r'[STY]?(\d+)', target_site)
                    if site_match:
                        # Use your site ID format: UniProtID_ResidueNumber
                        target_id = f"{target_uniprot}_{site_match.group(1)}"
                        structural_similar_sites.append(target_id)
            
            logger.info(f"Found {len(structural_similar_sites)} structurally-similar sites")
        else:
            logger.warning(f"No structural matches found for site {site_str}")
    except Exception as e:
        logger.error(f"Error finding structurally-similar sites: {e}")
    
    # Combine all sites, remove duplicates, and include the query site
    all_similar_sites = list(set([site_id] + sequence_similar_sites + structural_similar_sites))
    
    logger.info(f"Total similar sites for aggregation: {len(all_similar_sites)}")
    return all_similar_sites

def get_network_kinase_scores(site_id: str, score_type: str = 'structure',
                             similarity_threshold: float = 0.6, 
                             rmsd_threshold: float = 3.0) -> Dict[str, Dict[str, float]]:
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
    # Get similar sites
    similar_sites = get_similar_sites(site_id, 
                                     similarity_threshold=similarity_threshold, 
                                     rmsd_threshold=rmsd_threshold)
    
    # Get scores for all sites
    all_scores = {}
    for site in similar_sites:
        site_scores = get_site_kinase_scores(site, score_type)
        if site_scores and 'scores' in site_scores:
            all_scores[site] = site_scores['scores']
    
    return all_scores

def compute_aggregated_kinase_scores(site_id: str, score_type: str = 'structure',
                                   similarity_threshold: float = 0.6,
                                   rmsd_threshold: float = 3.0) -> List[Dict]:
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

def predict_kinases_network(site_id: str, top_n: int = 5, score_type: str = 'structure',
                           similarity_threshold: float = 0.6, 
                           rmsd_threshold: float = 3.0) -> List[Dict]:
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

def get_network_heatmap_data(site_id: str, top_n: int = 10, score_type: str = 'structure',
                          similarity_threshold: float = 0.6, 
                          rmsd_threshold: float = 3.0) -> Dict:
    """
    Get data for heatmap visualization of kinase scores across a network of similar sites.
    Enhanced with better error handling and data validation.
    
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
    
    # Log data for debugging
    logger.info(f"Creating heatmap with {len(heatmap_data['sites'])} sites and {len(top_kinases)} kinases")
    
    # Add scores for each site and kinase
    for site in network_scores.keys():
        site_scores = network_scores[site]
        
        for kinase in top_kinases:
            score = site_scores.get(kinase, 0)
            
            # Always add data point even if score is 0 to ensure complete heatmap
            heatmap_data['scores'].append({
                'site': site,
                'kinase': kinase,
                'score': float(score)
            })
    
    # Verify data integrity
    logger.info(f"Heatmap data has {len(heatmap_data['scores'])} score entries")
    
    # Ensure we have data for all site-kinase combinations
    expected_entries = len(heatmap_data['sites']) * len(heatmap_data['kinases'])
    if len(heatmap_data['scores']) != expected_entries:
        logger.warning(f"Heatmap data incomplete: expected {expected_entries} entries, got {len(heatmap_data['scores'])}")
        
        # Fill in missing entries with zeros
        existing_pairs = set((entry['site'], entry['kinase']) for entry in heatmap_data['scores'])
        
        for site in heatmap_data['sites']:
            for kinase in heatmap_data['kinases']:
                if (site, kinase) not in existing_pairs:
                    logger.info(f"Adding missing entry for site {site}, kinase {kinase}")
                    heatmap_data['scores'].append({
                        'site': site,
                        'kinase': kinase,
                        'score': 0.0
                    })
    
    return heatmap_data

def get_network_kinase_comparison(site_id: str, score_types: List[str] = ['structure', 'sequence'],
                               top_n: int = 5,
                               similarity_threshold: float = 0.6, 
                               rmsd_threshold: float = 3.0) -> Dict:
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

def get_kinase_family_distribution_network(site_id: str, score_type: str = 'structure',
                                         similarity_threshold: float = 0.6, 
                                         rmsd_threshold: float = 3.0) -> Dict:
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
    
    # Use existing function to categorize by family
    return categorize_kinases_by_family(kinases_for_categorization)