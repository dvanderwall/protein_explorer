"""
Functions for retrieving protein data from online databases.
"""

import os
import json
import requests
import logging
from typing import Dict, List, Optional, Union
import time

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# API endpoints
UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api"
STRING_API = "https://string-db.org/api"

# Cache directory - now with a function to ensure it exists
def get_cache_dir():
    """
    Get the cache directory path and ensure it exists.
    
    Returns:
        Path to the cache directory
    """
    cache_dir = os.path.expanduser("~/.protein_explorer/cache")
    
    # Create directory if it doesn't exist
    if not os.path.exists(cache_dir):
        try:
            print(f"Creating cache directory: {cache_dir}")
            os.makedirs(cache_dir, exist_ok=True)
        except Exception as e:
            # If we can't create the directory, use a temporary directory
            import tempfile
            cache_dir = os.path.join(tempfile.gettempdir(), "protein_explorer_cache")
            print(f"Failed to create cache at {cache_dir}, using temporary directory instead: {cache_dir}")
            os.makedirs(cache_dir, exist_ok=True)
    
    return cache_dir

def save_to_cache(filename, data):
    """
    Save data to the cache file.
    
    Args:
        filename: Name of the cache file
        data: Data to save (must be JSON serializable)
    """
    cache_dir = get_cache_dir()
    cache_file = os.path.join(cache_dir, filename)
    
    try:
        with open(cache_file, 'w') as f:
            json.dump(data, f)
        return True
    except Exception as e:
        logger.error(f"Error saving to cache: {e}")
        return False


###################################################################################################################
### CAN SOMEONE TAKE THE STRUCTURAL SIMILARITY DATABASE AND SET UP MONGODB/SEQUEL DATABASE FOR PROTEIN EXPLORER ###
##################################### CONSIDER MONGODB INTEGRATION ################################################
###################################################################################################################
def load_from_cache(filename):
    """
    Load data from the cache file.
    
    Args:
        filename: Name of the cache file
        
    Returns:
        Loaded data or None if file doesn't exist or error occurs
    """
    cache_dir = get_cache_dir()
    cache_file = os.path.join(cache_dir, filename)
    
    if not os.path.exists(cache_file):
        return None
        
    try:
        with open(cache_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Error loading from cache: {e}")
        return None

def get_uniprot_id_from_gene(gene_symbol: str, organism: str = "human") -> Optional[str]:
    """
    Convert a gene symbol to a UniProt ID for reviewed Homo sapiens entries.
    
    Args:
        gene_symbol: Gene symbol (e.g., "TP53")
        organism: Organism name (default: "human")
        
    Returns:
        Reviewed UniProt ID (in uppercase) or None if not found.
    """
    # Force gene symbol to uppercase
    gene_symbol = gene_symbol.upper()
    
    # Prepare cache key and filename
    cache_key = f"{gene_symbol}_{organism}"
    cache_filename = f"gene_{cache_key}.json"
    
    # Check cache first
    cache_data = load_from_cache(cache_filename)
    if cache_data:
        uniprot_id = cache_data.get('uniprot_id')
        if uniprot_id:
            return uniprot_id.upper()
    
    # Prepare search query
    # For human, use organism_id:9606
    if organism.lower() == "human":
        query = f"gene:{gene_symbol} AND organism_id:9606"
    else:
        query = f"gene:{gene_symbol} AND organism:{organism}"
    
    # Build the URL
    url = f"{UNIPROT_API}/search?query={query}&format=json&size=10"
    print(url)
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        print(data)
        
        # Iterate over results to find a reviewed Homo sapiens entry
        reviewed_result = None
        for r in data.get('results', []):
            organism_data = r.get('organism', {})
            entry_type = r.get('entryType', '').lower()
            if (organism_data.get('scientificName', '').lower() == 'homo sapiens' or organism_data.get('organismId') == 9606) \
               and ('reviewed' in entry_type):
                reviewed_result = r
                break
        
        if reviewed_result:
            print(reviewed_result)
            uniprot_id = reviewed_result.get('primaryAccession')
            if not uniprot_id:
                # Fallback: use first secondary accession if available
                secondary = reviewed_result.get('secondaryAccession')
                if secondary and len(secondary) > 0:
                    uniprot_id = secondary[0]
            if uniprot_id:
                uniprot_id = uniprot_id.upper()  # Ensure uppercase
                save_to_cache(cache_filename, {'uniprot_id': uniprot_id})
                return uniprot_id
            else:
                logger.warning(f"No UniProt ID found for gene {gene_symbol}")
                return None
        else:
            logger.warning(f"No reviewed Homo sapiens entry returned for gene {gene_symbol}")
            return None
    except requests.exceptions.RequestException as e:
        logger.error(f"Error converting gene symbol to UniProt ID: {e}")
        return None

def get_protein_by_id(uniprot_id: Optional[str] = None, 
                      gene_symbol: Optional[str] = None, 
                      organism: str = "human") -> Dict:
    """
    Retrieve protein data from UniProt by either UniProt ID or gene symbol.
    
    Args:
        uniprot_id: UniProt ID (e.g., "P53_HUMAN")
        gene_symbol: Gene symbol (e.g., "TP53")
        organism: Organism name (default: "human")
        
    Returns:
        Dictionary with protein metadata
    """
    if not uniprot_id and not gene_symbol:
        raise ValueError("Either UniProt ID or gene symbol must be provided")
        
    # If only gene symbol is provided, convert to UniProt ID
    if gene_symbol and not uniprot_id:
        # Force gene symbol to uppercase
        gene_symbol = gene_symbol.upper()
        uniprot_id = get_uniprot_id_from_gene(gene_symbol, organism)
        if not uniprot_id:
            raise ValueError(f"Could not find UniProt ID for gene {gene_symbol}")
    else:
        # If uniprot_id is provided directly, force it to uppercase
        uniprot_id = uniprot_id.upper()
    
    # Cache filename
    cache_filename = f"uniprot_{uniprot_id}.json"
    
    # Check cache
    metadata = load_from_cache(cache_filename)
    if not metadata:
        # Get protein metadata from UniProt
        url = f"{UNIPROT_API}/{uniprot_id}.json"
        try:
            response = requests.get(url)
            response.raise_for_status()
            metadata = response.json()
            
            # Cache the result
            save_to_cache(cache_filename, metadata)
                
        except requests.exceptions.RequestException as e:
            logger.error(f"Error retrieving UniProt data: {e}")
            raise ValueError(f"Failed to retrieve data for {uniprot_id}")
    
    # Check if AlphaFold structure exists
    has_structure = check_alphafold_exists(uniprot_id)
    
    # Prepare result
    result = {
        "uniprot_id": uniprot_id,
        "metadata": metadata,
        "has_structure": has_structure
    }

    # Add gene symbol if it wasn't provided
    if not gene_symbol:
        try:
            gene_names = metadata.get("genes", [])
            if gene_names and len(gene_names) > 0 and "geneName" in gene_names[0]:
                result["gene_symbol"] = gene_names[0]["geneName"]["value"]
        except (KeyError, IndexError):
            pass
    else:
        result["gene_symbol"] = gene_symbol
        
    return result

def check_alphafold_exists(uniprot_id: str, force_check: bool = False) -> bool:
    """
    Check if AlphaFold structure exists for a given UniProt ID.
    Tries multiple versions and formats to be more robust.
    
    Args:
        uniprot_id: UniProt ID
        force_check: If True, bypass cache and check directly
        
    Returns:
        Boolean indicating if structure exists
    """
    print(f"DEBUG: Checking if AlphaFold structure exists for {uniprot_id}")
    
    # Ensure cache directory exists and get cache file path
    cache_filename = f"af_exists_{uniprot_id}.json"
    cache_data = None if force_check else load_from_cache(cache_filename)
    
    # Check cache (unless force_check is True)
    if cache_data is not None:
        result = cache_data.get('exists', False)
        print(f"DEBUG: Cache indicates exists={result}")
        
        # If cache has a URL, return that information
        if 'url' in cache_data and cache_data['url']:
            return result
    
    print(f"DEBUG: Checking directly (bypass cache: {force_check})")
    
    # List of possible URL patterns to try
    urls_to_try = [
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb",
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v3.pdb",
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v2.pdb",
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v1.pdb"
    ]
    
    # Try each URL pattern
    for url in urls_to_try:
        print(f"DEBUG: Checking URL: {url}")
        try:
            response = requests.head(url, timeout=5)
            print(f"DEBUG: Response status code: {response.status_code}")
            
            if response.status_code == 200:
                # Cache the result
                print(f"DEBUG: Caching result: exists=True, url={url}")
                save_to_cache(cache_filename, {'exists': True, 'url': url})
                return True
        except requests.exceptions.RequestException as e:
            print(f"DEBUG: Request exception for {url}: {e}")
    
    # If none of the URLs worked, cache the negative result
    print(f"DEBUG: Caching result: exists=False")
    save_to_cache(cache_filename, {'exists': False, 'url': None})
    return False

def get_alphafold_structure(uniprot_id: str) -> Optional[str]:
    """
    Download the AlphaFold structure for a given UniProt ID.
    Tries multiple versions and formats for robustness.
    
    Args:
        uniprot_id: UniProt ID
        
    Returns:
        PDB format structure as string, or None if not available
    """
    print(f"DEBUG: Getting AlphaFold structure for {uniprot_id}")
    
    # Ensure cache directory exists and get cache file path
    cache_filename = f"alphafold_{uniprot_id}.pdb"
    cache_dir = get_cache_dir()
    cache_file = os.path.join(cache_dir, cache_filename)
    
    # Check cache
    if os.path.exists(cache_file):
        print(f"DEBUG: Structure cache file exists, reading from cache")
        try:
            with open(cache_file, 'r') as f:
                structure = f.read()
                return structure
        except Exception as e:
            print(f"DEBUG: Error reading from cache: {e}")
            # Continue to download if cache read fails
    
    print(f"DEBUG: No structure cache, downloading")
    
    # Check if structure exists and get the URL
    cache_info_filename = f"af_exists_{uniprot_id}.json"
    cache_data = load_from_cache(cache_info_filename)
    
    if cache_data:
        if cache_data.get('exists') and cache_data.get('url'):
            url = cache_data['url']
            print(f"DEBUG: Using cached URL: {url}")
        else:
            # Structure doesn't exist according to cache
            print(f"DEBUG: Structure does not exist according to cache")
            return None
    else:
        # No cache info, need to check if structure exists
        if not check_alphafold_exists(uniprot_id):
            print(f"DEBUG: Structure not found for {uniprot_id}")
            return None
        
        # Read the updated cache to get the URL
        cache_data = load_from_cache(cache_info_filename)
        url = cache_data.get('url') if cache_data else None
            
        if not url:
            # Still no URL, use default v4 as a last resort
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    
    # Download from AlphaFold
    print(f"DEBUG: Downloading from URL: {url}")
    
    try:
        response = requests.get(url)
        print(f"DEBUG: Response status code: {response.status_code}")
        
        if response.status_code != 200:
            print(f"DEBUG: Failed to download structure, status code: {response.status_code}")
            return None
            
        structure = response.text
        print(f"DEBUG: Downloaded {len(structure)} characters")
        
        # Cache the result
        print(f"DEBUG: Caching downloaded structure")
        try:
            with open(cache_file, 'w') as f:
                f.write(structure)
        except Exception as e:
            print(f"DEBUG: Error caching structure: {e}")
            # Continue even if caching fails
            
        return structure
    except requests.exceptions.RequestException as e:
        print(f"DEBUG: Request exception: {e}")
        logger.error(f"Error downloading AlphaFold structure: {e}")
        return None
    

###################################################################################################################
#################### STRINGS DATABASE BIOPLEX? ####################################################################
###################################################################################################################
def get_protein_interactions(uniprot_id: str, 
                           confidence_score: float = 0.7, 
                           limit: int = 100,
                           organism_id: int = 9606) -> Dict:
    """
    Retrieve protein-protein interactions from the STRING database.
    
    Args:
        uniprot_id: UniProt ID
        confidence_score: Minimum confidence score (0.0 to 1.0)
        limit: Maximum number of interactions to retrieve
        organism_id: NCBI taxonomy ID (default: 9606 for human)
        
    Returns:
        Dictionary of interacting proteins and confidence scores
    """
    cache_filename = f"string_{uniprot_id}_{confidence_score}.json"
    
    # Check cache
    interactions = load_from_cache(cache_filename)
    if interactions:
        return interactions
    
    # Prepare API request
    url = f"{STRING_API}/json/network"
    params = {
        "identifiers": uniprot_id,
        "species": organism_id,
        "caller_identity": "protein_explorer",
        "required_score": int(confidence_score * 1000),
        "limit": limit
    }
    
    try:
        response = requests.post(url, data=params)
        response.raise_for_status()
        data = response.json()
        
        # Process results
        interactions = {}
        for edge in data.get("edges", []):
            if edge["from"] == uniprot_id:
                target = edge["to"]
            else:
                target = edge["from"]
                
            score = edge["score"] / 1000.0  # Convert to 0.0-1.0 range
            interactions[target] = score
            
        # Cache the result
        save_to_cache(cache_filename, interactions)
            
        return interactions
    except requests.exceptions.RequestException as e:
        logger.error(f"Error retrieving protein interactions: {e}")
        return {}