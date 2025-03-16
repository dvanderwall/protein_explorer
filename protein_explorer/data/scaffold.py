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

# Cache directory
CACHE_DIR = os.path.expanduser("~/.protein_explorer/cache")
os.makedirs(CACHE_DIR, exist_ok=True)

def get_uniprot_id_from_gene(gene_symbol: str, organism: str = "human") -> Optional[str]:
    """
    Convert a gene symbol to a UniProt ID.
    
    Args:
        gene_symbol: Gene symbol (e.g., "TP53")
        organism: Organism name (default: "human")
        
    Returns:
        UniProt ID or None if not found
    """
    # Cache key
    cache_key = f"{gene_symbol}_{organism}"
    cache_file = os.path.join(CACHE_DIR, f"gene_{cache_key}.json")
    
    # Check cache
    if os.path.exists(cache_file):
        with open(cache_file, 'r') as f:
            return json.load(f).get('uniprot_id')
    
    # Prepare search query
    if organism.lower() == "human":
        query = f"gene:{gene_symbol} AND organism_id:9606"
    else:
        query = f"gene:{gene_symbol} AND organism:{organism}"
    
    # Make API request
    url = f"{UNIPROT_API}/search?query={query}&format=json&size=1"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        
        if data.get('results') and len(data['results']) > 0:
            uniprot_id = data['results'][0]['primaryAccession']
            
            # Cache the result
            with open(cache_file, 'w') as f:
                json.dump({'uniprot_id': uniprot_id}, f)
                
            return uniprot_id
        else:
            logger.warning(f"No UniProt ID found for gene {gene_symbol}")
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
        uniprot_id = get_uniprot_id_from_gene(gene_symbol, organism)
        if not uniprot_id:
            raise ValueError(f"Could not find UniProt ID for gene {gene_symbol}")
    
    # Cache key and file
    cache_file = os.path.join(CACHE_DIR, f"uniprot_{uniprot_id}.json")
    
    # Check cache
    if os.path.exists(cache_file):
        with open(cache_file, 'r') as f:
            metadata = json.load(f)
    else:
        # Get protein metadata from UniProt
        url = f"{UNIPROT_API}/{uniprot_id}.json"
        try:
            response = requests.get(url)
            response.raise_for_status()
            metadata = response.json()
            
            # Cache the result
            with open(cache_file, 'w') as f:
                json.dump(metadata, f)
                
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
        # Extract gene name from UniProt metadata
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
    
    Args:
        uniprot_id: UniProt ID
        force_check: If True, bypass cache and check directly
        
    Returns:
        Boolean indicating if structure exists
    """
    print(f"DEBUG: Checking if AlphaFold structure exists for {uniprot_id}")
    
    cache_file = os.path.join(CACHE_DIR, f"af_exists_{uniprot_id}.json")
    print(f"DEBUG: Cache file path: {cache_file}")
    
    # Check cache (unless force_check is True)
    if not force_check and os.path.exists(cache_file):
        print(f"DEBUG: Cache file exists, reading from cache")
        with open(cache_file, 'r') as f:
            result = json.load(f).get('exists', False)
            print(f"DEBUG: Cache indicates exists={result}")
            return result
    
    print(f"DEBUG: Checking directly (bypass cache: {force_check})")
    # Use direct URL for checking - try v4 first
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    print(f"DEBUG: Checking URL: {url}")
    
    try:
        response = requests.head(url, timeout=5)
        print(f"DEBUG: Response status code: {response.status_code}")
        exists = response.status_code == 200
        
        # Cache the result
        print(f"DEBUG: Caching result: exists={exists}")
        with open(cache_file, 'w') as f:
            json.dump({'exists': exists, 'url': url if exists else None}, f)
            
        return exists
    except requests.exceptions.RequestException as e:
        print(f"DEBUG: Request exception: {e}")
        return False

def get_alphafold_structure(uniprot_id: str) -> Optional[str]:
    """
    Download the AlphaFold structure for a given UniProt ID.
    
    Args:
        uniprot_id: UniProt ID
        
    Returns:
        PDB format structure as string, or None if not available
    """
    print(f"DEBUG: Getting AlphaFold structure for {uniprot_id}")
    
    cache_file = os.path.join(CACHE_DIR, f"alphafold_{uniprot_id}.pdb")
    print(f"DEBUG: Structure cache file path: {cache_file}")
    
    # Check cache
    if os.path.exists(cache_file):
        print(f"DEBUG: Structure cache file exists, reading from cache")
        with open(cache_file, 'r') as f:
            structure = f.read()
            print(f"DEBUG: Read {len(structure)} characters from cache")
            return structure
    
    print(f"DEBUG: No structure cache, downloading")
    # Download from AlphaFold
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    print(f"DEBUG: Downloading from URL: {url}")
    
    try:
        response = requests.get(url)
        print(f"DEBUG: Response status code: {response.status_code}")
        response.raise_for_status()
        structure = response.text
        print(f"DEBUG: Downloaded {len(structure)} characters")
        
        # Cache the result
        print(f"DEBUG: Caching downloaded structure")
        with open(cache_file, 'w') as f:
            f.write(structure)
            
        return structure
    except requests.exceptions.RequestException as e:
        print(f"DEBUG: Request exception: {e}")
        logger.error(f"Error downloading AlphaFold structure: {e}")
        return None

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
    cache_file = os.path.join(CACHE_DIR, f"string_{uniprot_id}_{confidence_score}.json")
    
    # Check cache
    if os.path.exists(cache_file):
        with open(cache_file, 'r') as f:
            return json.load(f)
    
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
        with open(cache_file, 'w') as f:
            json.dump(interactions, f)
            
        return interactions
    except requests.exceptions.RequestException as e:
        logger.error(f"Error retrieving protein interactions: {e}")
        return {}

def batch_retrieve_proteins(id_list: List[str], 
                          id_type: str = "uniprot",
                          organism: str = "human") -> Dict[str, Dict]:
    """
    Batch retrieve protein data for a list of IDs.
    
    Args:
        id_list: List of UniProt IDs or gene symbols
        id_type: Type of IDs provided ("uniprot" or "gene")
        organism: Organism name (default: "human")
        
    Returns:
        Dictionary of protein data keyed by ID
    """
    results = {}
    
    for identifier in id_list:
        try:
            if id_type.lower() == "uniprot":
                protein_data = get_protein_by_id(uniprot_id=identifier, organism=organism)
            else:
                protein_data = get_protein_by_id(gene_symbol=identifier, organism=organism)
                
            # Use the original ID as key
            results[identifier] = protein_data
            
            # Add a short delay to avoid overwhelming the API
            time.sleep(0.5)
            
        except Exception as e:
            logger.error(f"Error retrieving data for {identifier}: {e}")
            results[identifier] = {"error": str(e)}
            
    return results