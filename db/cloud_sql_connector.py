"""
Cloud SQL connector for KinoPlex.
Provides database connection and query functions for the phosphosite database.
"""

import os
import logging
import pandas as pd
import sqlalchemy
from sqlalchemy import create_engine, text
from sqlalchemy.pool import QueuePool
import pymysql
import json
from typing import Dict, List, Optional, Tuple, Union
from datetime import datetime, timedelta
import time
import threading
import functools

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Local memory cache for query results
_CACHE = {}
_CACHE_LOCK = threading.Lock()
_CACHE_EXPIRY = {}  # Store expiry times for cache entries
DEFAULT_CACHE_TTL = 3600  # Default TTL of 1 hour

# Configuration
class DBConfig:
    """Database configuration for Cloud SQL."""
    INSTANCE_CONNECTION_NAME = os.environ.get('CLOUD_SQL_CONNECTION_NAME', 'future-alcove-454817-e6:us-east4:kinoplex-db')
    DB_USER = os.environ.get('DB_USER', 'root')  # Replace with your database user if not root
    DB_PASS = os.environ.get('DB_PASS', '@Bismark6')  # Replace with your actual password
    DB_NAME = os.environ.get('DB_NAME', 'kinoplex-db')  # This should be the database name, not instance name
    DB_HOST = os.environ.get('DB_HOST', '35.245.113.195')  # Use the public IP for direct connections
    DB_PORT = os.environ.get('DB_PORT', '3306')
    
    # Use SSL for secure connection if cert path provided
    SSL_CA = os.environ.get('SSL_CA', '')
    SSL_CERT = os.environ.get('SSL_CERT', '')
    SSL_KEY = os.environ.get('SSL_KEY', '')
    
    # Connection pool config
    POOL_SIZE = int(os.environ.get('DB_POOL_SIZE', '5'))
    MAX_OVERFLOW = int(os.environ.get('DB_MAX_OVERFLOW', '10'))
    POOL_TIMEOUT = int(os.environ.get('DB_POOL_TIMEOUT', '30'))
    POOL_RECYCLE = int(os.environ.get('DB_POOL_RECYCLE', '1800'))  # 30 minutes

# Initialize the engine globally so we can reuse connections
_db_engine = None

def get_db_engine():
    """Get or create a SQLAlchemy engine for database connections."""
    global _db_engine
    
    if _db_engine is not None:
        return _db_engine
    
    # Check if running on Windows (which doesn't support Unix sockets)
    import platform
    import urllib.parse
    is_windows = platform.system() == 'Windows'
    
    # URL-encode the username and password to handle special characters
    encoded_user = urllib.parse.quote_plus(DBConfig.DB_USER)
    encoded_pass = urllib.parse.quote_plus(DBConfig.DB_PASS)  # This will encode @ as %40
    
    # Build connection string based on environment (Cloud SQL vs direct)
    if DBConfig.INSTANCE_CONNECTION_NAME and not is_windows:
        # Running on App Engine or other Google Cloud service (non-Windows)
        logger.info(f"Connecting to Cloud SQL instance: {DBConfig.INSTANCE_CONNECTION_NAME}")
        
        # Use Unix socket if on Google Cloud with encoded credentials
        connection_string = (
            f"mysql+pymysql://{encoded_user}:{encoded_pass}@/{DBConfig.DB_NAME}"
            f"?unix_socket=/cloudsql/{DBConfig.INSTANCE_CONNECTION_NAME}"
        )
    else:
        # Running locally, on Windows, or elsewhere with direct connection
        logger.info(f"Connecting to MySQL database at: {DBConfig.DB_HOST}:{DBConfig.DB_PORT}")
        
        # Use TCP connection with properly encoded credentials
        connection_string = (
            f"mysql+pymysql://{encoded_user}:{encoded_pass}@{DBConfig.DB_HOST}:{DBConfig.DB_PORT}/{DBConfig.DB_NAME}"
        )
        
        # Add SSL if configured
        if DBConfig.SSL_CA and DBConfig.SSL_CERT and DBConfig.SSL_KEY:
            ssl_args = {
                'ssl_ca': DBConfig.SSL_CA,
                'ssl_cert': DBConfig.SSL_CERT,
                'ssl_key': DBConfig.SSL_KEY
            }
            connection_string += '?' + '&'.join([f"{k}={v}" for k, v in ssl_args.items()])
    
    # Create the SQLAlchemy engine with connection pooling
    try:
        # Log the connection string with password masked for security
        masked_connection = connection_string.replace(encoded_pass, "********")
        logger.info(f"Creating database engine with connection string: {masked_connection}")
        
        _db_engine = create_engine(
            connection_string,
            poolclass=QueuePool,
            pool_size=DBConfig.POOL_SIZE,
            max_overflow=DBConfig.MAX_OVERFLOW,
            pool_timeout=DBConfig.POOL_TIMEOUT,
            pool_recycle=DBConfig.POOL_RECYCLE
        )
        
        # Test the connection
        with _db_engine.connect() as conn:
            conn.execute(text("SELECT 1"))
            logger.info("Database connection test successful")
            
        return _db_engine
    except Exception as e:
        import traceback
        logger.error(f"Error creating database engine: {e}")
        logger.error(traceback.format_exc())
        raise

def cache_result(ttl=DEFAULT_CACHE_TTL):
    """
    Decorator to cache function results in memory.
    
    Args:
        ttl: Time to live in seconds for cache entries
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Create a cache key from function name and arguments
            key_parts = [func.__name__]
            key_parts.extend([str(arg) for arg in args])
            key_parts.extend([f"{k}={v}" for k, v in sorted(kwargs.items())])
            cache_key = ":".join(key_parts)
            
            # Check if we have a cached result
            with _CACHE_LOCK:
                if cache_key in _CACHE:
                    # Check if cache has expired
                    if _CACHE_EXPIRY.get(cache_key, 0) > time.time():
                        logger.debug(f"Cache hit for {cache_key}")
                        return _CACHE[cache_key]
                    else:
                        # Expired, remove from cache
                        logger.debug(f"Cache expired for {cache_key}")
                        _CACHE.pop(cache_key, None)
                        _CACHE_EXPIRY.pop(cache_key, None)
            
            # Execute the function and cache the result
            result = func(*args, **kwargs)
            
            with _CACHE_LOCK:
                _CACHE[cache_key] = result
                _CACHE_EXPIRY[cache_key] = time.time() + ttl
            
            return result
        return wrapper
    return decorator

def clear_cache():
    """Clear the entire cache."""
    with _CACHE_LOCK:
        _CACHE.clear()
        _CACHE_EXPIRY.clear()
    logger.info("Cache cleared")

def execute_query(query: str, params: Dict = None) -> pd.DataFrame:
    """
    Execute a SQL query and return results as a pandas DataFrame.
    
    Args:
        query: SQL query string
        params: Parameters for the query
        
    Returns:
        Pandas DataFrame with query results
    """
    engine = get_db_engine()
    
    try:
        # Execute query with parameters if provided
        if params:
            df = pd.read_sql_query(text(query), engine, params=params)
        else:
            df = pd.read_sql_query(query, engine)
        
        return df
    except Exception as e:
        logger.error(f"Database query error: {e}")
        logger.error(f"Query: {query}")
        if params:
            logger.error(f"Params: {params}")
        raise


        logger.error(f"Error getting sequence matches: {e}")
        return pd.DataFrame()

@cache_result(ttl=3600)  # Cache for 1 hour
def get_phosphosite_data(site_id: str) -> Dict:
    """
    Get supplementary data for a phosphosite from the database.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        
    Returns:
        Dictionary with supplementary data
    """
    # Parse site_id to get uniprot_id and residue number
    parts = site_id.split('_')
    if len(parts) < 2:
        logger.warning(f"Invalid site_id format: {site_id}")
        return {}
    
    uniprot_id = parts[0]
    residue_number = parts[1]
    
    query = """
    SELECT * 
    FROM Phosphosite_Supplementary_Data 
    WHERE uniprot_id = :uniprot_id 
    AND Residue_Number = :residue_number
    LIMIT 1
    """
    
    try:
        df = execute_query(query, {'uniprot_id': uniprot_id, 'residue_number': residue_number})
        if len(df) > 0:
            # Convert first row to dictionary
            return df.iloc[0].to_dict()
        else:
            # Try a more flexible search if no exact match
            fallback_query = """
            SELECT * 
            FROM Phosphosite_Supplementary_Data 
            WHERE uniprot_id = :uniprot_id 
            AND (PhosphositeID = :site_id OR MOD_RSD LIKE :mod_residue)
            LIMIT 1
            """
            
            mod_residue_pattern = f"%{residue_number}%"
            df = execute_query(fallback_query, {
                'uniprot_id': uniprot_id, 
                'site_id': site_id,
                'mod_residue': mod_residue_pattern
            })
            
            if len(df) > 0:
                return df.iloc[0].to_dict()
            
            logger.warning(f"No phosphosite data found for {site_id}")
            return {}
    except Exception as e:
        logger.error(f"Error getting phosphosite data: {e}")
        return {}

@cache_result(ttl=3600)  # Cache for 1 hour
def get_kinase_scores(site_id: str, score_type: str = 'structure') -> Dict:
    """
    Get kinase prediction scores for a site from the database.
    
    Args:
        site_id: Site ID in format 'UniProtID_ResidueNumber'
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with kinase scores
    """
    # Use the appropriate table based on score_type
    table = f"{'Structure' if score_type.lower() == 'structure' else 'Sequence'}_Kinase_Scores"
    
    query = f"""
    SELECT *
    FROM {table}
    WHERE node = :site_id
    LIMIT 1
    """
    
    try:
        df = execute_query(query, {'site_id': site_id})
        if len(df) > 0:
            # Get all kinase score columns - these are all columns except 'node' and 'label'
            all_columns = df.columns.tolist()
            kinase_columns = [col for col in all_columns if col not in ['node', 'label']]
            
            # Extract scores into a dictionary
            scores_dict = {}
            for kinase in kinase_columns:
                try:
                    score_value = df.iloc[0][kinase]
                    if pd.notna(score_value):  # Skip NaN values
                        scores_dict[kinase] = float(score_value)
                except (ValueError, TypeError):
                    # Skip values that can't be converted to float
                    logger.debug(f"Skipping non-numeric value for {kinase}")
            
            # Get known kinase from supplementary data if available
            # This assumes the first KINASE_* field with a value is the known kinase
            phosphosite_data = get_phosphosite_data(site_id)
            known_kinase = 'unlabeled'
            
            for i in range(1, 23):  # KINASE_1 through KINASE_22
                kinase_field = f"KINASE_{i}"
                if kinase_field in phosphosite_data and phosphosite_data[kinase_field]:
                    known_kinase = phosphosite_data[kinase_field]
                    break
            
            return {
                'known_kinase': known_kinase,
                'scores': scores_dict
            }
        else:
            logger.warning(f"No {score_type} kinase scores found for {site_id}")
            return {}
    except Exception as e:
        logger.error(f"Error getting kinase scores: {e}")
        return {}

@cache_result(ttl=3600)  # Cache for 1 hour
def get_all_phosphosites(uniprot_id: str) -> List[Dict]:
    """
    Get all phosphosites for a protein from the database.
    
    Args:
        uniprot_id: UniProt ID of the protein
        
    Returns:
        List of dictionaries with phosphosite information
    """
    query = """
    SELECT *
    FROM Phosphosite_Supplementary_Data
    WHERE uniprot_id = :uniprot_id
    """
    
    try:
        df = execute_query(query, {'uniprot_id': uniprot_id})
        if len(df) > 0:
            # Process the data to match expected format
            sites = []
            for _, row in df.iterrows():
                site_data = row.to_dict()
                
                # Extract residue type and number
                residue = site_data.get('Residue', '')
                residue_number = site_data.get('Residue_Number')
                
                if residue and residue_number:
                    # Create a standardized site entry
                    site_entry = {
                        'site': residue + str(int(residue_number)),
                        'resno': int(residue_number),
                        'siteType': residue,
                        'motif': site_data.get('SITE_+/-7_AA', ''),
                        'mean_plddt': site_data.get('motif_plddt'),
                        'site_plddt': site_data.get('site_plddt'),
                        'nearby_count': site_data.get('nearby_count'),
                        'surface_accessibility': site_data.get('surface_accessibility'),
                        'is_known_phosphosite': site_data.get('is_known_phosphosite', 1),
                        'polar_aa_percent': site_data.get('polar_aa_percent'),
                        'nonpolar_aa_percent': site_data.get('nonpolar_aa_percent'),
                        'acidic_aa_percent': site_data.get('acidic_aa_percent'),
                        'basic_aa_percent': site_data.get('basic_aa_percent'),
                        'full_data': site_data  # Include all original data
                    }
                    sites.append(site_entry)
            
            return sites
        else:
            logger.warning(f"No phosphosites found for {uniprot_id}")
            return []
    except Exception as e:
        logger.error(f"Error getting phosphosites: {e}")
        return []

def get_similar_sites(site_id: str, similarity_threshold: float = 0.6, rmsd_threshold: float = 3.0) -> List[str]:
    """
    Get a list of sites similar to the query site based on both
    sequence similarity and structural similarity.
    
    Args:
        site_id: The query site ID (format: UniProtID_ResidueNumber)
        similarity_threshold: Minimum sequence similarity score to include
        rmsd_threshold: Maximum RMSD value to include for structural matches
        
    Returns:
        List of site IDs (including the query site) for aggregation
    """
    # Split site_id to get uniprot_id
    parts = site_id.split('_')
    if len(parts) < 2:
        logger.warning(f"Invalid site_id format: {site_id}")
        return [site_id]
    
    uniprot_id = parts[0]
    site_number = parts[1]
    
    # Get sequence-similar sites
    sequence_query = """
    SELECT target_id
    FROM sequence_similarity
    WHERE query_id = :site_id AND similarity >= :similarity_threshold
    """
    
    # Get structurally-similar sites
    structure_query = """
    SELECT target_uniprot, target_site
    FROM structural_similarity
    WHERE query_id = :site_id AND rmsd < :rmsd_threshold
    """
    
    try:
        # Get sequence matches
        seq_df = execute_query(sequence_query, {
            'site_id': site_id, 
            'similarity_threshold': similarity_threshold
        })
        
        # Get structural matches
        struct_df = execute_query(structure_query, {
            'site_id': site_id, 
            'rmsd_threshold': rmsd_threshold
        })
        
        # Extract site IDs from results
        sequence_similar_sites = seq_df['target_id'].tolist() if not seq_df.empty else []
        
        # Process structural matches to create site IDs
        structural_similar_sites = []
        if not struct_df.empty:
            for _, row in struct_df.iterrows():
                # Extract numeric part from target_site
                import re
                site_match = re.match(r'[STY]?(\d+)', str(row['target_site']))
                if site_match:
                    site_num = site_match.group(1)
                    target_id = f"{row['target_uniprot']}_{site_num}"
                    structural_similar_sites.append(target_id)
        
        # Combine all sites, remove duplicates, and include the query site
        all_similar_sites = list(set([site_id] + sequence_similar_sites + structural_similar_sites))
        
        logger.info(f"Found {len(all_similar_sites)} similar sites for {site_id}")
        return all_similar_sites
    
    except Exception as e:
        logger.error(f"Error finding similar sites: {e}")
        return [site_id]  # Return just the query site in case of error

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
        site_scores = get_kinase_scores(site, score_type)
        if site_scores and 'scores' in site_scores:
            all_scores[site] = site_scores['scores']
    
    return all_scores

def insert_or_update_record(table: str, data: Dict, unique_key: str, additional_keys: List[str] = None):
    """
    Insert a new record or update an existing one.
    
    Args:
        table: Table name
        data: Dictionary with column names and values
        unique_key: Column name that uniquely identifies the record
        additional_keys: Additional columns to use when determining if a record exists
    """
    engine = get_db_engine()
    
    # Create WHERE clause for unique identification
    where_clauses = [f"{unique_key} = :{unique_key}"]
    where_params = {unique_key: data[unique_key]}
    
    if additional_keys:
        for key in additional_keys:
            if key in data:
                where_clauses.append(f"{key} = :{key}")
                where_params[key] = data[key]
    
    where_clause = " AND ".join(where_clauses)
    
    # Check if record exists
    check_query = f"SELECT COUNT(*) as count FROM {table} WHERE {where_clause}"
    
    try:
        with engine.connect() as conn:
            result = conn.execute(text(check_query), where_params).fetchone()
            count = result[0] if result else 0
            
            if count > 0:
                # Update existing record
                set_clauses = []
                update_params = {}
                
                for key, value in data.items():
                    set_clauses.append(f"{key} = :{key}")
                    update_params[key] = value
                
                set_clause = ", ".join(set_clauses)
                update_query = f"UPDATE {table} SET {set_clause} WHERE {where_clause}"
                
                conn.execute(text(update_query), update_params)
                conn.commit()
                logger.debug(f"Updated record in {table} where {where_clause}")
            else:
                # Insert new record
                columns = ", ".join(data.keys())
                value_params = ", ".join([f":{key}" for key in data.keys()])
                insert_query = f"INSERT INTO {table} ({columns}) VALUES ({value_params})"
                
                conn.execute(text(insert_query), data)
                conn.commit()
                logger.debug(f"Inserted new record into {table}")
    
    except Exception as e:
        logger.error(f"Error inserting/updating record: {e}")
        raise

def get_heatmap_data(site_ids: List[str], top_n: int = 10, score_type: str = 'structure') -> Dict:
    """
    Get data for heatmap visualization of kinase scores directly from database.
    
    Args:
        site_ids: List of site IDs
        top_n: Number of top kinases to include
        score_type: Type of scores - 'structure' or 'sequence'
        
    Returns:
        Dictionary with heatmap data
    """
    if not site_ids:
        return {'sites': [], 'kinases': [], 'scores': []}
    
    # Use the appropriate table based on score_type
    table = f"{'Structure' if score_type.lower() == 'structure' else 'Sequence'}_Kinase_Scores"
    
    # Create query to get scores for all sites
    placeholders = ', '.join([f"'{site}'" for site in site_ids])
    query = f"""
    SELECT * 
    FROM {table}
    WHERE node IN ({placeholders})
    """
    
    try:
        # Execute query
        df = execute_query(query)
        
        if df.empty:
            logger.warning(f"No {score_type} kinase scores found for any of the specified sites")
            return {'sites': [], 'kinases': [], 'scores': []}
        
        # Exclude non-kinase columns
        excluded_cols = ['node', 'label', 'id']
        all_columns = df.columns.tolist()
        kinase_columns = [col for col in all_columns if col not in excluded_cols]
        
        # Convert data to numeric, handling errors
        for col in kinase_columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Fill NaN values with 0
        df = df.fillna(0)
        
        # Calculate mean score for each kinase across all sites
        mean_scores = df[kinase_columns].mean()
        
        # Get top N kinases by mean score
        top_kinases = mean_scores.sort_values(ascending=False).head(top_n).index.tolist()
        
        # Prepare heatmap data
        heatmap_data = {
            'sites': df['node'].tolist(),
            'kinases': top_kinases,
            'scores': []
        }
        
        # Add scores for each site and kinase
        for idx, row in df.iterrows():
            site_id = row['node']
            for kinase in top_kinases:
                score_value = row.get(kinase, 0)
                
                # Ensure we have a valid numeric value
                try:
                    score = float(score_value) if pd.notna(score_value) else 0.0
                except (ValueError, TypeError):
                    score = 0.0
                
                heatmap_data['scores'].append({
                    'site': site_id,
                    'kinase': kinase,
                    'score': score
                })
        
        return heatmap_data
    
    except Exception as e:
        logger.error(f"Error getting heatmap data: {e}")
        return {'sites': [], 'kinases': [], 'scores': []}

def health_check() -> Dict:
    """
    Perform a health check on the database connection.
    
    Returns:
        Dictionary with health check results
    """
    start_time = time.time()
    
    try:
        engine = get_db_engine()
        with engine.connect() as conn:
            result = conn.execute(text("SELECT 1")).fetchone()
            
        elapsed_time = time.time() - start_time
        
        return {
            'status': 'healthy',
            'latency_ms': round(elapsed_time * 1000, 2),
            'timestamp': datetime.now().isoformat()
        }
    except Exception as e:
        elapsed_time = time.time() - start_time
        
        return {
            'status': 'unhealthy',
            'error': str(e),
            'latency_ms': round(elapsed_time * 1000, 2),
            'timestamp': datetime.now().isoformat()
        }