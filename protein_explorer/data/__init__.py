"""
Data scaffolding and processing modules for protein data.
"""

from protein_explorer.data.scaffold import (
    get_protein_by_id,
    get_alphafold_structure
)

from protein_explorer.data.processor import (
    parse_pdb_structure,
    extract_coordinates,
    parse_interaction_data
)