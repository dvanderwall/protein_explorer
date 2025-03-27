"""
Database integration module for Protein Explorer.

This module provides database connectivity and query functions
for accessing phosphosite data from the database.
"""

# Set version information
__version__ = '0.1.0'

# Avoid importing from db_integration here to prevent circular imports
# Functions will be accessible through the module import path
# e.g., from protein_explorer.db.db_integration import use_database