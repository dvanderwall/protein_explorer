"""
Functions for analyzing potential phosphorylation sites in proteins.
"""

import numpy as np
import requests
import json
import os
from typing import Dict, List, Optional, Tuple
import requests

def analyze_phosphosites(sequence: str, 
                        structure_data: str, 
                        uniprot_id: str = None) -> List[Dict]:
    """
    Analyze potential phosphorylation sites (S, T, Y) in a protein sequence.
    
    Args:
        sequence: Protein sequence string
        structure_data: PDB format data as string
        uniprot_id: UniProt ID for checking PhosphositePlus (optional)
        
    Returns:
        List of dictionaries with phosphosite information
    """
    from Bio.PDB import PDBParser, NeighborSearch, Selection
    import io
    
    # Initialize known phosphosites dictionary
    known_sites = {}
    
    # Attempt to fetch known phosphosites if UniProt ID is provided
    if uniprot_id:
        try:
            # Use updated UniProt REST API URL
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
            
            # Perform the request
            response = requests.get(url)
            
            # Ensure successful response
            if response.status_code == 200:
                # Split the response into lines
                lines = response.text.split('\n')
                
                # Iterate through lines
                for i in range(len(lines)):
                    # Look for MOD_RES lines
                    if lines[i].startswith('FT   MOD_RES'):
                        # Extract site number
                        parts = lines[i].split()
                        
                        if len(parts) >= 3:
                            current_site = int(parts[2])
                            
                            # Check next line for phosphorylation
                            if i + 1 < len(lines):
                                next_line = lines[i + 1]
                                
                                # Check for phosphorylation types
                                if 'Phosphothreonine' in next_line:
                                    known_sites[f'T{current_site}'] = True
                                elif 'Phosphoserine' in next_line:
                                    known_sites[f'S{current_site}'] = True
                                elif 'Phosphotyrosine' in next_line:
                                    known_sites[f'Y{current_site}'] = True
        except Exception as e:
            print(f"Error fetching known phosphosites for {uniprot_id}: {e}")
    
    # Parse PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", io.StringIO(structure_data))
    
    # Extract coordinates and B-factors (PLDDT scores in AlphaFold)
    atoms = list(structure.get_atoms())
    ns = NeighborSearch(atoms)
    
    # Get B-factors (PLDDT) by residue
    plddt_by_residue = {}
    for residue in structure.get_residues():
        b_factors = [atom.get_bfactor() for atom in residue]
        if b_factors:
            plddt_by_residue[residue.get_id()[1]] = np.mean(b_factors)
    
    # Find all S, T, Y residues
    phosphosites = []
    for i, aa in enumerate(sequence):
        if aa in ['S', 'T', 'Y']:
            resno = i + 1  # 1-based residue numbering
            
            # Extract motif (-7 to +7)
            motif_start = max(0, i - 7)
            motif_end = min(len(sequence), i + 8)
            motif = sequence[motif_start:motif_end]
            
            # Calculate mean pLDDT for the motif
            motif_resno_range = range(motif_start + 1, motif_end + 1)
            motif_plddts = [plddt_by_residue.get(j, 0) for j in motif_resno_range]
            mean_plddt = np.mean(motif_plddts) if motif_plddts else 0
            
            # Count residues within 10Å
            nearby_count = 0
            try:
                # Get center atom for the residue
                for res in structure.get_residues():
                    if res.get_id()[1] == resno:
                        center_atom = res['CA'] if 'CA' in res else next(res.get_atoms())
                        nearby_atoms = ns.search(center_atom.get_coord(), 10)  # 10Å radius
                        
                        # Use a set to automatically track unique residues
                        nearby_residues = set()
                        for atom in nearby_atoms:
                            parent = atom.get_parent()
                            full_id = parent.get_id()
                            res_id = full_id[1]
                            #print(f"Atom: {atom.get_name()}, Parent: {parent}, Full ID: {full_id}, Res ID: {res_id}")
                            if res_id != resno:
                                nearby_residues.add(res_id)
                                
                        nearby_count = len(nearby_residues)
                        break
            except Exception as e:
                print(f"Error calculating neighbors for {aa}{resno}: {str(e)}")
                nearby_count = 0
            
            # Check if known in PhosphositePlus
            site_key = f"{aa}{resno}"
            is_known = site_key in known_sites
            
            # Add to results
            phosphosites.append({
                'site': site_key,
                'resno': resno,
                'motif': motif,
                'mean_plddt': f"{mean_plddt:.1f}",
                'nearby_count': nearby_count,
                'is_known': is_known
            })
    
    return phosphosites