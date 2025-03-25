"""
Flask web application for the Protein Explorer.
"""

from flask import Flask, render_template, request, jsonify, redirect, url_for
import os
import sys
import logging
import networkx as nx
import requests
import numpy as np
from Bio.PDB import PDBParser, Selection, NeighborSearch
import io
import re
import json
from typing import Dict, List, Optional

# Add the parent directory to the path to import protein_explorer
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import protein_explorer as pe

# Import phosphosite analysis functions
from protein_explorer.analysis.phospho import analyze_phosphosites
from protein_explorer.analysis.phospho_analyzer import (
    preload_structural_data, 
    get_phosphosite_data,
    enhance_phosphosite, 
    enhance_structural_matches,
    find_structural_matches, 
    analyze_phosphosite_context,
    enhance_site_visualization, 
    create_comparative_motif_visualization,
    analyze_residue_distributions
)

# Import visualization functions
from protein_explorer.visualization.network import create_phosphosite_network

from protein_explorer.analysis.sequence_analyzer import (
    preload_sequence_data,
    find_sequence_matches,
    analyze_motif_conservation,
    create_sequence_network_data,
    get_motif_enrichment,
    create_sequence_motif_visualization
)

from protein_explorer.analysis.kinase_predictor import (
    load_kinase_scores, 
    get_site_kinase_scores,
    predict_kinases,
    get_heatmap_data,
    get_kinase_comparison_data,
    get_known_kinase_info,
    categorize_kinases_by_family
)

# Import the new network-based functions
from protein_explorer.analysis.network_kinase_predictor import (
    get_similar_sites, compute_aggregated_kinase_scores, predict_kinases_network,
    get_network_heatmap_data, get_network_kinase_comparison, 
    get_kinase_family_distribution_network
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
# Preload structural similarity data at application startup
try:
    print("Preloading structural similarity data...")
    from protein_explorer.analysis.phospho_analyzer import preload_structural_data
    
    # Check for feather file first, fall back to parquet
    feather_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 
                             'Combined_Kinome_10A_Master_Filtered_2.feather')
    if os.path.exists(feather_file):
        preload_structural_data(feather_file)
    else:
        parquet_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 
                                 'Combined_Kinome_10A_Master_Filtered_2.parquet')
        preload_structural_data(parquet_file)
    
    print("Structural similarity data preloaded successfully")

    print("Preloading sequence similarity data...")
    seq_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 
                         'Sequence_Similarity_Edges.parquet')
    preload_sequence_data(seq_file)
    print("Sequence similarity data preloaded successfully")
except Exception as e:
    print(f"Warning: Failed to preload structural similarity data: {e}")
    print("Data will be loaded on first request (may cause delay)")

"""
Add this code near the top of app.py, just after the imports and before the Flask app initialization
"""

# Ensure cache directory exists
def ensure_cache_directory():
    """Ensure the cache directory exists at application startup."""
    # Paths for cache
    cache_dir = os.path.expanduser("~/.protein_explorer/cache")
    
    # Create cache directory if it doesn't exist
    if not os.path.exists(cache_dir):
        try:
            print(f"Creating cache directory: {cache_dir}")
            os.makedirs(cache_dir, exist_ok=True)
            print(f"Cache directory created successfully: {cache_dir}")
        except Exception as e:
            # If we can't create in home directory, use a temporary directory
            import tempfile
            alt_cache_dir = os.path.join(tempfile.gettempdir(), "protein_explorer_cache")
            print(f"Failed to create cache in home directory: {e}")
            print(f"Using alternative cache directory: {alt_cache_dir}")
            
            try:
                os.makedirs(alt_cache_dir, exist_ok=True)
            except Exception as e3:
                print(f"Failed to create alternative cache directory: {e3}")
                print("Application may have issues with caching")

# Run cache directory initialization
ensure_cache_directory()


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)

@app.route('/')
def index():
    """Render the home page."""
    return render_template('index.html')

@app.route('/search', methods=['GET', 'POST'])
def search():
    """Search for a protein by UniProt ID or gene symbol."""
    if request.method == 'POST':
        # Get search parameters
        identifier = request.form.get('identifier', '')
        id_type = request.form.get('id_type', 'uniprot')
        
        if not identifier:
            return render_template('search.html', error="Please enter an identifier")
        
        try:
            # Redirect to protein page
            return redirect(url_for('protein', identifier=identifier, id_type=id_type))
        except Exception as e:
            logger.error(f"Error in search: {e}")
            return render_template('search.html', error=str(e))
    
    # GET request
    return render_template('search.html')

@app.route('/protein/<identifier>')
def protein(identifier):
    """Display protein information and structure with full supplementary structural data."""
    try:
        print(f"DEBUG: Processing protein view for identifier: {identifier}")
        id_type = request.args.get('id_type', 'uniprot')
        
        # Get protein data
        if id_type.lower() == 'uniprot':
            print(f"DEBUG: Getting protein by UniProt ID")
            protein_data = pe.data.get_protein_by_id(uniprot_id=identifier)
        else:
            print(f"DEBUG: Getting protein by gene symbol")
            protein_data = pe.data.get_protein_by_id(gene_symbol=identifier)
        
        print(f"DEBUG: Protein data retrieved, has_structure: {protein_data.get('has_structure', False)}")
        
        # TEMPORARY FIX: For specific proteins, manually set has_structure=True
        if identifier in ['P04637', 'P53_HUMAN']:
            protein_data['has_structure'] = True
            print(f"DEBUG: Manually set has_structure=True for {identifier}")
        
        # Create cache directory if it doesn't exist
        import os
        cache_dir = os.path.expanduser("~/.protein_explorer/cache")
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir, exist_ok=True)
            print(f"DEBUG: Created cache directory: {cache_dir}")
        
        # Check if protein has a structure
        structure_html = None
        structure = None
        phosphosite_html = None
        phosphosites = None
        
        if protein_data.get('has_structure', False):
            print(f"DEBUG: Protein has structure, retrieving...")
            
            # Get structure 
            try:
                # Try direct download first for known proteins
                if identifier in ['P04637', 'P53_HUMAN']:
                    import requests
                    url = f"https://alphafold.ebi.ac.uk/files/AF-{protein_data['uniprot_id']}-F1-model_v4.pdb"
                    print(f"DEBUG: Trying direct URL: {url}")
                    response = requests.get(url)
                    if response.status_code == 200:
                        structure = response.text
                        print(f"DEBUG: Got structure directly, length: {len(structure)}")
                    else:
                        # Fall back to normal method
                        structure = pe.data.get_alphafold_structure(protein_data['uniprot_id'])
                else:
                    # Use normal method
                    structure = pe.data.get_alphafold_structure(protein_data['uniprot_id'])
                
                # Get the protein sequence from metadata if available
                sequence = None
                if protein_data.get('metadata', {}).get('sequence', {}).get('value'):
                    sequence = protein_data['metadata']['sequence']['value']
                
                # Create structure visualization
                if structure:
                    print(f"DEBUG: Creating structure visualization")
                    structure_html = pe.visualization.visualize_structure(
                        structure,
                        sequence=sequence
                    )
                    
                    # Analyze phosphorylation sites
                    if sequence:
                        print(f"DEBUG: Analyzing phosphorylation sites")
                        try:
                            # Use enhanced function that includes supplementary data
                            from protein_explorer.analysis.phospho_analyzer import get_phosphosites, get_phosphosite_data
                            phosphosites = get_phosphosites(protein_data['uniprot_id'])
                            print(f"DEBUG: Found {len(phosphosites)} potential phosphorylation sites")
                            
                            # Enhance phosphosites with additional metrics for visualization
                            for site in phosphosites:
                                # Check for any additional supplementary data not already added
                                if 'resno' in site:
                                    site_id = f"{protein_data['uniprot_id']}_{site['resno']}"
                                    supp_data = get_phosphosite_data(site_id)
                                    
                                    if supp_data:
                                        # Add any supplementary fields not already present
                                        for key in ['surface_accessibility', 'site_plddt', 
                                                'polar_aa_percent', 'nonpolar_aa_percent', 
                                                'acidic_aa_percent', 'basic_aa_percent']:
                                            if key in supp_data and supp_data[key] is not None and key not in site:
                                                site[key] = supp_data[key]
                            
                            # Generate enhanced phosphosite HTML with the enhanced_table function
                            from protein_explorer.analysis.enhanced_table import enhance_phosphosite_table
                            phosphosite_html = enhance_phosphosite_table(phosphosites, protein_data['uniprot_id'])
                            
                        except Exception as e:
                            print(f"DEBUG: Error analyzing phosphosites: {e}")
                            import traceback
                            print(traceback.format_exc())
                            
                            # Fall back to basic analysis without supplementary data
                            phosphosites = pe.analysis.phospho.analyze_phosphosites(
                                sequence, structure, uniprot_id=protein_data.get('uniprot_id')
                            )
                            
                            # Use the enhanced table function even with basic data
                            try:
                                from protein_explorer.analysis.enhanced_table import enhance_phosphosite_table
                                phosphosite_html = enhance_phosphosite_table(phosphosites, protein_data['uniprot_id'])
                            except Exception as e2:
                                print(f"DEBUG: Error using enhanced table: {e2}")
                                # Fall back to original HTML generation if enhanced_table fails
                                
                                # Generate phosphosite HTML with basic data (original version)
                                phosphosite_html = f"""
                                <div class="card mt-4">
                                    <div class="card-header">
                                        <h5 class="mb-0">Phosphorylation Site Analysis</h5>
                                    </div>
                                    <div class="card-body p-0">
                                        <div class="table-responsive">
                                            <table class="table table-striped table-hover">
                                                <thead class="thead-light">
                                                    <tr>
                                                        <th>Site</th>
                                                        <th>Motif (-7 to +7)</th>
                                                        <th>Mean pLDDT</th>
                                                        <th>Nearby Residues (10Å)</th>
                                                        <th>Known in PhosphositePlus</th>
                                                    </tr>
                                                </thead>
                                                <tbody id="phosphosite-table">
                                """
                                
                                # Add rows for each phosphosite with basic data
                                for site in phosphosites:
                                    phosphosite_html += f"""
                                    <tr>
                                        <td><a href="/site/{protein_data['uniprot_id']}/{site['site']}" class="site-link" data-resno="{site['resno']}">{site['site']}</a></td>
                                        <td><code class="motif-sequence">{site['motif']}</code></td>
                                        <td>{site['mean_plddt']}</td>
                                        <td>{site['nearby_count']}</td>
                                        <td>{"Yes" if site.get('is_known', False) else "No"}</td>
                                    </tr>
                                    """
                                
                                # Close the table (original version)
                                phosphosite_html += """
                                                </tbody>
                                            </table>
                                        </div>
                                    </div>
                                </div>
                                
                                <script>
                                    document.addEventListener('DOMContentLoaded', function() {
                                        // Add click handlers to site links
                                        const siteLinks = document.querySelectorAll('.site-link');
                                        siteLinks.forEach(link => {
                                            link.addEventListener('click', function(e) {
                                                e.preventDefault();
                                                const resno = this.getAttribute('data-resno');
                                                
                                                // Find the span in the sequence viewer
                                                const sequenceSpans = document.querySelectorAll('.sequence-viewer span');
                                                if (sequenceSpans.length > 0) {
                                                    // Find and click the span for this residue
                                                    const index = parseInt(resno) - 1;
                                                    if (index >= 0 && index < sequenceSpans.length) {
                                                        sequenceSpans[index].click();
                                                    }
                                                }
                                            });
                                        });
                                    });
                                </script>
                                """
            except Exception as e:
                print(f"DEBUG: Error getting structure: {e}")
                structure_html = f'<div class="alert alert-danger">Error loading structure: {str(e)}</div>'
        
        print(f"DEBUG: Rendering template with phosphosite_html: {'Present' if phosphosite_html else 'None'}")
        print(f"DEBUG: phosphosite_html content length: {len(phosphosite_html) if phosphosite_html else 'None'}")
        
        return render_template(
            'protein.html',
            protein=protein_data,
            structure_html=structure_html,
            phosphosites=phosphosites,  # Pass the phosphosites data
            phosphosite_html=phosphosite_html  # Pass the HTML for rendering
        )
    except Exception as e:
        print(f"DEBUG: Exception in protein view: {e}")
        logger.error(f"Error in protein view: {e}")
        return render_template('error.html', error=str(e))


@app.route('/analyze', methods=['GET', 'POST'])
def analyze():
    """Analyze multiple proteins."""
    if request.method == 'POST':
        # Get proteins from form
        proteins_text = request.form.get('proteins', '')
        analysis_type = request.form.get('analysis_type', 'network')
        
        # Parse protein list
        protein_list = [p.strip() for p in proteins_text.split(',') if p.strip()]
        
        if not protein_list:
            return render_template('analyze.html', error="Please enter at least one protein")
        
        try:
            # Determine ID type (assume all are the same type)
            id_type = 'gene' if '_' not in protein_list[0] else 'uniprot'
            
            # Convert to UniProt IDs if needed
            if id_type == 'gene':
                uniprot_ids = []
                for gene in protein_list:
                    try:
                        protein = pe.data.get_protein_by_id(gene_symbol=gene)
                        uniprot_ids.append(protein['uniprot_id'])
                    except Exception as e:
                        logger.warning(f"Could not convert gene {gene} to UniProt ID: {e}")
                protein_list = uniprot_ids
            
            # Perform the analysis
            if analysis_type == 'network':
                # Build network
                network = pe.navigation.build_interaction_network(protein_list, max_depth=1)
                
                # Create network visualization
                network_html = pe.visualization.visualize_network(
                    network,
                    highlight_nodes=protein_list,
                    title=f"Protein Interaction Network"
                )
                
                # Find common interactors
                common_interactors = pe.navigation.find_common_interactors(network, protein_list)
                
                return render_template(
                    'analyze_results.html',
                    analysis_type=analysis_type,
                    proteins=protein_list,
                    network_html=network_html,
                    common_interactors=common_interactors,
                    node_count=network.number_of_nodes(),
                    edge_count=network.number_of_edges()
                )
            elif analysis_type == 'structure':
                # Get structures and compare
                structures = {}
                for uniprot_id in protein_list:
                    try:
                        structure = pe.data.get_alphafold_structure(uniprot_id)
                        if structure:
                            structures[uniprot_id] = structure
                    except Exception as e:
                        logger.warning(f"Could not get structure for {uniprot_id}: {e}")
                
                if not structures:
                    return render_template('analyze.html', 
                                          error="Could not find structures for any of the proteins")
                
                # Create structure visualization
                structure_html = None
                if len(structures) == 1:
                    # Single structure
                    uniprot_id, structure = next(iter(structures.items()))
                    structure_html = pe.visualization.visualize_structure(structure)
                else:
                    # Compare structures
                    structure_html = pe.visualization.compare_structures(list(structures.values()))
                
                return render_template(
                    'analyze_results.html',
                    analysis_type=analysis_type,
                    proteins=protein_list,
                    structure_html=structure_html,
                    structures_found=list(structures.keys())
                )
            else:
                return render_template('analyze.html', error=f"Unknown analysis type: {analysis_type}")
                
        except Exception as e:
            logger.error(f"Error in analysis: {e}")
            return render_template('analyze.html', error=str(e))
    
    # GET request
    return render_template('analyze.html')

@app.route('/api/phosphosites/<uniprot_id>', methods=['GET'])
def api_phosphosites(uniprot_id):
    """API endpoint for phosphorylation site analysis with supplementary data."""
    try:
        # Import our enhanced phosphosite functions
        from protein_explorer.analysis.phospho_analyzer import get_phosphosites, get_phosphosite_data

        # Get protein data
        protein_data = pe.data.get_protein_by_id(uniprot_id=uniprot_id)
        
        # Get sequence
        sequence = protein_data.get('metadata', {}).get('sequence', {}).get('value')
        if not sequence:
            return jsonify({'error': 'Protein sequence not found'}), 404
            
        # Get structure
        structure = pe.data.get_alphafold_structure(uniprot_id)
        if not structure:
            return jsonify({'error': 'Protein structure not found'}), 404
            
        # Get phosphosites with enhanced data
        try:
            # Use the enhanced function first
            phosphosites = get_phosphosites(uniprot_id)
            
            # Add any additional supplementary data not already included
            for site in phosphosites:
                if 'resno' in site:
                    site_id = f"{uniprot_id}_{site['resno']}"
                    supp_data = get_phosphosite_data(site_id)
                    
                    if supp_data:
                        for key, value in supp_data.items():
                            # Only add fields not already present
                            if key not in site and value is not None:
                                site[key] = value
        except Exception as e:
            # Fall back to basic analysis without supplementary data
            logger.warning(f"Enhanced phosphosite analysis failed, using basic analysis: {e}")
            phosphosites = pe.analysis.phospho.analyze_phosphosites(sequence, structure, uniprot_id)
        
        # Get structural matches for each site
        try:
            from protein_explorer.analysis.phospho_analyzer import find_structural_matches, enhance_structural_matches
            
            # For each site, find and enhance matches
            for site in phosphosites:
                site_matches = find_structural_matches(uniprot_id, [site], top_n=10)
                if site['site'] in site_matches:
                    raw_matches = site_matches[site['site']]
                    # Enhance matches with supplementary data
                    enhanced_matches = enhance_structural_matches(raw_matches, site['site'])
                    # Add to site
                    site['structural_matches'] = enhanced_matches
        except Exception as e:
            logger.warning(f"Error finding structural matches: {e}")
            # Continue without matches
        
        return jsonify(phosphosites)
    except Exception as e:
        logger.error(f"Error in API endpoint: {e}")
        return jsonify({'error': str(e)}), 400

@app.route('/api/structure/<uniprot_id>', methods=['GET'])
def api_structure(uniprot_id):
    """API endpoint for protein structure."""
    try:
        structure = pe.data.get_alphafold_structure(uniprot_id)
        if not structure:
            return jsonify({'error': 'Structure not found'}), 404
        
        # Return basic info about the structure
        structure_info = pe.data.parse_pdb_structure(structure)
        return jsonify(structure_info)
    except Exception as e:
        return jsonify({'error': str(e)}), 400

@app.route('/api/network/<uniprot_id>', methods=['GET'])
def api_network(uniprot_id):
    """API endpoint for protein network."""
    try:
        depth = int(request.args.get('depth', 1))
        network = pe.navigation.build_interaction_network([uniprot_id], max_depth=depth)
        
        # Convert network to dict for JSON serialization
        network_data = {
            'nodes': list(network.nodes()),
            'edges': list(network.edges()),
            'node_count': network.number_of_nodes(),
            'edge_count': network.number_of_edges()
        }
        
        return jsonify(network_data)
    except Exception as e:
        return jsonify({'error': str(e)}), 400

    
@app.route('/phosphosite', methods=['GET', 'POST'])
def phosphosite_analysis():
    """Phosphosite structural analysis page with improved error handling and supplementary data."""
    import os
    from flask import request, render_template
    
    # Import phospho_analyzer functions
    from protein_explorer.analysis.phospho_analyzer import (
        analyze_protein, get_phosphosite_data, enhance_phosphosite, enhance_structural_matches
    )
    
    # Import enhanced table function
    from protein_explorer.analysis.enhanced_table import enhance_phosphosite_table
    
    # Initialize variables
    results = None
    error = None
    
    if request.method == 'POST':
        # Get identifier from form
        identifier = request.form.get('identifier', '')
        id_type = request.form.get('id_type', 'uniprot')
        
        if not identifier:
            return render_template('phosphosite.html', error="Please enter an identifier")
        
        try:
            # Find the data file
            parquet_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 
                         'Combined_Kinome_10A_Master_Filtered_2.feather')

            # Run the full analysis with supplementary data
            results = analyze_protein(identifier, id_type, parquet_file)
            
            # Use enhanced table for the phosphosites
            if results and 'phosphosites' in results and results['phosphosites']:
                # Create enhanced HTML for the phosphosites
                phosphosites_html = enhance_phosphosite_table(
                    results['phosphosites'], 
                    results['protein_info']['uniprot_id']
                )
                
                # Add the HTML to the results
                results['phosphosites_html'] = phosphosites_html
            
            # Return results including enhanced data
            return render_template('phosphosite.html', 
                                  protein_info=results['protein_info'],
                                  phosphosites=results['phosphosites'],
                                  phosphosites_html=results.get('phosphosites_html'),
                                  structural_matches=results['structural_matches'],
                                  error=results.get('error'))
                
        except Exception as e:
            print(f"Error in phosphosite analysis: {e}")
            error = str(e)
            return render_template('phosphosite.html', error=error)
    
    # GET request - show empty form
    return render_template('phosphosite.html')

# Add this to app.py
@app.route('/api/sequence_matches/<site_id>', methods=['GET'])
def api_sequence_matches(site_id):
    """API endpoint for sequence similarity matches."""
    try:
        # Get query parameters
        top_n = request.args.get('top_n', default=100, type=int)
        min_similarity = request.args.get('min_similarity', default=0.4, type=float)
        
        # Validate site_id format (e.g., UniProtID_ResidueNumber)
        if '_' not in site_id:
            return jsonify({'error': 'Invalid site ID format. Expected: UniProtID_ResidueNumber'}), 400
            
        # Get sequence matches
        matches = find_sequence_matches(site_id, top_n=top_n, min_similarity=min_similarity)
        
        # If no matches, return empty list with a message
        if not matches:
            return jsonify({
                'matches': [],
                'count': 0,
                'message': f'No sequence similarity matches found for {site_id}'
            })
        
        # Create a complete response with metadata
        response = {
            'site_id': site_id,
            'matches': matches,
            'count': len(matches),
            'params': {
                'top_n': top_n,
                'min_similarity': min_similarity
            }
        }
        
        return jsonify(response)
    except Exception as e:
        logger.error(f"Error in sequence matches API: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/sequence_conservation/<site_id>', methods=['GET'])
def api_sequence_conservation(site_id):
    """API endpoint for sequence conservation analysis."""
    try:
        # Get query parameters
        top_n = request.args.get('top_n', default=200, type=int)
        min_similarity = request.args.get('min_similarity', default=0.4, type=float)
        
        # Get the query motif if possible
        query_motif = None
        
        # Parse site_id to extract uniprot_id and site_number
        parts = site_id.split('_')
        if len(parts) >= 2:
            uniprot_id = parts[0]
            site_number = int(parts[1])
            
            # Try to get site data with motif
            try:
                # Get all sites for the protein
                from protein_explorer.analysis.phospho_analyzer import get_phosphosites
                sites = get_phosphosites(uniprot_id)
                
                # Find the specific site
                for site in sites:
                    if site['resno'] == site_number:
                        query_motif = site.get('motif')
                        break
            except:
                # If that fails, try direct motif lookup
                pass
        
        # Get sequence matches
        matches = find_sequence_matches(site_id, top_n=top_n, min_similarity=min_similarity)
        
        # Analyze conservation
        conservation = analyze_motif_conservation(matches, query_motif=query_motif)
        
        # Add metadata to response
        response = {
            'site_id': site_id,
            'analysis': conservation,
            'match_count': len(matches),
            'query_motif': query_motif
        }
        
        return jsonify(response)
    except Exception as e:
        logger.error(f"Error in sequence conservation API: {e}")
        return jsonify({'error': str(e)}), 500


@app.route('/site/<uniprot_id>/<site>')
def site_detail(uniprot_id, site):
    """Display detailed information about a specific phosphorylation site with enhanced supplementary data."""
    try:
        # Get protein data
        protein_data = pe.data.get_protein_by_id(uniprot_id=uniprot_id)
        
        # Get the sequence from metadata
        sequence = protein_data.get('metadata', {}).get('sequence', {}).get('value')
        if not sequence:
            return render_template('error.html', error=f"Protein sequence not found for {uniprot_id}")
            
        # Get structure
        structure = pe.data.get_alphafold_structure(uniprot_id)
        if not structure:
            return render_template('error.html', error=f"Protein structure not found for {uniprot_id}")
        
        # Parse the site string to get type and residue number
        site_match = re.match(r'([A-Z])(\d+)', site)
        if not site_match:
            return render_template('error.html', error=f"Invalid site format: {site}")
            
        site_type = site_match.group(1)
        site_number = int(site_match.group(2))
        
        # Generate site_id in the format UniProtID_ResidueNumber
        site_id = f"{uniprot_id}_{site_number}"
        
        # Get supplementary data for this site
        from protein_explorer.analysis.phospho_analyzer import get_phosphosite_data, enhance_phosphosite
        site_data = get_phosphosite_data(site_id)
        
        # If site_data is None, create a minimal site data dictionary
        if not site_data:
            site_data = {
                'site_id': site_id,
                'site': site,
                'resno': site_number,
                'siteType': site_type
            }
        
        # Analyze phosphosites
        from protein_explorer.analysis.phospho import analyze_phosphosites
        all_phosphosites = analyze_phosphosites(sequence, structure, uniprot_id)
        
        # Find the specific site
        site_data_from_analysis = next((s for s in all_phosphosites if s['site'] == site), None)
        if site_data_from_analysis:
            # Merge with existing site_data
            for key, value in site_data_from_analysis.items():
                if key not in site_data or site_data[key] is None:
                    site_data[key] = value
        
        # Find structural matches
        from protein_explorer.analysis.phospho_analyzer import find_structural_matches, enhance_structural_matches
        structural_matches = []
        try:
            # Get raw matches
            all_matches = find_structural_matches(uniprot_id, [site_data], top_n=None)
            raw_matches = all_matches.get(site, [])
            
            # Filter out self-matches
            raw_matches = [match for match in raw_matches if match.get('rmsd', 0) > 0.01]
            
            # Enhance matches with supplementary data
            structural_matches = enhance_structural_matches(raw_matches, site)
            
            # Sort matches by RMSD
            structural_matches = sorted(structural_matches, key=lambda x: x.get('rmsd', float('inf')))
            
        except Exception as e:
            logger.error(f"Error finding structural matches: {e}")
            structural_matches = []
        
        # Create basic structure visualization
        structure_html = pe.visualization.visualize_structure(
            structure,
            sequence=sequence
        )
        
        # Create specialized visualizations
        from protein_explorer.visualization.network import create_phosphosite_network
        from protein_explorer.analysis.phospho_analyzer import (
            create_comparative_motif_visualization,
            enhance_site_visualization,
            analyze_residue_distributions,
            analyze_phosphosite_context
        )
        
        # 1. Network visualization
        try:
            network_html = create_phosphosite_network(site, structural_matches, site_data)
        except Exception as e:
            logger.error(f"Error creating network visualization: {e}")
            network_html = f"<div class='alert alert-warning'>Error creating network visualization: {str(e)}</div>"
            
        # 2. Comparative motif visualization
        try:
            if structural_matches and 'motif' in site_data:
                motif_html = create_comparative_motif_visualization(site_data, structural_matches)
            else:
                motif_html = "<div class='alert alert-info'>No motif data available for comparison.</div>"
        except Exception as e:
            logger.error(f"Error creating motif visualization: {e}")
            motif_html = f"<div class='alert alert-warning'>Error creating motif visualization: {str(e)}</div>"
            
        # 3. Residue distribution analysis
        try:
            if structural_matches:
                distribution_data = analyze_residue_distributions(structural_matches)
            else:
                distribution_data = None
        except Exception as e:
            logger.error(f"Error analyzing residue distributions: {e}")
            distribution_data = None
            
        # 4. Enhanced 3D visualization
        try:
            enhanced_3d_html = enhance_site_visualization(uniprot_id, site, site_data)
        except Exception as e:
            logger.error(f"Error creating enhanced 3D visualization: {e}")
            enhanced_3d_html = f"<div class='alert alert-warning'>Error creating enhanced 3D visualization: {str(e)}</div>"
            
        # 5. Structural context analysis
        try:
            context_data = analyze_phosphosite_context(structure, site_number, site_type)
        except Exception as e:
            logger.error(f"Error analyzing structural context: {e}")
            context_data = None
        
        # Find sequence similarity matches
        sequence_matches = []
        sequence_network_data = None
        sequence_conservation = None
        sequence_motif_html = None
        
        try:
            # Get sequence similarity matches
            logger.info(f"Finding sequence matches for {site_id}")
            
            from protein_explorer.analysis.sequence_analyzer import (
                find_sequence_matches,
                analyze_motif_conservation,
                create_sequence_network_data,
                create_sequence_motif_visualization
            )
            
            sequence_matches = find_sequence_matches(site_id, top_n=200, min_similarity=0.4)
            logger.info(f"Found {len(sequence_matches)} sequence matches")
            
            # Log a sample of the matches to see what they contain
            if sequence_matches:
                logger.info(f"Sample match: {sequence_matches[0]}")
                
                # Create sequence network data
                sequence_network_data = create_sequence_network_data(
                    site_id, 
                    sequence_matches,
                    query_motif=site_data.get('motif')
                )
                logger.info(f"Created network data with {len(sequence_network_data['nodes'])} nodes")
                
                # Create sequence conservation analysis
                sequence_conservation = analyze_motif_conservation(
                    sequence_matches,
                    query_motif=site_data.get('motif')
                )
                logger.info(f"Created conservation analysis with {sequence_conservation['motif_count']} motifs")
                
                # Create sequence motif visualization
                sequence_motif_html = create_sequence_motif_visualization(
                    site_id,
                    site_data.get('motif', ''),
                    sequence_matches
                )
                logger.info(f"Created motif visualization HTML of length {len(sequence_motif_html) if sequence_motif_html else 0}")
            else:
                logger.warning(f"No sequence matches found for {site_id}")
        except Exception as e:
            logger.error(f"Error analyzing sequence similarity: {e}")
            import traceback
            logger.error(traceback.format_exc())
            sequence_matches = []
        
        # Standard kinase prediction data
        structure_kinase_data = {}
        sequence_kinase_data = {}
        
        # Network-based kinase prediction data
        structure_network_kinase_data = {}
        sequence_network_kinase_data = {}
        
        try:
            # Import kinase prediction functions
            from protein_explorer.analysis.kinase_predictor import (
                get_site_kinase_scores,
                predict_kinases,
                get_heatmap_data,
                get_known_kinase_info,
                categorize_kinases_by_family
            )
            
            # Import network-based kinase prediction functions
            from protein_explorer.analysis.network_kinase_predictor import (
                get_similar_sites, 
                compute_aggregated_kinase_scores, 
                predict_kinases_network,
                get_network_heatmap_data, 
                get_network_kinase_comparison, 
                get_kinase_family_distribution_network
            )
            
            # Standard structural kinase prediction
            struct_scores = get_site_kinase_scores(site_id, 'structure')
            if struct_scores and 'scores' in struct_scores:
                # Get top kinases
                top_kinases = predict_kinases(site_id, top_n=10, score_type='structure')
                
                # Get kinase families
                kinase_families = categorize_kinases_by_family(top_kinases)
                
                # Get known kinase info
                known_kinase = get_known_kinase_info(site_id, 'structure')
                
                # Prepare data for heatmap visualization
                struct_match_ids = [f"{match['target_uniprot']}_{match['target_site'].replace('S', '').replace('T', '').replace('Y', '')}" 
                                   for match in structural_matches[:20] if match.get('rmsd', 10) < 5.0]
                struct_match_ids.insert(0, site_id)  # Add the query site first
                
                heatmap_data = get_heatmap_data(struct_match_ids, top_n=10, score_type='structure')
                
                # Bundle all data
                structure_kinase_data = {
                    'site_id': site_id,
                    'top_kinases': top_kinases,
                    'kinase_families': kinase_families,
                    'known_kinase': known_kinase,
                    'heatmap': heatmap_data
                }
            
            # Standard sequence kinase prediction
            seq_scores = get_site_kinase_scores(site_id, 'sequence')
            if seq_scores and 'scores' in seq_scores:
                # Get top kinases
                top_kinases = predict_kinases(site_id, top_n=10, score_type='sequence')
                
                # Get kinase families
                kinase_families = categorize_kinases_by_family(top_kinases)
                
                # Get known kinase info
                known_kinase = get_known_kinase_info(site_id, 'sequence')
                
                # Prepare data for heatmap visualization
                seq_match_ids = []
                if sequence_matches:
                    seq_match_ids = [match['target_id'] for match in sequence_matches[:20] 
                                    if match.get('similarity', 0) > 0.6]
                    seq_match_ids.insert(0, site_id)  # Add the query site first
                
                heatmap_data = get_heatmap_data(seq_match_ids, top_n=10, score_type='sequence')
                
                # Bundle all data
                sequence_kinase_data = {
                    'site_id': site_id,
                    'top_kinases': top_kinases,
                    'kinase_families': kinase_families,
                    'known_kinase': known_kinase,
                    'heatmap': heatmap_data
                }
            
            # Network-based kinase predictions
            try:
                # 1. Structure network prediction
                structure_network_kinases = predict_kinases_network(
                    site_id, top_n=10, score_type='structure',
                    similarity_threshold=0.6, rmsd_threshold=3.0
                )
                
                if structure_network_kinases:
                    # Get similar sites
                    struct_similar_sites = get_similar_sites(
                        site_id, uniprot_id, site_data,
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Get heatmap data with enhanced function
                    struct_heatmap_data = get_network_heatmap_data(
                        site_id, top_n=10, score_type='structure',
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Log heatmap data for debugging
                    logger.info(f"Structure heatmap data: {len(struct_heatmap_data['sites'])} sites, " 
                                f"{len(struct_heatmap_data['kinases'])} kinases, "
                                f"{len(struct_heatmap_data['scores'])} scores")
                    
                    # Get family distribution
                    struct_family_dist = get_kinase_family_distribution_network(
                        site_id, score_type='structure',
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Combine into data structure
                    structure_network_kinase_data = {
                        'site_id': site_id,
                        'top_kinases': structure_network_kinases,
                        'heatmap': struct_heatmap_data,
                        'kinase_families': struct_family_dist,
                        'site_count': len(struct_similar_sites),
                        'rmsd_threshold': 3.0,
                        'similarity_threshold': 0.6
                    }
                
                # 2. Sequence network prediction
                sequence_network_kinases = predict_kinases_network(
                    site_id, top_n=10, score_type='sequence',
                    similarity_threshold=0.6, rmsd_threshold=3.0
                )
                
                if sequence_network_kinases:
                    # Get similar sites
                    seq_similar_sites = get_similar_sites(
                        site_id, uniprot_id, site_data,
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Get heatmap data with enhanced function
                    seq_heatmap_data = get_network_heatmap_data(
                        site_id, top_n=10, score_type='sequence',
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Log heatmap data for debugging
                    logger.info(f"Sequence heatmap data: {len(seq_heatmap_data['sites'])} sites, " 
                                f"{len(seq_heatmap_data['kinases'])} kinases, "
                                f"{len(seq_heatmap_data['scores'])} scores")
                    
                    # Get family distribution
                    seq_family_dist = get_kinase_family_distribution_network(
                        site_id, score_type='sequence',
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Combine into data structure
                    sequence_network_kinase_data = {
                        'site_id': site_id,
                        'top_kinases': sequence_network_kinases,
                        'heatmap': seq_heatmap_data,
                        'kinase_families': seq_family_dist,
                        'site_count': len(seq_similar_sites),
                        'rmsd_threshold': 3.0,
                        'similarity_threshold': 0.6
                    }
            except Exception as network_e:
                logger.error(f"Error in network kinase prediction: {network_e}")
                import traceback
                logger.error(traceback.format_exc())
                
                # Provide minimal fallback data to avoid template errors
                structure_network_kinase_data = {
                    'site_id': site_id,
                    'top_kinases': [],
                    'heatmap': {'sites': [], 'kinases': [], 'scores': []},
                    'kinase_families': {},
                    'site_count': 0,
                    'rmsd_threshold': 3.0,
                    'similarity_threshold': 0.6,
                    'error': str(network_e)
                }
                
                sequence_network_kinase_data = {
                    'site_id': site_id,
                    'top_kinases': [],
                    'heatmap': {'sites': [], 'kinases': [], 'scores': []},
                    'kinase_families': {},
                    'site_count': 0,
                    'rmsd_threshold': 3.0,
                    'similarity_threshold': 0.6,
                    'error': str(network_e)
                }

        except Exception as e:
            logger.error(f"Error processing kinase data: {e}")
            import traceback
            logger.error(traceback.format_exc())
        
        # Render the template with all the data
        return render_template(
            'site.html',
            protein=protein_data,
            site=site,
            site_data=site_data,
            site_id=site_id,  # Make sure site_id is passed to the template
            structure_html=structure_html,
            structural_matches=structural_matches,
            supplementary_data=site_data,
            network_html=network_html,
            motif_html=motif_html,
            distribution_data=distribution_data,
            enhanced_3d_html=enhanced_3d_html,
            context_data=context_data,
            # Sequence similarity data:
            sequence_matches=sequence_matches,
            sequence_network_data=sequence_network_data,
            sequence_conservation=sequence_conservation,
            sequence_motif_html=sequence_motif_html,
            # Standard kinase prediction data:
            structure_kinase_data=structure_kinase_data,
            sequence_kinase_data=sequence_kinase_data,
            # Network kinase prediction data:
            structure_network_kinase_data=structure_network_kinase_data,
            sequence_network_kinase_data=sequence_network_kinase_data
        )
    except Exception as e:
        logger.error(f"Error in site detail view: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return render_template('error.html', error=str(e))
    

def analyze_motif(motif: str, site_type: str, site_number: int) -> Dict:
    """
    Analyze a phosphosite motif sequence for additional insights.
    
    Args:
        motif: The motif sequence string
        site_type: The site type (S, T, or Y)
        site_number: The residue number
        
    Returns:
        Dictionary with motif analysis
    """
    # Find the position of the phosphosite in the motif
    motif_length = len(motif)
    
    # Identify the center residue position in the motif
    center_pos = None
    for i, aa in enumerate(motif):
        # Look for site_type (S, T, or Y) in appropriate position
        if aa == site_type:
            # Calculate expected position from start of motif
            expected_pos = site_number - 8  # If site is the center of a -7 to +7 motif
            if i == 7:  # Expected position for -7 to +7 motif (0-indexed)
                center_pos = i
                break
    
    if center_pos is None:
        # If we can't find it exactly, assume it's in the middle
        center_pos = motif_length // 2
    
    # Count amino acid types
    aa_groups = {
        'polar': 'STYCNQ',
        'nonpolar': 'AVILMFWPG',
        'acidic': 'DE',
        'basic': 'KRH',
        'other': 'X'
    }
    
    counts = {group: 0 for group in aa_groups}
    upstream_counts = {group: 0 for group in aa_groups}
    downstream_counts = {group: 0 for group in aa_groups}
    
    for i, aa in enumerate(motif):
        # Skip if it's the phosphosite
        if i == center_pos:
            continue
            
        # Determine group
        group = 'other'
        for g, aas in aa_groups.items():
            if aa in aas:
                group = g
                break
                
        # Increment total count
        counts[group] += 1
        
        # Increment upstream/downstream count
        if i < center_pos:
            upstream_counts[group] += 1
        else:
            downstream_counts[group] += 1
    
    # Calculate percentages
    total_aas = motif_length - 1  # Exclude phosphosite
    upstream_total = center_pos
    downstream_total = motif_length - center_pos - 1
    
    percentages = {group: count/total_aas*100 for group, count in counts.items()}
    upstream_percentages = {f"upstream_{group}": count/upstream_total*100 if upstream_total > 0 else 0 
                          for group, count in upstream_counts.items()}
    downstream_percentages = {f"downstream_{group}": count/downstream_total*100 if downstream_total > 0 else 0 
                            for group, count in downstream_counts.items()}
    
    # Extract -3 to +3 submotif (important region for kinase recognition)
    narrow_motif = motif[max(0, center_pos-3):min(motif_length, center_pos+4)]
    narrow_motif_highlighted = f"{narrow_motif[:3]}<strong>{narrow_motif[3:4]}</strong>{narrow_motif[4:]}" if len(narrow_motif) >= 7 else narrow_motif
    
    # Create motif visualization with color coding
    motif_vis = []
    for i, aa in enumerate(motif):
        position = i - center_pos
        # Determine AA type for coloring
        aa_class = "highlighted" if i == center_pos else "other"
        for group, aas in aa_groups.items():
            if aa in aas:
                aa_class = group
                break
                
        motif_vis.append({
            "aa": aa,
            "position": position,
            "class": aa_class
        })
    
    # Check for known kinase recognition motifs
    motif_patterns = [
        {"name": "CDK", "pattern": r"[ST]P.?[KR]", "kinases": ["CDK1", "CDK2", "CDK5"]},
        {"name": "MAPK/ERK", "pattern": r"P.[ST]P", "kinases": ["ERK1", "ERK2", "p38", "JNK"]},
        {"name": "PKA/PKG", "pattern": r"[RK][RK].?[ST]", "kinases": ["PKA", "PKG"]},
        {"name": "CK2", "pattern": r"[ST]..[DE]", "kinases": ["CK2"]},
        {"name": "GSK3", "pattern": r"[ST]..[ST]", "kinases": ["GSK3α", "GSK3β"]},
        {"name": "CK1", "pattern": r"[DE]..[ST].[ST]", "kinases": ["CK1"]},
        {"name": "PLK", "pattern": r"[DE]..[ST]", "kinases": ["PLK1"]},
        {"name": "Akt/PKB", "pattern": r"R.R..[ST]", "kinases": ["Akt1", "Akt2", "Akt3"]},
        {"name": "AGC", "pattern": r"R.[ST]", "kinases": ["PKA", "PKC", "PKG"]}
    ]
    
    potential_kinases = []
    for pattern_info in motif_patterns:
        pattern = pattern_info["pattern"]
        match = re.search(pattern, motif)
        if match:
            potential_kinases.append({
                "pattern_name": pattern_info["name"],
                "pattern": pattern,
                "kinases": pattern_info["kinases"],
                "match": match.group(0)
            })
    
    return {
        "motif_vis": motif_vis,
        "center_position": center_pos,
        "narrow_motif": narrow_motif,
        "narrow_motif_highlighted": narrow_motif_highlighted,
        "aa_counts": counts,
        "aa_percentages": percentages,
        "upstream_percentages": upstream_percentages,
        "downstream_percentages": downstream_percentages,
        "potential_kinases": potential_kinases
    }


def create_site_focused_model(structure_data, residue_number, residue_type):
    """
    Create a 3D model visualization focused on a specific residue.
    Highlights the residue and its surrounding environment.
    Enhanced version with better color coding and structural context.
    """
    # Use the existing visualize_structure function but add a highlight for the specific residue
    import base64
    
    # Base64 encode the PDB data
    pdb_base64 = base64.b64encode(structure_data.encode()).decode()
    
    # Create a customized version of the NGL viewer script that focuses on the site
    js_code = f"""
    <style>
        .site-viewer {{
            width: 100%;
            height: 400px;
            position: relative;
        }}
        .site-info-panel {{
            position: absolute;
            top: 10px;
            right: 10px;
            background-color: rgba(255, 255, 255, 0.9);
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 8px;
            font-size: 14px;
            z-index: 100;
        }}
        .aa-legend {{
            position: absolute;
            bottom: 10px;
            right: 10px;
            background-color: rgba(255, 255, 255, 0.9);
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 8px;
            font-size: 12px;
            z-index: 100;
        }}
        .aa-group {{
            display: flex;
            align-items: center;
            margin-bottom: 3px;
        }}
        .aa-color {{
            width: 12px;
            height: 12px;
            border-radius: 50%;
            margin-right: 5px;
        }}
    </style>
    
    <div class="site-viewer">
        <div id="site-viewer-container" style="width: 100%; height: 100%;"></div>
        <div class="site-info-panel">
            <strong>Site:</strong> {residue_type}{residue_number}<br>
            <strong>View:</strong> <span id="view-mode" style="cursor:pointer;">Site Environment</span><br>
            <strong>Mode:</strong> <span id="color-mode" style="cursor:pointer;">Element</span>
        </div>
        <div class="aa-legend">
            <div class="aa-group"><div class="aa-color" style="background-color:#FF4500;"></div> Target Site</div>
            <div class="aa-group"><div class="aa-color" style="background-color:#87CEFA;"></div> Polar</div>
            <div class="aa-group"><div class="aa-color" style="background-color:#FFD700;"></div> Non-polar</div>
            <div class="aa-group"><div class="aa-color" style="background-color:#FF6347;"></div> Acidic</div>
            <div class="aa-group"><div class="aa-color" style="background-color:#98FB98;"></div> Basic</div>
        </div>
    </div>
    
    <script src="https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.37/dist/ngl.js"></script>
    <script type="text/javascript">
        document.addEventListener('DOMContentLoaded', function() {{
            // Create NGL Stage object
            var stage = new NGL.Stage('site-viewer-container', {{backgroundColor: "white"}});
            
            // Handle window resizing
            window.addEventListener('resize', function() {{
                stage.handleResize();
            }});
            
            // Define residue type groupings
            const aminoAcidGroups = {{
                polar: ["SER", "THR", "TYR", "CYS", "ASN", "GLN"],
                nonpolar: ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "GLY"],
                acidic: ["ASP", "GLU"],
                basic: ["LYS", "ARG", "HIS"]
            }};
            
            // Define colors for each group
            const groupColors = {{
                polar: [135/255, 206/255, 250/255],     // Light blue
                nonpolar: [255/255, 215/255, 0/255],    // Gold
                acidic: [255/255, 99/255, 71/255],      // Tomato
                basic: [152/255, 251/255, 152/255]      // Pale green
            }};
            
            // UI elements
            const viewModeText = document.getElementById('view-mode');
            const colorModeText = document.getElementById('color-mode');
            
            // State variables
            let isFullView = false;
            let colorMode = "element";  // "element", "group", "plddt"
            
            // Function to determine AA group from residue name
            function getAminoAcidGroup(resname) {{
                for (const [group, residues] of Object.entries(aminoAcidGroups)) {{
                    if (residues.includes(resname)) {{
                        return group;
                    }}
                }}
                return "other";
            }}
            
            // Function to get color for residue based on its group
            function colorByType(atom) {{
                const residue = atom.resname;
                
                // Special case for target residue
                if (atom.resno === {residue_number} && atom.resname.includes("{residue_type}")) {{
                    return [1.0, 0.27, 0.0];  // #FF4500 orange-red
                }}
                
                // Determine group
                const group = getAminoAcidGroup(residue);
                if (group in groupColors) {{
                    return groupColors[group];
                }}
                
                // Default color for "other"
                return [0.5, 0.5, 0.5];  // Grey
            }}
            
            // Load PDB data from base64 string
            var pdbBlob = new Blob([atob('{pdb_base64}')], {{type: 'text/plain'}});
            
            // Load the structure
            stage.loadFile(pdbBlob, {{ext: 'pdb'}}).then(function(component) {{
                // Get target selection
                const siteSelection = "{residue_number} and .{residue_type}";
                const environmentSelection = siteSelection + " or (" + siteSelection + " around 5)";
                
                // State variables
                let isFullView = false;
                let colorMode = "element";  // "element" or "type"
                
                // Reset button
                document.getElementById('reset-view-btn')?.addEventListener('click', function() {{
                    updateRepresentations();
                }});
                
                // Toggle view button
                viewModeText.addEventListener('click', function() {{
                    isFullView = !isFullView;
                    this.textContent = isFullView ? 'Site Focus' : 'Full Protein';
                    updateRepresentations();
                }});
                
                // Toggle color button
                colorModeText.addEventListener('click', function() {{
                    colorMode = colorMode === "element" ? "type" : "element";
                    this.textContent = colorMode === "element" ? 'Color by Type' : 'Color by Element';
                    updateRepresentations();
                }});
                
                // Update all representations based on current state
                function updateRepresentations() {{
                    // Remove all existing representations
                    component.removeAllRepresentations();
                    
                    // Add cartoon representation for entire protein
                    component.addRepresentation("cartoon", {{
                        color: colorMode === "type" ? colorByType : "chainid",
                        opacity: 0.7,
                        smoothSheet: true
                    }});
                    
                    // Add ball and stick for target residue
                    component.addRepresentation("ball+stick", {{
                        sele: siteSelection,
                        color: colorMode === "type" ? colorByType : "element",
                        aspectRatio: 1.5,
                        scale: 1.2
                    }});
                    
                    // Add licorice for environment (if not full view)
                    if (!isFullView) {{
                        component.addRepresentation("licorice", {{
                            sele: environmentSelection + " and not " + siteSelection,
                            color: colorMode === "type" ? colorByType : "element",
                            opacity: 0.8,
                            scale: 0.8
                        }});
                        
                        // Add labels
                        component.addRepresentation("label", {{
                            sele: environmentSelection,
                            color: "#333333",
                            labelType: "format",
                            labelFormat: "{{resname}}{{resno}}",
                            labelGrouping: "residue",
                            attachment: "middle-center",
                            showBackground: true,
                            backgroundColor: "white",
                            backgroundOpacity: 0.5
                        }});
                    }}
                    
                    // Set view
                    if (isFullView) {{
                        component.autoView();
                    }} else {{
                        component.autoView(environmentSelection, 2000);
                    }}
                }}
                
                // Initial setup
                updateRepresentations();
            }});
        }});
    </script>
    """
    
    return js_code

@app.route('/site-search', methods=['GET', 'POST'])
def site_search():
    """Search for a specific phosphorylation site."""
    if request.method == 'POST':
        # Get form data
        uniprot_id = request.form.get('uniprot_id', '')
        site = request.form.get('site', '')
        
        if not uniprot_id:
            return render_template('site-search.html', error="Please enter a UniProt ID")
        
        if not site:
            return render_template('site-search.html', error="Please enter a site identifier")
            
        # Validate site format (S, T, or Y followed by a number)
        import re
        if not re.match(r'^[STY]\d+$', site):
            return render_template('site-search.html', error="Invalid site format. Expected format: S123, T45, Y678, etc.")
            
        # Redirect to site detail page
        return redirect(url_for('site_detail', uniprot_id=uniprot_id, site=site))
    
    # GET request - show search form
    return render_template('site-search.html')

# Add a new API endpoint for kinase data
@app.route('/api/kinases/<site_id>', methods=['GET'])
def api_kinases(site_id):
    """API endpoint for kinase prediction scores."""
    try:
        score_type = request.args.get('type', 'structure')
        top_n = int(request.args.get('top_n', 10))
        
        # Get kinase scores
        scores = get_site_kinase_scores(site_id, score_type)
        if not scores or 'scores' not in scores:
            return jsonify({'error': f'No {score_type} kinase scores found for {site_id}'}), 404
        
        # Get top predicted kinases
        top_kinases = predict_kinases(site_id, top_n, score_type)
        
        # Get known kinase info
        known_kinase = get_known_kinase_info(site_id, score_type)
        
        # Get comparison with the other score type
        other_type = 'sequence' if score_type == 'structure' else 'structure'
        comparison = get_kinase_comparison_data(site_id, [score_type, other_type], top_n)
        
        # Return all data
        return jsonify({
            'site_id': site_id,
            'score_type': score_type,
            'known_kinase': known_kinase,
            'top_kinases': top_kinases,
            'comparison': comparison
        })
    except Exception as e:
        logger.error(f"Error in kinase API: {e}")
        return jsonify({'error': str(e)}), 500

# Add another API endpoint for comparing kinases between sites
@app.route('/api/kinases/compare', methods=['POST'])
def api_compare_kinases():
    """API endpoint for comparing kinase predictions across multiple sites."""
    try:
        data = request.get_json()
        if not data or 'site_ids' not in data:
            return jsonify({'error': 'Missing site_ids parameter'}), 400
        
        site_ids = data['site_ids']
        score_type = data.get('type', 'structure')
        top_n = int(data.get('top_n', 10))
        
        # Get heatmap data
        heatmap_data = get_heatmap_data(site_ids, top_n, score_type)
        
        # Return the data
        return jsonify({
            'heatmap': heatmap_data
        })
    except Exception as e:
        logger.error(f"Error in kinase comparison API: {e}")
        return jsonify({'error': str(e)}), 500
    

# Add a new API endpoint for network-based kinase predictions
@app.route('/api/network-kinases/<site_id>', methods=['GET'])
def api_network_kinases(site_id):
    """API endpoint for network-based kinase prediction scores."""
    try:
        score_type = request.args.get('type', 'structure')
        top_n = int(request.args.get('top_n', 10))
        similarity_threshold = float(request.args.get('similarity_threshold', 0.6))
        rmsd_threshold = float(request.args.get('rmsd_threshold', 3.0))
        
        # Get network-based kinase predictions
        network_kinases = predict_kinases_network(
            site_id, top_n=top_n, score_type=score_type,
            similarity_threshold=similarity_threshold, 
            rmsd_threshold=rmsd_threshold
        )
        
        if not network_kinases:
            return jsonify({'error': f'No {score_type} network kinase scores found for {site_id}'}), 404
        
        # Get similar sites
        similar_sites = get_similar_sites(
            site_id,
            similarity_threshold=similarity_threshold, 
            rmsd_threshold=rmsd_threshold
        )
        
        # Get heatmap data
        heatmap_data = get_network_heatmap_data(
            site_id, top_n=top_n, score_type=score_type,
            similarity_threshold=similarity_threshold, 
            rmsd_threshold=rmsd_threshold
        )
        
        # Get kinase families
        kinase_families = get_kinase_family_distribution_network(
            site_id, score_type=score_type,
            similarity_threshold=similarity_threshold, 
            rmsd_threshold=rmsd_threshold
        )
        
        # Return all data
        return jsonify({
            'site_id': site_id,
            'score_type': score_type,
            'top_kinases': network_kinases,
            'heatmap': heatmap_data,
            'kinase_families': kinase_families,
            'site_count': len(similar_sites),
            'rmsd_threshold': rmsd_threshold,
            'similarity_threshold': similarity_threshold
        })
    except Exception as e:
        logger.error(f"Error in network kinase API: {e}")
        return jsonify({'error': str(e)}), 500

# Add an endpoint to update network predictions with new thresholds
@app.route('/api/update-network-kinases/<site_id>', methods=['POST'])
def update_network_kinases(site_id):
    """Update network-based kinase predictions with new thresholds with enhanced error handling."""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'error': 'No data provided'}), 400
        
        score_type = data.get('score_type', 'structure')
        top_n = int(data.get('top_n', 10))
        
        # Get appropriate threshold based on score_type
        if score_type == 'structure':
            threshold_param = 'rmsd_threshold'
            threshold_default = 3.0
        else:  # sequence
            threshold_param = 'similarity_threshold'
            threshold_default = 0.6
            
        threshold = float(data.get(threshold_param, threshold_default))
        
        logger.info(f"Updating {score_type} network kinase predictions with {threshold_param}={threshold}")
        
        # Common parameters for all function calls
        kwargs = {
            'site_id': site_id,
            'top_n': top_n,
            'score_type': score_type
        }
        
        # Add the appropriate threshold parameter
        if score_type == 'structure':
            kwargs['rmsd_threshold'] = threshold
            # Use default for similarity
            kwargs['similarity_threshold'] = 0.6
        else:
            kwargs['similarity_threshold'] = threshold
            # Use default for RMSD
            kwargs['rmsd_threshold'] = 3.0
        
        # Get updated predictions
        network_kinases = predict_kinases_network(**kwargs)
        
        # Get updated similar sites
        similar_sites = get_similar_sites(site_id, **{k: v for k, v in kwargs.items() if k != 'top_n'})
        
        # Get updated heatmap data with enhanced function
        try:
            heatmap_data = get_network_heatmap_data(**kwargs)
            
            # Log heatmap data stats to debug any issues
            logger.info(f"Heatmap data: {len(heatmap_data['sites'])} sites, {len(heatmap_data['kinases'])} kinases, {len(heatmap_data['scores'])} scores")
            
            # Add family distribution for radar chart
            family_distribution = get_kinase_family_distribution_network(
                site_id, score_type=score_type,
                similarity_threshold=kwargs.get('similarity_threshold', 0.6),
                rmsd_threshold=kwargs.get('rmsd_threshold', 3.0)
            )
            
            # Return updated data with full response
            result = {
                'site_id': site_id,
                'score_type': score_type,
                'top_kinases': network_kinases,
                'heatmap': heatmap_data,
                'kinase_families': family_distribution,
                'site_count': len(similar_sites),
                'rmsd_threshold': kwargs.get('rmsd_threshold'),
                'similarity_threshold': kwargs.get('similarity_threshold')
            }
            
            return jsonify(result)
            
        except Exception as heatmap_error:
            # Return basic data if heatmap fails
            logger.error(f"Error generating heatmap: {heatmap_error}")
            
            return jsonify({
                'site_id': site_id,
                'score_type': score_type,
                'top_kinases': network_kinases,
                'heatmap': {
                    'sites': [site_id],  # Include at least the query site
                    'kinases': [k['kinase'] for k in network_kinases[:top_n]],
                    'scores': []  # Empty scores as fallback
                },
                'site_count': len(similar_sites),
                'rmsd_threshold': kwargs.get('rmsd_threshold'),
                'similarity_threshold': kwargs.get('similarity_threshold'),
                'error': f"Heatmap generation failed: {str(heatmap_error)}"
            })
            
    except Exception as e:
        logger.error(f"Error updating network kinase predictions: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return jsonify({
            'error': str(e),
            'traceback': traceback.format_exc()
        }), 500


@app.route('/faq')
def faq():
    """Render the FAQ page."""
    return render_template('faq.html')

if __name__ == '__main__':
    app.run(debug=True)