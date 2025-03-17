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
from protein_explorer.analysis.phospho import analyze_phosphosites

# Add the parent directory to the path to import protein_explorer
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import protein_explorer as pe


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
    """Display protein information and structure."""
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
                            phosphosites = analyze_phosphosites(
                                sequence, 
                                structure, 
                                uniprot_id=protein_data.get('uniprot_id')
                            )
                            print(f"DEBUG: Found {len(phosphosites)} potential phosphorylation sites")
                            
                            # Generate phosphosite HTML
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
                            
                            # Add rows for each phosphosite
                            for site in phosphosites:
                                phosphosite_html += f"""
                                <tr>
                                    <td><a href="#" class="site-link" data-resno="{site['resno']}">{site['site']}</a></td>
                                    <td><code class="motif-sequence">{site['motif']}</code></td>
                                    <td>{site['mean_plddt']}</td>
                                    <td>{site['nearby_count']}</td>
                                    <td>{"Yes" if site['is_known'] else "No"}</td>
                                </tr>
                                """
                            
                            # Close the table
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
                            print(f"DEBUG: Error analyzing phosphosites: {e}")
                            phosphosite_html = f"""
                            <div class="card mt-4">
                                <div class="card-header">
                                    <h5 class="mb-0">Phosphorylation Site Analysis</h5>
                                </div>
                                <div class="card-body">
                                    <div class="alert alert-warning">
                                        Error analyzing phosphorylation sites: {str(e)}
                                    </div>
                                </div>
                            </div>
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
            phosphosites=phosphosites  # Directly pass the phosphosites list
        )
    except Exception as e:
        print(f"DEBUG: Exception in protein view: {e}")
        logger.error(f"Error in protein view: {e}")
        return render_template('error.html', error=str(e))

@app.route('/network/<uniprot_id>')
def network(uniprot_id):
    """Display protein interaction network."""
    try:
        # Get network depth from query parameters
        depth = int(request.args.get('depth', 1))
        
        # Build network
        network = pe.navigation.build_interaction_network([uniprot_id], max_depth=depth)
        
        # Create network visualization
        network_html = pe.visualization.visualize_network(
            network,
            highlight_nodes=[uniprot_id],
            title=f"Protein Interaction Network for {uniprot_id}"
        )
        
        # Get key proteins
        key_proteins = pe.analysis.identify_key_proteins(network)
        
        return render_template(
            'network.html',
            uniprot_id=uniprot_id,
            network_html=network_html,
            key_proteins=key_proteins,
            current_depth=depth,
            node_count=network.number_of_nodes(),
            edge_count=network.number_of_edges()
        )
    except Exception as e:
        logger.error(f"Error in network view: {e}")
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

@app.route('/api/protein/<uniprot_id>', methods=['GET'])
def api_protein(uniprot_id):
    """API endpoint for protein data."""
    try:
        protein_data = pe.data.get_protein_by_id(uniprot_id=uniprot_id)
        return jsonify(protein_data)
    except Exception as e:
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

@app.route('/api/phosphosites/<uniprot_id>', methods=['GET'])
def api_phosphosites(uniprot_id):
    """API endpoint for phosphorylation site analysis."""
    try:
        protein_data = pe.data.get_protein_by_id(uniprot_id=uniprot_id)
        sequence = protein_data.get('metadata', {}).get('sequence', {}).get('value')
        
        if not sequence:
            return jsonify({'error': 'Protein sequence not found'}), 404
            
        structure = pe.data.get_alphafold_structure(uniprot_id)
        if not structure:
            return jsonify({'error': 'Protein structure not found'}), 404
            
        phosphosites = analyze_phosphosites(sequence, structure, uniprot_id)
        return jsonify(phosphosites)
    except Exception as e:
        return jsonify({'error': str(e)}), 400
    
@app.route('/phosphosite', methods=['GET', 'POST'])
def phosphosite_analysis():
    """Phosphosite structural analysis page with improved error handling."""
    import os
    from flask import request, render_template
    
    # Import our phospho_analyzer module
    from protein_explorer.analysis.phospho_analyzer import analyze_protein
    
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
            # Get the parquet file path
            parquet_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 
                                     'Combined_Kinome_10A_Master_Filtered_2.parquet')
            
            # Run the full analysis
            results = analyze_protein(identifier, id_type, parquet_file)
            
            # Return results - now including potential error message
            return render_template('phosphosite.html', 
                                  protein_info=results['protein_info'],
                                  phosphosites=results['phosphosites'],
                                  structural_matches=results['structural_matches'],
                                  error=results.get('error'))
                
        except Exception as e:
            print(f"Error in phosphosite analysis: {e}")
            error = str(e)
            return render_template('phosphosite.html', error=error)
    
    # GET request - show empty form
    return render_template('phosphosite.html')

if __name__ == '__main__':
    app.run(debug=True)