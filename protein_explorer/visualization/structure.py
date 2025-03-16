"""
Functions for visualizing protein structures.
"""
import numpy as np
import base64
import plotly.graph_objects as go
import plotly.express as px
import logging
from typing import Dict, List, Tuple, Optional

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def visualize_contact_map(contact_map: np.ndarray, 
                        title: str = "Protein Contact Map") -> go.Figure:
    """
    Visualize a protein contact map.
    
    Args:
        contact_map: NxN numpy array of contacts
        title: Plot title
        
    Returns:
        Plotly figure object
    """
    fig = go.Figure(data=go.Heatmap(
        z=contact_map,
        colorscale=[[0, 'white'], [1, 'blue']],
        showscale=False
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="Residue Index",
        yaxis_title="Residue Index",
        width=600,
        height=600,
        xaxis=dict(scaleanchor="y", scaleratio=1),
        yaxis=dict(autorange='reversed')
    )
    
    return fig

def visualize_structure(pdb_data: str, 
                       container_id: str = "ngl-viewer",
                       height: int = 500,
                       sequence: str = None) -> str:
    """
    Create an enhanced 3D visualization of protein structure using NGL Viewer.
    Shows side chains only for selected residue and residues within 5Å.
    Colors protein by PLDDT score for confidence visualization (red=high, blue=low).
    
    Args:
        pdb_data: PDB format data as string
        container_id: HTML container ID
        height: Height of the viewer in pixels
        sequence: Protein sequence string (optional)
        
    Returns:
        HTML and JavaScript code for the viewer
    """
    # Base64 encode the PDB data
    pdb_base64 = base64.b64encode(pdb_data.encode()).decode()
    
    # Get sequence from protein data if not provided
    if not sequence:
        import re
        try:
            # Try to extract from PDB data
            sequence_match = re.search(r'SEQRES.*?\n(.*?)ENDMDL', pdb_data, re.DOTALL)
            if sequence_match:
                sequence = ''.join(sequence_match.group(1).split())
            else:
                sequence = ""
        except:
            sequence = ""
    
    # Create the JavaScript code for NGL Viewer with PLDDT coloring and sequence navigation
    js_code = f"""
    <style>
        .sequence-viewer {{
            font-family: monospace;
            line-height: 1.5;
            font-size: 14px;
            overflow-x: auto;
            white-space: nowrap;
            margin-bottom: 10px;
            padding: 8px;
            border: 1px solid #ddd;
            border-radius: 5px;
            background-color: #f8f9fa;
        }}
        .sequence-viewer span {{
            cursor: pointer;
            padding: 2px 1px;
        }}
        .sequence-viewer span:hover {{
            background-color: #f0f0f0;
        }}
        .sequence-viewer span.highlighted {{
            background-color: #4285f4;
            color: white;
            font-weight: bold;
        }}
        .sequence-viewer span.temp-highlight {{
            background-color: #ffcc00;
        }}
        .aa-info-panel {{
            position: absolute;
            top: 10px;
            right: 10px;
            background-color: rgba(255, 255, 255, 0.9);
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 8px;
            font-size: 14px;
            font-family: sans-serif;
            z-index: 100;
            display: none;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .controls-panel {{
            position: absolute;
            bottom: 10px;
            left: 10px;
            z-index: 100;
        }}
        .controls-panel button {{
            background-color: rgba(255, 255, 255, 0.9);
            border: 1px solid #ddd;
            border-radius: 3px;
            padding: 5px 10px;
            margin-right: 5px;
            font-size: 12px;
            cursor: pointer;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }}
        .controls-panel button:hover {{
            background-color: rgba(240, 240, 240, 0.9);
        }}
        .color-legend {{
            position: absolute;
            bottom: 10px;
            right: 10px;
            background-color: rgba(255, 255, 255, 0.9);
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 8px;
            font-size: 12px;
            z-index: 100;
            width: 170px;
        }}
        .legend-gradient {{
            height: 15px;
            background: linear-gradient(to right, #0000FF, #00FFFF, #00FF00, #FFFF00, #FF0000);
            margin-top: 5px;
            margin-bottom: 2px;
        }}
        .legend-labels {{
            display: flex;
            justify-content: space-between;
            font-size: 10px;
        }}
    </style>
    
    <div class="sequence-viewer" id="sequence-{container_id}">
    {sequence if sequence else "Loading sequence..."}
    </div>
    
    <div style="position: relative; width: 100%; height: {height}px;">
        <div id="{container_id}" style="width: 100%; height: 100%;"></div>
        <div class="aa-info-panel" id="info-panel-{container_id}">
            <strong>Residue:</strong> <span id="residue-info"></span><br>
            <strong>Chain:</strong> <span id="chain-info"></span><br>
            <strong>PLDDT:</strong> <span id="plddt-info">N/A</span>
        </div>
        <div class="controls-panel">
            <button id="reset-view-btn">Reset View</button>
            <button id="toggle-sidechains-btn">Hide Sidechains</button>
            <button id="toggle-color-mode-btn">Color by Chain</button>
        </div>
        <div class="color-legend">
            <strong>PLDDT Score</strong>
            <div class="legend-gradient"></div>
            <div class="legend-labels">
                <span>High</span>
                <span>Medium</span>
                <span>Low</span>
            </div>
        </div>
    </div>
    
    <script src="https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.37/dist/ngl.js"></script>
    <script type="text/javascript">
        document.addEventListener('DOMContentLoaded', function() {{
            // Create NGL Stage object
            var stage = new NGL.Stage('{container_id}', {{backgroundColor: "white"}});
            
            // Handle window resizing
            window.addEventListener('resize', function() {{
                stage.handleResize();
            }}, false);
            
            // Get DOM elements
            var infoPanel = document.getElementById('info-panel-{container_id}');
            var residueInfo = document.getElementById('residue-info');
            var chainInfo = document.getElementById('chain-info');
            var plddtInfo = document.getElementById('plddt-info');
            var resetViewBtn = document.getElementById('reset-view-btn');
            var toggleSidechainsBtn = document.getElementById('toggle-sidechains-btn');
            var toggleColorModeBtn = document.getElementById('toggle-color-mode-btn');
            var sequenceViewer = document.getElementById('sequence-{container_id}');
            
            // Variables to track state
            var component = null;
            var sidechainsVisible = true;
            var mainRepresentation = null;
            var colorMode = "plddt"; // Start with PLDDT coloring
            var currentHighlightedResidues = [];
            var sequenceSpans = [];
            
            // Amino acid 3-letter to 1-letter code mapping
            var aaMap = {{
                'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
                'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
            }};
            
            // Function to remove highlights
            function removeHighlights() {{
                // Remove all highlighted side chains
                component.removeAllRepresentations();
                
                // Add back main cartoon representation with current color mode
                updateMainRepresentation();
                
                // Reset sequence highlighting
                sequenceSpans.forEach(function(span) {{
                    span.classList.remove('highlighted');
                    span.classList.remove('temp-highlight');
                }});
                
                currentHighlightedResidues = [];
            }}
            
            // Function to update main representation with current color mode
            function updateMainRepresentation() {{
                if (colorMode === "plddt") {{
                    // Color by PLDDT (using B-factor field in PDB)
                    // Red is high confidence, blue is low confidence
                    mainRepresentation = component.addRepresentation("cartoon", {{
                        color: "bfactor",
                        colorScale: "rwb",
                        colorReverse: false,
                        colorDomain: [50, 100], // PLDDT range is typically 50-100
                        smoothSheet: true
                    }});
                }} else {{
                    // Color by chain
                    mainRepresentation = component.addRepresentation("cartoon", {{
                        color: "chainid",
                        smoothSheet: true
                    }});
                }}
            }}
            
            // Function to toggle color mode
            function toggleColorMode() {{
                colorMode = colorMode === "plddt" ? "chain" : "plddt";
                toggleColorModeBtn.innerText = colorMode === "plddt" ? "Color by Chain" : "Color by PLDDT";
                removeHighlights();
            }}
            
            // Function to toggle sidechains visibility
            function toggleSidechains() {{
                sidechainsVisible = !sidechainsVisible;
                toggleSidechainsBtn.innerText = sidechainsVisible ? 'Hide Sidechains' : 'Show Sidechains';
                if (!sidechainsVisible) {{
                    removeHighlights();
                }}
            }}
            
            // Function to highlight residues in the sequence
            function highlightSequenceResidue(resno, isTemp = false) {{
                // Convert from PDB residue number to zero-based index if needed
                var index = resno - 1; // Adjust if your PDB doesn't start at 1
                
                if (index >= 0 && index < sequenceSpans.length) {{
                    sequenceSpans[index].classList.add(isTemp ? 'temp-highlight' : 'highlighted');
                    if (!isTemp) {{
                        currentHighlightedResidues.push(resno);
                        
                        // Scroll to make highlighted residue visible
                        sequenceSpans[index].scrollIntoView({{
                            behavior: 'smooth',
                            block: 'nearest',
                            inline: 'center'
                        }});
                    }}
                }}
            }}
            
            // Format sequence with numbers
            function formatSequence(sequence) {{
                var html = '';
                for (var i = 0; i < sequence.length; i++) {{
                    var resno = i + 1;
                    html += '<span data-resno="' + resno + '">' + sequence[i] + '</span>';
                    
                    // Add residue number every 10 residues
                    if ((i + 1) % 10 === 0) {{
                        html += '<sub>' + (i + 1) + '</sub> ';
                    }}
                }}
                return html;
            }}
            
            // Load PDB data from base64 string
            var pdbBlob = new Blob([atob('{pdb_base64}')], {{type: 'text/plain'}});
            
            // Load the structure
            stage.loadFile(pdbBlob, {{ext: 'pdb'}}).then(function (comp) {{
                component = comp;
                var structure = component.structure;
                
                // Extract PLDDT data (stored as B-factors in AlphaFold PDBs)
                var plddtByResidue = {{}};
                structure.eachResidue(function(rp) {{
                    var bFactorSum = 0;
                    var atomCount = 0;
                    
                    // Calculate average B-factor for the residue
                    structure.eachAtom(function(ap) {{
                        if (ap.resno === rp.resno && ap.chainname === rp.chainname) {{
                            bFactorSum += ap.bfactor;
                            atomCount++;
                        }}
                    }});
                    
                    if (atomCount > 0) {{
                        plddtByResidue[rp.resno] = bFactorSum / atomCount;
                    }}
                }});
                
                // Add main representation with PLDDT coloring
                updateMainRepresentation();
                
                // Process structure to extract sequence and set up the sequence viewer
                var proteinSequence = "{sequence}";
                
                if (!proteinSequence) {{
                    // Extract sequence from structure
                    var chainSeq = {{}};
                    
                    structure.eachResidue(function(rp) {{
                        var chain = rp.chainname || rp.chainid;
                        if (!chainSeq[chain]) chainSeq[chain] = [];
                        
                        var resno = rp.resno;
                        var resname = rp.resname;
                        var aa = aaMap[resname] || 'X';
                        
                        chainSeq[chain][resno] = aa;
                    }});
                    
                    // Use the first chain's sequence
                    var mainChain = Object.keys(chainSeq)[0];
                    proteinSequence = Object.values(chainSeq[mainChain]).join('').replace(/undefined/g, '');
                }}
                
                // Format and display the sequence
                sequenceViewer.innerHTML = formatSequence(proteinSequence);
                
                // Get all sequence spans
                sequenceSpans = sequenceViewer.querySelectorAll('span');
                
                // Add click event to sequence spans
                sequenceSpans.forEach(function(span) {{
                    span.addEventListener('click', function() {{
                        var resno = parseInt(this.getAttribute('data-resno'));
                        showSideChainsForResidue(resno);
                    }});
                }});
                
                // Function to show side chains for a residue and nearby residues
                function showSideChainsForResidue(resno) {{
                    if (!sidechainsVisible) return;
                    
                    // Remove previous highlights
                    removeHighlights();
                    
                    // Add back main representation with current color mode
                    updateMainRepresentation();
                    
                    // Selection for the clicked residue
                    var selection = resno;
                    
                    // Create a selection for residues within 5Å
                    var withinSelection = selection + " or (" + selection + " around 5)";
                    
                    // Add side chain representation for selected residue and nearby residues
                    component.addRepresentation("licorice", {{
                        sele: withinSelection + " and sidechainAttached",
                        color: "element",
                        opacity: 1.0,
                        multipleBond: "symmetric"
                    }});
                    
                    // Highlight in sequence
                    highlightSequenceResidue(resno);
                    
                    // Find and highlight nearby residues in sequence
                    var withinSel = new NGL.Selection(withinSelection);
                    structure.eachResidue(function(rp) {{
                        if (withinSel.test(rp) && rp.resno !== resno) {{
                            highlightSequenceResidue(rp.resno);
                        }}
                    }});
                    
                    // Focus view on the selection
                    component.autoView(withinSelection, 2000);
                }}
                
                // Handle hover events to show residue info
                stage.signals.hovered.add(function(pickingData) {{
                    // Remove all temporary highlights
                    sequenceSpans.forEach(span => span.classList.remove('temp-highlight'));
                    
                    if (pickingData.atom) {{
                        var atom = pickingData.atom;
                        var resname = atom.resname;
                        var resno = atom.resno;
                        var chain = atom.chainname || atom.chainid;
                        var aa = aaMap[resname] || resname;
                        
                        // Get PLDDT score for residue
                        var plddt = plddtByResidue[resno];
                        var plddtText = plddt ? plddt.toFixed(1) : "N/A";
                        
                        // Show info panel
                        residueInfo.innerText = aa + ' ' + resno;
                        chainInfo.innerText = chain;
                        plddtInfo.innerText = plddtText;
                        infoPanel.style.display = 'block';
                        
                        // Update cursor
                        document.body.style.cursor = 'pointer';
                        
                        // Highlight in sequence with temp highlight
                        highlightSequenceResidue(resno, true);
                        
                    }} else {{
                        // Hide info panel when not hovering over an atom
                        infoPanel.style.display = 'none';
                        document.body.style.cursor = 'default';
                    }}
                }});
                
                // Handle click events
                stage.signals.clicked.add(function(pickingData) {{
                    if (pickingData.atom) {{
                        var atom = pickingData.atom;
                        var resno = atom.resno;
                        
                        showSideChainsForResidue(resno);
                    }}
                }});
                
                // Reset view button
                resetViewBtn.addEventListener('click', function() {{
                    removeHighlights();
                    component.autoView();
                }});
                
                // Toggle sidechains button
                toggleSidechainsBtn.addEventListener('click', function() {{
                    toggleSidechains();
                }});
                
                // Toggle color mode button
                toggleColorModeBtn.addEventListener('click', function() {{
                    toggleColorMode();
                }});
                
                // Zoom to full structure
                component.autoView();
            }});
        }});
    </script>
    """
    
    return js_code

def compare_structures(structure_list: List[str], 
                     container_id: str = "molstar-comparison",
                     height: int = 500) -> str:
    """
    Create a visualization comparing multiple protein structures.
    
    Args:
        structure_list: List of PDB format structures as strings
        container_id: HTML container ID
        height: Height of the viewer in pixels
        
    Returns:
        HTML and JavaScript code for the viewer
    """
    # Base64 encode the PDB data
    pdb_base64_list = [base64.b64encode(s.encode()).decode() for s in structure_list]
    
    # Create the JavaScript code for MolStar
    js_code = f"""
    <div id="{container_id}" style="width:100%; height:{height}px;"></div>
    <script src="https://www.alphafold.ebi.ac.uk/assets/js/molstar.js"></script>
    <script type="text/javascript">
        document.addEventListener('DOMContentLoaded', function () {{
            // Make sure molstar is loaded
            if (typeof molstar === 'undefined') {{
                console.error('MolStar library not loaded');
                document.getElementById('{container_id}').innerHTML = '<div class="alert alert-danger">Error loading MolStar viewer. Please check console for details.</div>';
                return;
            }}
            
            console.log('Creating MolStar comparison viewer');
            molstar.Viewer.create('{container_id}', {{
                layoutIsExpanded: false,
                layoutShowControls: true,
                layoutShowRemoteState: false,
                layoutShowSequence: true,
                layoutShowLog: false,
                layoutShowLeftPanel: true,
                layoutShowStructureSourceControls: true,
                viewportShowExpand: true,
                viewportShowSelectionMode: true,
                viewportShowAnimation: true,
            }}).then(viewer => {{
                // Load all structures
                const structures = {pdb_base64_list};
                
                // Load the first structure
                const data1 = new Uint8Array(atob(structures[0]).split('').map(c => c.charCodeAt(0)));
                console.log('Loading first structure, length: ' + data1.length);
                viewer.loadStructureFromData(data1, 'pdb', {{
                    representationParams: {{
                        theme: {{ globalName: 'chain-id' }},
                        type: 'cartoon',
                    }}
                }}).then(() => {{
                    // Load additional structures
                    for (let i = 1; i < structures.length; i++) {{
                        const data = new Uint8Array(atob(structures[i]).split('').map(c => c.charCodeAt(0)));
                        console.log('Loading additional structure ' + i + ', length: ' + data.length);
                        viewer.loadStructureFromData(data, 'pdb', {{
                            representationParams: {{
                                theme: {{ globalName: 'chain-id' }},
                                type: 'cartoon',
                            }}
                        }}).catch(error => {{
                            console.error('Error loading structure ' + i + ':', error);
                        }});
                    }}
                }}).catch(error => {{
                    console.error('Error loading first structure:', error);
                }});
            }}).catch(error => {{
                console.error('Error creating comparison viewer:', error);
                document.getElementById('{container_id}').innerHTML = '<div class="alert alert-danger">Error creating MolStar viewer: ' + error.message + '</div>';
            }});
        }});
    </script>
    """
    
    return js_code

def visualize_pca_results(eigenvalues: np.ndarray, 
                        projected_coords: np.ndarray) -> Dict:
    """
    Visualize PCA results for a protein structure.
    
    Args:
        eigenvalues: Eigenvalues from PCA
        projected_coords: Coordinates projected onto principal components
        
    Returns:
        Dictionary with Plotly figure objects
    """
    # Calculate explained variance ratio
    explained_var = eigenvalues / np.sum(eigenvalues)
    
    # Create scree plot (eigenvalues)
    scree_fig = go.Figure(data=[
        go.Bar(
            x=list(range(1, len(eigenvalues) + 1)),
            y=explained_var,
            marker_color='royalblue'
        )
    ])
    
    scree_fig.update_layout(
        title="PCA Scree Plot",
        xaxis_title="Principal Component",
        yaxis_title="Explained Variance Ratio",
        xaxis=dict(tickmode='linear'),
        yaxis=dict(range=[0, max(explained_var) * 1.1])
    )
    
    # Create 2D projection plot (PC1 vs PC2)
    projection_fig = go.Figure(data=[
        go.Scatter(
            x=projected_coords[:, 0],
            y=projected_coords[:, 1],
            mode='markers',
            marker=dict(
                size=8,
                color=list(range(len(projected_coords))),
                colorscale='Viridis',
                colorbar=dict(title="Residue Index")
            )
        )
    ])
    
    projection_fig.update_layout(
        title="PCA Projection (PC1 vs PC2)",
        xaxis_title=f"PC1 ({explained_var[0]:.2%})",
        yaxis_title=f"PC2 ({explained_var[1]:.2%})"
    )
    
    # Create 3D projection plot (PC1 vs PC2 vs PC3) if we have at least 3 components
    if projected_coords.shape[1] >= 3:
        projection_3d_fig = go.Figure(data=[
            go.Scatter3d(
                x=projected_coords[:, 0],
                y=projected_coords[:, 1],
                z=projected_coords[:, 2],
                mode='markers',
                marker=dict(
                    size=4,
                    color=list(range(len(projected_coords))),
                    colorscale='Viridis',
                    colorbar=dict(title="Residue Index")
                )
            )
        ])
        
        projection_3d_fig.update_layout(
            title="PCA 3D Projection",
            scene=dict(
                xaxis_title=f"PC1 ({explained_var[0]:.2%})",
                yaxis_title=f"PC2 ({explained_var[1]:.2%})",
                zaxis_title=f"PC3 ({explained_var[2]:.2%})"
            )
        )
    else:
        projection_3d_fig = None
    
    return {
        "scree_plot": scree_fig,
        "projection_2d": projection_fig,
        "projection_3d": projection_3d_fig if projected_coords.shape[1] >= 3 else None
    }


def visualize_phosphosites(phosphosites: List[Dict], 
                          container_id: str = "phosphosite-table") -> str:
    """
    Create an HTML visualization of phosphorylation sites.
    
    Args:
        phosphosites: List of dictionaries with phosphosite information
        container_id: HTML container ID
        
    Returns:
        HTML for displaying the phosphorylation sites table
    """
    # Generate HTML for the phosphorylation site table
    html = f"""
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
                    <tbody id="{container_id}">
    """
    
    # Add rows for each phosphosite
    for site in phosphosites:
        html += f"""
        <tr>
            <td><a href="#" class="site-link" data-resno="{site['resno']}">{site['site']}</a></td>
            <td><code class="motif-sequence">{site['motif']}</code></td>
            <td>{site['mean_plddt']}</td>
            <td>{site['nearby_count']}</td>
            <td>{"Yes" if site['is_known'] else "No"}</td>
        </tr>
        """
    
    # Close the table and add JavaScript for interaction
    html += """
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
    
    return html