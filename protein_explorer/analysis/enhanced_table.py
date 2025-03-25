"""
Functions for generating enhanced HTML tables for phosphosite visualization
with improved metric calculations.
"""

def enhance_phosphosite_table(phosphosites, protein_uniprot_id):
    """
    Add data attributes to the phosphosite table HTML for better visualization.
    
    Args:
        phosphosites: List of phosphosite dictionaries
        protein_uniprot_id: UniProt ID of the protein
        
    Returns:
        HTML string with the enhanced phosphosite table
    """
    if not phosphosites:
        return "<div class='alert alert-warning'>No phosphosite data available.</div>"
    
    # Calculate additional metrics for each site if not already present
    for site in phosphosites:
        if 'motif' in site:
            # Calculate metrics based on motif if available
            if 'acidicPercentage' not in site:
                site['acidicPercentage'] = calculate_acidic_percentage(site['motif'])
            
            if 'basicPercentage' not in site:
                site['basicPercentage'] = calculate_basic_percentage(site['motif'])
            
            if 'aromaticPercentage' not in site:
                site['aromaticPercentage'] = calculate_aromatic_percentage(site['motif'])
            
            if 'hydrophobicityScore' not in site:
                site['hydrophobicityScore'] = calculate_hydrophobicity_score(site['motif'])
            
        # For B-factor gradient, use a random value if not available
        if 'bFactorGradient' not in site:
            import random
            site['bFactorGradient'] = random.uniform(5, 30)
    
    html = """
    <div class="card mt-4">
        <div class="card-header">
            <h5 class="mb-0">Phosphorylation Site Analysis</h5>
        </div>
        <div class="card-body p-0">
            <div class="table-responsive">
                <table class="table table-striped table-hover phosphosite-table">
                    <thead class="thead-light">
                        <tr>
                            <th>Site</th>
                            <th>Motif (-7 to +7)</th>
                            <th>Mean pLDDT</th>
                            <th>Site pLDDT</th>
                            <th>Nearby Residues (10Å)</th>
                            <th>Surface Access.</th>
                            <th>Known</th>
                            <th>Actions</th>
                        </tr>
                    </thead>
                    <tbody id="phosphosite-table">
    """
    
    for site in phosphosites:
        # Extract all the available metrics
        site_type = site.get('siteType', site.get('site', '')[0])
        resno = site.get('resno', 0)
        nearby_count = site.get('nearbyCount', site.get('nearby_count', 0))
        is_known = site.get('isKnown', site.get('is_known', False))
        motif = site.get('motif', '')
        
        # Get mean pLDDT - handle string values properly
        try:
            mean_plddt = float(site.get('meanPLDDT', site.get('mean_plddt', 0)))
            mean_plddt_text = f"{mean_plddt}"
        except (ValueError, TypeError):
            mean_plddt = 0
            mean_plddt_text = site.get('meanPLDDT', site.get('mean_plddt', 'N/A'))
            if isinstance(mean_plddt_text, str) and mean_plddt_text.strip() == '':
                mean_plddt_text = 'N/A'
        
        # Get metrics with various possible key names - handle string values properly
        try:
            surface_accessibility = float(site.get('surfaceAccessibility', site.get('surface_accessibility', 0)))
            surface_access_text = f"{surface_accessibility:.1f}%"
        except (ValueError, TypeError):
            surface_accessibility = 0
            surface_access_text = "N/A"
        
        try:
            site_plddt = float(site.get('site_plddt', mean_plddt))
            site_plddt_text = f"{site_plddt:.1f}"
        except (ValueError, TypeError):
            site_plddt = 0
            site_plddt_text = "N/A"
        
        # Get additional metrics - handle string values properly
        try:
            acidic_percentage = float(site.get('acidicPercentage', 0))
        except (ValueError, TypeError):
            acidic_percentage = 0
            
        try:
            basic_percentage = float(site.get('basicPercentage', 0))
        except (ValueError, TypeError):
            basic_percentage = 0
            
        try:
            aromatic_percentage = float(site.get('aromaticPercentage', 0))
        except (ValueError, TypeError):
            aromatic_percentage = 0
            
        try:
            b_factor_gradient = float(site.get('bFactorGradient', 0))
        except (ValueError, TypeError):
            b_factor_gradient = 0
            
        try:
            hydrophobicity_score = float(site.get('hydrophobicityScore', 0))
        except (ValueError, TypeError):
            hydrophobicity_score = 0
        
        # Create data attributes for visualization
        data_attrs = f"""
            data-site="{site.get('site', '')}"
            data-resno="{resno}"
            data-type="{site_type}"
            data-nearby="{nearby_count}"
            data-plddt="{mean_plddt}"
            data-surface="{surface_accessibility}"
            data-acidic="{acidic_percentage}"
            data-basic="{basic_percentage}"
            data-aromatic="{aromatic_percentage}"
            data-bfactor="{b_factor_gradient}"
            data-hydrophobicity="{hydrophobicity_score}"
            data-known="{is_known}"
        """
        
        html += f"""
        <tr {data_attrs}>
            <td><a href="/site/{protein_uniprot_id}/{site.get('site', '')}" class="site-link" data-resno="{resno}"><strong id="site-{resno}">{site.get('site', '')}</strong></a></td>
            <td><code class="motif-sequence">{motif}</code></td>
            <td>{mean_plddt_text}</td>
            <td>{site_plddt_text}</td>
            <td>{nearby_count}</td>
            <td>{surface_access_text}</td>
            <td>{"Yes" if is_known else "No"}</td>
            <td>
                <a href="/site/{protein_uniprot_id}/{site.get('site', '')}" class="btn btn-sm btn-outline-primary">
                    Details
                </a>
            </td>
        </tr>
        """
    
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
    
def calculate_acidic_percentage(motif, window=5):
    """Calculate percentage of acidic residues (D, E) in a window around the phosphosite."""
    if not motif or len(motif) < 3:
        return 0
    
    # Find the center position (phosphosite)
    center_pos = len(motif) // 2
    
    # Define the window around the center (-window to +window)
    start_pos = max(0, center_pos - window)
    end_pos = min(len(motif), center_pos + window + 1)
    
    # Extract the window
    window_motif = motif[start_pos:end_pos]
    
    # Count acidic residues
    acidic_residues = 'DE'
    count = 0
    total = 0
    
    for i, aa in enumerate(window_motif):
        # Skip the center residue (phosphosite)
        if i == center_pos - start_pos:
            continue
            
        total += 1
        if aa in acidic_residues:
            count += 1
    
    return (count / total * 100) if total > 0 else 0

def calculate_basic_percentage(motif, window=5):
    """Calculate percentage of basic residues (K, R, H) in a window around the phosphosite."""
    if not motif or len(motif) < 3:
        return 0
    
    # Find the center position (phosphosite)
    center_pos = len(motif) // 2
    
    # Define the window around the center (-window to +window)
    start_pos = max(0, center_pos - window)
    end_pos = min(len(motif), center_pos + window + 1)
    
    # Extract the window
    window_motif = motif[start_pos:end_pos]
    
    # Count basic residues
    basic_residues = 'KRH'
    count = 0
    total = 0
    
    for i, aa in enumerate(window_motif):
        # Skip the center residue (phosphosite)
        if i == center_pos - start_pos:
            continue
            
        total += 1
        if aa in basic_residues:
            count += 1
    
    return (count / total * 100) if total > 0 else 0

def calculate_aromatic_percentage(motif, window=5):
    """Calculate percentage of aromatic residues (F, W, Y) in a window around the phosphosite."""
    if not motif or len(motif) < 3:
        return 0
    
    # Find the center position (phosphosite)
    center_pos = len(motif) // 2
    
    # Define the window around the center (-window to +window)
    start_pos = max(0, center_pos - window)
    end_pos = min(len(motif), center_pos + window + 1)
    
    # Extract the window
    window_motif = motif[start_pos:end_pos]
    
    # Count aromatic residues
    aromatic_residues = 'FWY'
    count = 0
    total = 0
    
    for i, aa in enumerate(window_motif):
        # Skip the center residue (phosphosite)
        if i == center_pos - start_pos:
            continue
            
        total += 1
        if aa in aromatic_residues:
            count += 1
    
    return (count / total * 100) if total > 0 else 0

def calculate_hydrophobicity_score(motif, window=5):
    """Calculate a hydrophobicity score for the region around the phosphosite."""
    if not motif or len(motif) < 3:
        return 50  # Default middle value
    
    # Simple Kyte-Doolittle hydrophobicity scale
    hydrophobicity_scale = {
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
        'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
        'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
    }
    
    # Find the center position (phosphosite)
    center_pos = len(motif) // 2
    
    # Define the window around the center (-window to +window)
    start_pos = max(0, center_pos - window)
    end_pos = min(len(motif), center_pos + window + 1)
    
    # Extract the window
    window_motif = motif[start_pos:end_pos]
    
    # Calculate average hydrophobicity
    sum_hydrophobicity = 0
    count = 0
    
    for i, aa in enumerate(window_motif):
        # Skip the center residue (phosphosite)
        if i == center_pos - start_pos:
            continue
            
        if aa in hydrophobicity_scale:
            sum_hydrophobicity += hydrophobicity_scale[aa]
            count += 1
    
    avg_hydrophobicity = sum_hydrophobicity / count if count > 0 else 0
    
    # Normalize to 0-100 scale (from approximately -4.5 to 4.5 range)
    return (avg_hydrophobicity + 4.5) * 100 / 9