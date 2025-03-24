// Structural Analysis Tabs Script
document.addEventListener('DOMContentLoaded', function() {
    // Initialize tab switching
    document.querySelectorAll('#structuralTab button').forEach(button => {
        button.addEventListener('click', function(e) {
            e.preventDefault();
            const tabTarget = this.getAttribute('data-bs-target');
            
            // Hide all tabs
            document.querySelectorAll('#structuralTabContent .tab-pane').forEach(tab => {
                tab.classList.remove('show', 'active');
            });
            
            // Show target tab
            document.querySelector(tabTarget).classList.add('show', 'active');
            
            // Update active tab button
            document.querySelectorAll('#structuralTab button').forEach(btn => {
                btn.classList.remove('active');
                btn.setAttribute('aria-selected', 'false');
            });
            
            this.classList.add('active');
            this.setAttribute('aria-selected', 'true');
            
            // If switching to network tab, trigger a resize to ensure proper rendering
            if (tabTarget === '#network' && window.dispatchEvent) {
                window.dispatchEvent(new Event('resize'));
            }
            
            // If switching to enrichment tab, make sure the charts are rendered
            if (tabTarget === '#enrichment') {
                // Wait a bit for the tab to be fully visible before rendering charts
                setTimeout(() => {
                    const event = new Event('resize');
                    window.dispatchEvent(event);
                    
                    // Re-render terminus enrichment charts if needed
                    if (typeof analyzeMotifs === 'function') {
                        try {
                            const dataScript = document.getElementById('terminus-enrichment-data');
                            if (dataScript) {
                                const structuralMatches = JSON.parse(dataScript.textContent);
                                analyzeMotifs(structuralMatches);
                            }
                        } catch (e) {
                            console.error("Error re-rendering terminus enrichment charts:", e);
                        }
                    }
                }, 200);
            }
        });
    });
    
    // Format motifs when tabs change
    document.querySelectorAll('#structuralTab button').forEach(function(button) {
        button.addEventListener('shown.bs.tab', function() {
            // Small delay to ensure content is rendered
            setTimeout(formatMotifBlocks, 200);
        });
    });
});

// Function to format motif blocks with color coding
function formatMotifBlocks() {
    document.querySelectorAll('.motif-sequence:not(.formatted)').forEach(function(element) {
        // Skip elements already formatted or those containing formatted sub-elements
        if (element.classList.contains('formatted') || element.querySelector('.motif-aa')) {
            return;
        }
        
        // Get the motif text 
        const motifText = element.textContent.trim();
        if (!motifText || motifText.length < 7) return;
        
        // Only format strings that look like amino acid sequences
        if (!/^[ACDEFGHIKLMNPQRSTVWYX]+$/.test(motifText)) return;
        
        // Create a flex container
        const container = document.createElement('div');
        container.className = 'motif-flex-container';
        container.style.display = 'flex';
        container.style.flexWrap = 'nowrap';
        
        // Determine center position (phosphosite)
        const centerPos = Math.floor(motifText.length / 2);
        
        // Process each character in the motif
        for (let i = 0; i < motifText.length; i++) {
            const aa = motifText[i];
            const box = document.createElement('div');
            
            // Basic styling
            box.className = 'motif-aa';
            box.style.width = '24px';
            box.style.height = '24px';
            box.style.display = 'flex';
            box.style.alignItems = 'center';
            box.style.justifyContent = 'center';
            box.style.margin = '0 1px';
            box.style.borderRadius = '3px';
            box.textContent = aa;
            
            // Apply colors based on amino acid type
            if (i === centerPos) {
                // Phosphosite (center position)
                box.style.backgroundColor = '#ff5722';
                box.style.color = 'white';
                box.style.fontWeight = 'bold';
                box.classList.add('highlighted');
            } else if (aa === 'X') {
                // Placeholder X
                box.style.backgroundColor = '#e0e0e0';
                box.style.color = '#9e9e9e';
                box.classList.add('aa-x');
            } else if ('STY'.includes(aa)) {
                // STY group
                box.style.backgroundColor = '#bbdefb';
                box.classList.add('sty');
            } else if ('NQ'.includes(aa)) {
                // NQ group
                box.style.backgroundColor = '#b39ddb';
                box.classList.add('nq');
            } else if (aa === 'C') {
                // Cysteine
                box.style.backgroundColor = '#ffcc80';
                box.classList.add('cys');
            } else if (aa === 'P') {
                // Proline
                box.style.backgroundColor = '#81c784';
                box.classList.add('proline');
            } else if ('AVILMFWG'.includes(aa)) {
                // Other nonpolar
                box.style.backgroundColor = '#ffecb3';
                box.classList.add('nonpolar');
            } else if ('DE'.includes(aa)) {
                // Acidic
                box.style.backgroundColor = '#ffcdd2';
                box.classList.add('acidic');
            } else if ('KRH'.includes(aa)) {
                // Basic
                box.style.backgroundColor = '#c8e6c9';
                box.classList.add('basic');
            } else {
                // Special/other
                box.style.backgroundColor = '#e1bee7';
                box.classList.add('special');
            }
            
            container.appendChild(box);
        }
        
        // Replace the original content
        element.innerHTML = '';
        element.appendChild(container);
        element.classList.add('formatted');
    });
}