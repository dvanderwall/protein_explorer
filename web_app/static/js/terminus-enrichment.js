// Terminal Residue Enrichment Analysis Script
document.addEventListener('DOMContentLoaded', function() {
    // Check if we're on the appropriate tab/page
    if (!document.getElementById('terminus-enrichment-data') || 
        !document.getElementById('n-terminus-chart') || 
        !document.getElementById('c-terminus-chart')) {
        return;
    }
    
    // Get structural matches data
    let structuralMatches = [];
    try {
        const dataScript = document.getElementById('terminus-enrichment-data');
        structuralMatches = JSON.parse(dataScript.textContent);
        console.log(`Loaded ${structuralMatches.length} structural matches for enrichment analysis`);
    } catch (e) {
        console.error("Error loading structural matches data:", e);
        return;
    }
    
    if (!structuralMatches || structuralMatches.length === 0) {
        console.log("No structural matches available for enrichment analysis");
        return;
    }
    
    // Analyze motifs from matches
    analyzeMotifs(structuralMatches);
});

function analyzeMotifs(matches) {
    // Extract motifs from matches
    const motifs = matches
        .filter(match => match.motif && typeof match.motif === 'string')
        .map(match => match.motif);
    
    if (motifs.length === 0) {
        console.log("No motifs available for analysis");
        document.getElementById('n-terminus-chart').innerHTML = 
            '<div class="alert alert-info">No motif data available for enrichment analysis.</div>';
        document.getElementById('c-terminus-chart').innerHTML = 
            '<div class="alert alert-info">No motif data available for enrichment analysis.</div>';
        return;
    }
    
    console.log(`Found ${motifs.length} motifs for enrichment analysis`);
    
    // Amino acid group classification
    const aaGroups = {
        'polar': 'STYCNQ',
        'nonpolar': 'AVILMFWPG',
        'acidic': 'DE',
        'basic': 'KRH'
    };
    
    // Background frequencies in the proteome
    const backgroundFreq = {
        'polar': 0.30,      // ~30% of amino acids are polar
        'nonpolar': 0.50,   // ~50% of amino acids are nonpolar
        'acidic': 0.10,     // ~10% of amino acids are acidic
        'basic': 0.10       // ~10% of amino acids are basic
    };
    
    // Process each motif to extract N and C terminal regions
    let nTermData = { 'polar': 0, 'nonpolar': 0, 'acidic': 0, 'basic': 0 };
    let cTermData = { 'polar': 0, 'nonpolar': 0, 'acidic': 0, 'basic': 0 };
    let nTermTotal = 0;
    let cTermTotal = 0;
    
    motifs.forEach(motif => {
        // Standardize motif length by assuming phosphosite is in the middle
        const centerPos = Math.floor(motif.length / 2);
        
        // Extract N-terminal region (-5 to -1)
        const nTermRegion = motif.substring(Math.max(0, centerPos - 5), centerPos);
        // Extract C-terminal region (+1 to +5)
        const cTermRegion = motif.substring(centerPos + 1, Math.min(motif.length, centerPos + 6));
        
        // Process N-terminal region
        for (let i = 0; i < nTermRegion.length; i++) {
            const aa = nTermRegion[i];
            if (aa === 'X') continue; // Skip placeholder X
            
            nTermTotal++;
            
            // Classify amino acid
            for (const [group, aas] of Object.entries(aaGroups)) {
                if (aas.includes(aa)) {
                    nTermData[group]++;
                    break;
                }
            }
        }
        
        // Process C-terminal region
        for (let i = 0; i < cTermRegion.length; i++) {
            const aa = cTermRegion[i];
            if (aa === 'X') continue; // Skip placeholder X
            
            cTermTotal++;
            
            // Classify amino acid
            for (const [group, aas] of Object.entries(aaGroups)) {
                if (aas.includes(aa)) {
                    cTermData[group]++;
                    break;
                }
            }
        }
    });
    
    // Calculate enrichment scores
    const nTermEnrichment = {};
    const cTermEnrichment = {};
    
    for (const group of Object.keys(aaGroups)) {
        // N-terminal enrichment
        const nTermFreq = nTermData[group] / nTermTotal;
        nTermEnrichment[group] = nTermFreq / backgroundFreq[group];
        
        // C-terminal enrichment
        const cTermFreq = cTermData[group] / cTermTotal;
        cTermEnrichment[group] = cTermFreq / backgroundFreq[group];
    }
    
    // Create N-term enrichment chart
    createEnrichmentChart('n-terminus-chart', nTermEnrichment, nTermData);
    
    // Create C-term enrichment chart
    createEnrichmentChart('c-terminus-chart', cTermEnrichment, cTermData);
}

function createEnrichmentChart(elementId, enrichmentData, countData) {
    const colors = {
        'polar': '#bbdefb',     // Light blue
        'nonpolar': '#ffecb3',  // Light yellow/gold
        'acidic': '#ffcdd2',    // Light red
        'basic': '#c8e6c9'      // Light green
    };
    
    const labels = {
        'polar': 'Polar',
        'nonpolar': 'Non-polar',
        'acidic': 'Acidic',
        'basic': 'Basic'
    };
    
    // Prepare data for the chart
    const data = [];
    for (const [group, score] of Object.entries(enrichmentData)) {
        data.push({
            group: labels[group],
            score: score,
            count: countData[group],
            color: colors[group]
        });
    }
    
    // Sort by enrichment score (descending)
    data.sort((a, b) => b.score - a.score);
    
    // Create SVG
    const container = d3.select(`#${elementId}`);
    container.selectAll("*").remove(); // Clear existing content
    
    const margin = {top: 20, right: 20, bottom: 50, left: 60};
    const width = container.node().getBoundingClientRect().width - margin.left - margin.right;
    const height = 250 - margin.top - margin.bottom;
    
    const svg = container.append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
      .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);
    
    // X scale
    const x = d3.scaleBand()
        .domain(data.map(d => d.group))
        .range([0, width])
        .padding(0.2);
    
    // Y scale
    const y = d3.scaleLinear()
        .domain([0, Math.max(2, Math.ceil(d3.max(data, d => d.score)))])
        .range([height, 0]);
    
    // Add X axis
    svg.append("g")
        .attr("transform", `translate(0,${height})`)
        .call(d3.axisBottom(x))
        .selectAll("text")
        .attr("transform", "translate(-10,0)rotate(-45)")
        .style("text-anchor", "end");
    
    // Add Y axis
    svg.append("g")
        .call(d3.axisLeft(y));
    
    // Add Y axis label
    svg.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", -margin.left + 15)
        .attr("x", -height / 2)
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text("Enrichment Score");
    
    // Add reference line at 1.0 (no enrichment)
    svg.append("line")
        .attr("x1", 0)
        .attr("x2", width)
        .attr("y1", y(1))
        .attr("y2", y(1))
        .attr("stroke", "#888")
        .attr("stroke-dasharray", "4,4")
        .attr("stroke-width", 1);
        
    // Add bars
    svg.selectAll(".bar")
        .data(data)
        .enter()
        .append("rect")
        .attr("class", "bar")
        .attr("x", d => x(d.group))
        .attr("y", d => y(Math.max(0, d.score)))
        .attr("width", x.bandwidth())
        .attr("height", d => Math.abs(y(d.score) - y(1)))
        .attr("fill", d => d.color)
        .attr("opacity", 0.8)
        .attr("transform", d => d.score < 1 ? `translate(0,${y(1)})` : "translate(0,0)");
    
    // Add count labels
    svg.selectAll(".count-label")
        .data(data)
        .enter()
        .append("text")
        .attr("class", "count-label")
        .attr("x", d => x(d.group) + x.bandwidth()/2)
        .attr("y", d => d.score >= 1 ? y(d.score) - 5 : y(d.score) + 15)
        .attr("text-anchor", "middle")
        .style("font-size", "12px")
        .style("fill", d => d.score >= 1 ? "#333" : "#333")
        .text(d => `n=${d.count}`);
        
    // Add value labels
    svg.selectAll(".value-label")
        .data(data)
        .enter()
        .append("text")
        .attr("class", "value-label")
        .attr("x", d => x(d.group) + x.bandwidth()/2)
        .attr("y", d => d.score >= 1 ? y(d.score) - 20 : y(d.score) + 30)
        .attr("text-anchor", "middle")
        .style("font-size", "11px")
        .style("font-weight", "bold")
        .style("fill", "#333")
        .text(d => d.score.toFixed(2));
}