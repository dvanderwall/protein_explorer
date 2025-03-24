// Kinase Prediction Visualization
document.addEventListener('DOMContentLoaded', function() {
    // Initialize structural kinase chart
    initStructuralKinaseCharts();
    
    // Initialize sequence kinase chart
    initSequenceKinaseCharts();
    
    // Initialize combined kinase analysis
    initCombinedKinaseAnalysis();
    
    // Add tab change event listeners to ensure charts are properly sized
    document.querySelectorAll('#structuralTab button, #sequenceTab button').forEach(function(button) {
        button.addEventListener('click', function() {
            setTimeout(function() {
                // Force Chart.js to resize if we're on a kinase tab
                if (button.id === 'struct-kinase-tab' || button.id === 'seq-kinase-tab') {
                    window.dispatchEvent(new Event('resize'));
                }
            }, 100);
        });
    });
});

// Initialize structural kinase charts
function initStructuralKinaseCharts() {
    const dataScript = document.getElementById('struct-kinase-data');
    if (!dataScript) return;
    
    try {
        // Parse kinase data
        const kinaseData = JSON.parse(dataScript.textContent);
        console.log('Structural kinase data:', kinaseData);
        
        // Create bar chart for top kinases
        createKinaseBarChart('struct-kinase-chart', kinaseData.top_kinases, 'structure');
        
        // Create kinase family chart
        createKinaseFamilyChart('struct-kinase-family-chart', kinaseData.top_kinases);
        
        // Create motif analysis
        createMotifAnalysis('struct-motif-analysis', kinaseData.top_kinases);
        
        // Create heatmap for kinase comparison if available
        if (kinaseData.heatmap && document.getElementById('struct-kinase-heatmap')) {
            createKinaseHeatmap('struct-kinase-heatmap', kinaseData.heatmap, 'structure');
        }
    } catch (e) {
        console.error('Error initializing structural kinase charts:', e);
    }
}

// Initialize sequence kinase charts
function initSequenceKinaseCharts() {
    const dataScript = document.getElementById('seq-kinase-data');
    if (!dataScript) return;
    
    try {
        // Parse kinase data
        const kinaseData = JSON.parse(dataScript.textContent);
        console.log('Sequence kinase data:', kinaseData);
        
        // Create bar chart for top kinases
        createKinaseBarChart('seq-kinase-chart', kinaseData.top_kinases, 'sequence');
        
        // Create kinase family chart
        createKinaseFamilyChart('seq-kinase-family-chart', kinaseData.top_kinases);
        
        // Create motif analysis
        createMotifAnalysis('seq-motif-analysis', kinaseData.top_kinases);
        
        // Create heatmap for kinase comparison if available
        if (kinaseData.heatmap && document.getElementById('seq-kinase-heatmap')) {
            createKinaseHeatmap('seq-kinase-heatmap', kinaseData.heatmap, 'sequence');
        }
    } catch (e) {
        console.error('Error initializing sequence kinase charts:', e);
    }
}

// Initialize combined kinase analysis
function initCombinedKinaseAnalysis() {
    const structDataScript = document.getElementById('struct-kinase-data');
    const seqDataScript = document.getElementById('seq-kinase-data');
    if (!structDataScript || !seqDataScript) return;
    
    try {
        // Parse kinase data
        const structData = JSON.parse(structDataScript.textContent);
        const seqData = JSON.parse(seqDataScript.textContent);
        
        if (structData.top_kinases && seqData.top_kinases) {
            // Create comparison chart
            createKinaseComparisonChart(
                'kinase-comparison-chart', 
                structData.top_kinases, 
                seqData.top_kinases
            );
            
            // Fill combined kinase table
            fillCombinedKinaseTable(
                'combined-kinase-table',
                structData.top_kinases, 
                seqData.top_kinases
            );
        }
    } catch (e) {
        console.error('Error initializing combined kinase analysis:', e);
    }
}

// Create bar chart for top kinases
function createKinaseBarChart(canvasId, kinaseData, dataType) {
    const canvas = document.getElementById(canvasId);
    if (!canvas || !kinaseData || kinaseData.length === 0) return;
    
    // Prepare chart data
    const labels = kinaseData.map(k => k.kinase);
    const scores = kinaseData.map(k => k.score);
    
    // Define colors based on score
    const colors = scores.map(score => {
        if (score >= 0.8) return 'rgba(40, 167, 69, 0.7)';  // Green for high confidence
        if (score >= 0.5) return 'rgba(255, 193, 7, 0.7)';  // Yellow for medium confidence
        return 'rgba(220, 53, 69, 0.7)';  // Red for low confidence
    });
    
    // Create chart
    const ctx = canvas.getContext('2d');
    
    // Destroy existing chart if it exists
    if (canvas.chart) {
        canvas.chart.destroy();
    }
    
    canvas.chart = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: labels,
            datasets: [{
                label: `${dataType.charAt(0).toUpperCase() + dataType.slice(1)} Kinase Prediction Scores`,
                data: scores,
                backgroundColor: colors,
                borderColor: colors.map(c => c.replace('0.7', '1')),
                borderWidth: 1
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    max: 1,
                    title: {
                        display: true,
                        text: 'Prediction Score'
                    }
                },
                x: {
                    title: {
                        display: true,
                        text: 'Kinase'
                    }
                }
            },
            plugins: {
                tooltip: {
                    callbacks: {
                        label: function(context) {
                            const score = context.raw;
                            let confidence = 'Low Confidence';
                            if (score >= 0.8) confidence = 'High Confidence';
                            else if (score >= 0.5) confidence = 'Medium Confidence';
                            return [`Score: ${score.toFixed(3)}`, `Confidence: ${confidence}`];
                        }
                    }
                }
            }
        }
    });
}

// Create heatmap for kinase comparison
function createKinaseHeatmap(containerId, heatmapData, dataType) {
    const container = document.getElementById(containerId);
    if (!container || !heatmapData || !heatmapData.scores || heatmapData.scores.length === 0) return;
    
    // Prepare data for heatmap
    const sites = [...new Set(heatmapData.scores.map(d => d.site))];
    const kinases = [...new Set(heatmapData.scores.map(d => d.kinase))];
    
    // Format scores for heatmap
    const formattedData = [];
    heatmapData.scores.forEach(score => {
        formattedData.push(score);
    });
    
    // Clear container
    container.innerHTML = '';
    
    // Set dimensions
    const margin = { top: 30, right: 50, bottom: 120, left: 150 };
    const width = container.clientWidth - margin.left - margin.right;
    const height = 400 - margin.top - margin.bottom;
    
    // Create SVG
    const svg = d3.select(container)
        .append('svg')
        .attr('width', width + margin.left + margin.right)
        .attr('height', height + margin.top + margin.bottom)
        .append('g')
        .attr('transform', `translate(${margin.left},${margin.top})`);
    
    // Scales
    const x = d3.scaleBand()
        .range([0, width])
        .domain(kinases)
        .padding(0.05);
    
    const y = d3.scaleBand()
        .range([height, 0])
        .domain(sites)
        .padding(0.05);
    
    // Add X axis label
    svg.append('g')
        .style('font-size', 12)
        .attr('transform', `translate(0,${height})`)
        .call(d3.axisBottom(x))
        .selectAll('text')
        .attr('transform', 'translate(-10,0)rotate(-45)')
        .style('text-anchor', 'end')
        .style('font-size', '10px');
    
    // Add Y axis label
    svg.append('g')
        .style('font-size', 12)
        .call(d3.axisLeft(y))
        .selectAll('text')
        .style('font-size', '10px');
    
    // Build color scale
    const colorScale = d3.scaleSequential()
        .interpolator(d3.interpolateViridis)
        .domain([0, 1]);
    
    // Create tooltip
    const tooltip = d3.select(container)
        .append('div')
        .style('opacity', 0)
        .attr('class', 'tooltip')
        .style('background-color', 'white')
        .style('border', 'solid')
        .style('border-width', '1px')
        .style('border-radius', '5px')
        .style('padding', '10px')
        .style('position', 'absolute');
    
    // Add squares
    svg.selectAll('rect')
        .data(formattedData)
        .enter()
        .append('rect')
        .attr('x', d => x(d.kinase))
        .attr('y', d => y(d.site))
        .attr('width', x.bandwidth())
        .attr('height', y.bandwidth())
        .style('fill', d => colorScale(d.score))
        .style('stroke-width', 4)
        .style('stroke', 'none')
        .style('opacity', 0.8)
        .on('mouseover', function(event, d) {
            // Show tooltip
            tooltip.style('opacity', 1)
                .html(`Site: ${d.site}<br>Kinase: ${d.kinase}<br>Score: ${d.score.toFixed(3)}`)
                .style('left', (event.pageX + 10) + 'px')
                .style('top', (event.pageY + 10) + 'px');
            
            // Highlight cell
            d3.select(this)
                .style('stroke', 'black');
        })
        .on('mousemove', function(event) {
            tooltip.style('left', (event.pageX + 10) + 'px')
                .style('top', (event.pageY + 10) + 'px');
        })
        .on('mouseleave', function() {
            tooltip.style('opacity', 0);
            d3.select(this).style('stroke', 'none');
        });
    
    // Add title
    svg.append('text')
        .attr('x', 0)
        .attr('y', -10)
        .attr('text-anchor', 'left')
        .style('font-size', '14px')
        .style('fill', '#333')
        .text(`${dataType.charAt(0).toUpperCase() + dataType.slice(1)} Kinase Prediction Heatmap`);
    
    // Add color scale legend
    const legend = svg.append('g')
        .attr('transform', `translate(0,${height + 50})`);
    
    const legendWidth = width;
    const legendHeight = 20;
    
    // Create gradient for legend
    const defs = svg.append('defs');
    const gradient = defs.append('linearGradient')
        .attr('id', `${dataType}-gradient`)
        .attr('x1', '0%')
        .attr('x2', '100%')
        .attr('y1', '0%')
        .attr('y2', '0%');
    
    // Add color stops
    const stops = [0, 0.2, 0.4, 0.6, 0.8, 1];
    stops.forEach(stop => {
        gradient.append('stop')
            .attr('offset', `${stop * 100}%`)
            .attr('stop-color', colorScale(stop));
    });
    
    // Draw legend rectangle
    legend.append('rect')
        .attr('width', legendWidth)
        .attr('height', legendHeight)
        .style('fill', `url(#${dataType}-gradient)`);
    
    // Add legend axis
    const legendScale = d3.scaleLinear()
        .domain([0, 1])
        .range([0, legendWidth]);
    
    legend.append('g')
        .attr('transform', `translate(0,${legendHeight})`)
        .call(d3.axisBottom(legendScale)
            .tickFormat(d3.format('.1f'))
            .ticks(5));
    
    // Add legend title
    legend.append('text')
        .attr('x', legendWidth / 2)
        .attr('y', -5)
        .attr('text-anchor', 'middle')
        .style('font-size', '12px')
        .text('Prediction Score');
}

// Create a comparison chart for structural vs sequence predictions
function createKinaseComparisonChart(canvasId, structuralData, sequenceData) {
    const canvas = document.getElementById(canvasId);
    if (!canvas || !structuralData || !sequenceData) return;
    
    // Get all kinases from both datasets
    const allKinases = new Set();
    structuralData.forEach(k => allKinases.add(k.kinase));
    sequenceData.forEach(k => allKinases.add(k.kinase));
    
    // Prepare data for chart
    const kinases = Array.from(allKinases);
    const structScores = [];
    const seqScores = [];
    
    // Get scores for each kinase
    kinases.forEach(kinase => {
        const structKinase = structuralData.find(k => k.kinase === kinase);
        structScores.push(structKinase ? structKinase.score : 0);
        
        const seqKinase = sequenceData.find(k => k.kinase === kinase);
        seqScores.push(seqKinase ? seqKinase.score : 0);
    });
    
    // Create chart
    const ctx = canvas.getContext('2d');
    
    // Destroy existing chart if it exists
    if (canvas.chart) {
        canvas.chart.destroy();
    }
    
    canvas.chart = new Chart(ctx, {
        type: 'radar',
        data: {
            labels: kinases,
            datasets: [
                {
                    label: 'Structural Score',
                    data: structScores,
                    backgroundColor: 'rgba(54, 162, 235, 0.2)',
                    borderColor: 'rgb(54, 162, 235)',
                    pointBackgroundColor: 'rgb(54, 162, 235)',
                    pointBorderColor: '#fff',
                    pointHoverBackgroundColor: '#fff',
                    pointHoverBorderColor: 'rgb(54, 162, 235)'
                },
                {
                    label: 'Sequence Score',
                    data: seqScores,
                    backgroundColor: 'rgba(255, 99, 132, 0.2)',
                    borderColor: 'rgb(255, 99, 132)',
                    pointBackgroundColor: 'rgb(255, 99, 132)',
                    pointBorderColor: '#fff',
                    pointHoverBackgroundColor: '#fff',
                    pointHoverBorderColor: 'rgb(255, 99, 132)'
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            elements: {
                line: {
                    borderWidth: 3
                }
            },
            scales: {
                r: {
                    angleLines: {
                        display: true
                    },
                    suggestedMin: 0,
                    suggestedMax: 1
                }
            }
        }
    });
}

// Create a kinase family analysis chart
function createKinaseFamilyChart(canvasId, kinaseData) {
    const canvas = document.getElementById(canvasId);
    if (!canvas || !kinaseData || kinaseData.length === 0) return;
    
    // Kinase family mapping
    const kinaseFamilies = {
        'CDK': ['CDK1', 'CDK2', 'CDK4', 'CDK5', 'CDK6', 'CDK7', 'CDK8', 'CDK9'],
        'MAPK': ['ERK1', 'ERK2', 'p38', 'JNK1', 'JNK2', 'JNK3', 'MAPK'],
        'GSK': ['GSK3', 'GSK3A', 'GSK3B'],
        'CK': ['CK1', 'CK2', 'CSNK1', 'CSNK2'],
        'PKC': ['PKC', 'PKCALPHA', 'PKCBETA', 'PKCDELTA', 'PKCEPSILON', 'PKCGAMMA', 'PKCZETA'],
        'PKA': ['PKA', 'PKACA', 'PKACB', 'PKACG'],
        'AKT': ['AKT', 'AKT1', 'AKT2', 'AKT3'],
        'SRC': ['SRC', 'FYN', 'LCK', 'LYN', 'HCK', 'FGR', 'BLK', 'YES'],
        'CAMK': ['CAMK', 'CAMK1', 'CAMK2', 'CAMK4'],
        'ATM/ATR': ['ATM', 'ATR', 'DNAPK'],
        'PLK': ['PLK1', 'PLK2', 'PLK3', 'PLK4'],
        'AURORA': ['AURKA', 'AURKB', 'AURKC'],
        'Other': []
    };
    
    // Assign each kinase to a family
    const familyScores = {};
    kinaseData.forEach(kinase => {
        let assigned = false;
        
        // Find the family for this kinase
        for (const [family, members] of Object.entries(kinaseFamilies)) {
            if (members.some(member => kinase.kinase.toUpperCase().includes(member))) {
                if (!familyScores[family]) {
                    familyScores[family] = 0;
                }
                familyScores[family] += kinase.score;
                assigned = true;
                break;
            }
        }
        
        // If not assigned to any family, put in 'Other'
        if (!assigned) {
            if (!familyScores['Other']) {
                familyScores['Other'] = 0;
            }
            familyScores['Other'] += kinase.score;
        }
    });
    
    // Convert to array for chart
    const labels = Object.keys(familyScores);
    const scores = Object.values(familyScores);
    
    // Generate colors
    const colors = [
        'rgba(255, 99, 132, 0.7)',
        'rgba(54, 162, 235, 0.7)',
        'rgba(255, 206, 86, 0.7)',
        'rgba(75, 192, 192, 0.7)',
        'rgba(153, 102, 255, 0.7)',
        'rgba(255, 159, 64, 0.7)',
        'rgba(199, 199, 199, 0.7)',
        'rgba(83, 102, 255, 0.7)',
        'rgba(255, 99, 255, 0.7)',
        'rgba(99, 255, 132, 0.7)',
        'rgba(255, 77, 77, 0.7)',
        'rgba(128, 128, 255, 0.7)',
        'rgba(204, 204, 0, 0.7)'
    ];
    
    // Create chart
    const ctx = canvas.getContext('2d');
    
    // Destroy existing chart if it exists
    if (canvas.chart) {
        canvas.chart.destroy();
    }
    
    canvas.chart = new Chart(ctx, {
        type: 'pie',
        data: {
            labels: labels,
            datasets: [{
                data: scores,
                backgroundColor: colors.slice(0, labels.length),
                borderColor: colors.map(c => c.replace('0.7', '1')).slice(0, labels.length),
                borderWidth: 1
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                legend: {
                    position: 'right',
                    labels: {
                        font: {
                            size: 12
                        }
                    }
                },
                tooltip: {
                    callbacks: {
                        label: function(context) {
                            const label = context.label || '';
                            const value = context.raw || 0;
                            const total = context.dataset.data.reduce((a, b) => a + b, 0);
                            const percentage = total > 0 ? ((value / total) * 100).toFixed(1) : 0;
                            return `${label}: ${value.toFixed(2)} (${percentage}%)`;
                        }
                    }
                }
            }
        }
    });
}

// Create a motif analysis visualization
function createMotifAnalysis(containerId, kinaseData) {
    const container = document.getElementById(containerId);
    if (!container || !kinaseData || kinaseData.length === 0) return;
    
    // Kinase motif mapping
    const kinaseMotifs = {
        'CDK1': 'S/T-P-X-K/R',
        'CDK2': 'S/T-P-X-K/R',
        'CDK5': 'S/T-P-X-K/R',
        'ERK1': 'P-X-S/T-P',
        'ERK2': 'P-X-S/T-P',
        'MAPK': 'P-X-S/T-P',
        'p38': 'P-X-S/T-P',
        'JNK': 'P-X-S/T-P',
        'PKA': 'R-R-X-S/T',
        'PKG': 'R-R-X-S/T',
        'PKC': 'S/T-X-K/R',
        'CK1': 'pS-X-X-S/T',
        'CK2': 'S/T-X-X-E/D',
        'GSK3': 'S/T-X-X-X-pS',
        'PLK': 'D/E-X-S/T-Φ',
        'AKT': 'R-X-R-X-X-S/T',
        'CAMK2': 'R-X-X-S/T',
        'AURORA': 'R-X-S/T',
        'ATM': 'S/T-Q',
        'ATR': 'S/T-Q',
        'DNAPK': 'S/T-Q',
        'SRC': 'X-E/D-E/D-X-Y-Φ-X-X',
        'FYN': 'X-E/D-E/D-X-Y-Φ-X-X',
        'LCK': 'X-E/D-E/D-X-Y-Φ-X-X'
    };
    
    // Take top 5 kinases
    const topKinases = kinaseData.slice(0, 5);
    
    // Create HTML for motif analysis
    let html = '<div class="row">';
    
    topKinases.forEach((kinase, index) => {
        // Find matching motif by looking for matches in the kinase name
        let motif = 'Unknown';
        for (const [kinaseRegex, kinaseMotif] of Object.entries(kinaseMotifs)) {
            if (kinase.kinase.toUpperCase().includes(kinaseRegex.toUpperCase())) {
                motif = kinaseMotif;
                break;
            }
        }
        
        // Create a colored motif representation
        html += `
        <div class="col-lg-6 mb-3">
            <div class="card h-100">
                <div class="card-header">
                    <h6 class="mb-0">${kinase.kinase} (Score: ${kinase.score.toFixed(3)})</h6>
                </div>
                <div class="card-body">
                    <p><strong>Consensus Motif:</strong> ${motif}</p>
                    <div class="motif-visualization">
                        ${createMotifVisualization(motif)}
                    </div>
                    <p class="mt-3 small text-muted">
                        This kinase recognizes sites with the specified motif pattern.
                        S/T indicates serine or threonine, X is any amino acid, 
                        pS is phosphorylated serine, and Φ is a hydrophobic residue.
                    </p>
                </div>
            </div>
        </div>
        `;
        
        // Create two columns per row
        if (index % 2 === 1 && index < topKinases.length - 1) {
            html += '</div><div class="row">';
        }
    });
    
    html += '</div>';
    container.innerHTML = html;
}

// Helper function to create a visual representation of a motif
function createMotifVisualization(motif) {
    // Split the motif into individual positions
    const positions = motif.split('');
    
    // Create the visualization
    let html = '<div style="display: flex; margin-top: 10px;">';
    
    positions.forEach(pos => {
        let color = '#e0e0e0'; // Default gray
        let textColor = '#000000';
        
        // Color based on position type
        if (pos === 'S' || pos === 'T' || pos.includes('S/T')) {
            color = '#bbdefb'; // Light blue for S/T
        } else if (pos === 'P') {
            color = '#81c784'; // Light green for P
        } else if (pos === 'K' || pos === 'R' || pos.includes('K/R')) {
            color = '#c8e6c9'; // Green for basic
        } else if (pos === 'D' || pos === 'E' || pos.includes('D/E')) {
            color = '#ffcdd2'; // Light red for acidic
        } else if (pos === 'Φ') {
            color = '#ffecb3'; // Light yellow for hydrophobic
        } else if (pos.includes('pS') || pos.includes('pT')) {
            color = '#ff5722'; // Orange for phosphorylated
            textColor = '#ffffff';
        }
        
        html += `<div style="width: 30px; height: 30px; display: flex; align-items: center; justify-content: center; margin: 0 2px; background-color: ${color}; border-radius: 4px; color: ${textColor};">${pos}</div>`;
    });
    
    html += '</div>';
    return html;
}

// Function to fill the combined kinase table
function fillCombinedKinaseTable(tableId, structData, seqData) {
    const tableBody = document.getElementById(tableId);
    if (!tableBody || !structData || !seqData) return;
    
    // Get all unique kinases
    const allKinases = new Set();
    structData.forEach(k => allKinases.add(k.kinase));
    seqData.forEach(k => allKinases.add(k.kinase));
    
    // Create map for easy lookups
    const structMap = new Map(structData.map(k => [k.kinase, k.score]));
    const seqMap = new Map(seqData.map(k => [k.kinase, k.score]));
    
    // Calculate combined scores
    const combinedScores = [];
    allKinases.forEach(kinase => {
        const structScore = structMap.get(kinase) || 0;
        const seqScore = seqMap.get(kinase) || 0;
        const avgScore = (structScore + seqScore) / 2;
        
        combinedScores.push({
            kinase: kinase,
            structScore: structScore,
            seqScore: seqScore,
            avgScore: avgScore
        });
    });
    
    // Sort by average score (descending)
    combinedScores.sort((a, b) => b.avgScore - a.avgScore);
    
    // Create table rows
    let html = '';
    combinedScores.forEach(score => {
        html += `
        <tr>
            <td>${score.kinase}</td>
            <td>${score.structScore.toFixed(3)}</td>
            <td>${score.seqScore.toFixed(3)}</td>
            <td>
                <strong>${score.avgScore.toFixed(3)}</strong>
                <div class="progress">
                    <div class="progress-bar 
                        ${score.avgScore >= 0.8 ? 'bg-success' : 
                          score.avgScore >= 0.5 ? 'bg-warning' : 'bg-danger'}"
                        style="width: ${score.avgScore * 100}%">
                    </div>
                </div>
            </td>
        </tr>
        `;
    });
    
    tableBody.innerHTML = html;
}