// Network-Based Kinase Prediction Visualization
document.addEventListener('DOMContentLoaded', function() {
    // Initialize network-based kinase charts
    initNetworkKinaseCharts();
    
    // Add tab change event listeners to ensure charts are properly sized
    document.querySelectorAll('.network-kinase-tab-btn').forEach(function(button) {
        button.addEventListener('click', function() {
            setTimeout(function() {
                window.dispatchEvent(new Event('resize'));
            }, 100);
        });
    });
    
    // Initialize RMSD threshold slider
    initThresholdSliders();
});

// Initialize network-based kinase charts
function initNetworkKinaseCharts() {
    const structDataScript = document.getElementById('struct-network-kinase-data');
    const seqDataScript = document.getElementById('seq-network-kinase-data');
    
    // Initialize structural charts if data exists
    if (structDataScript) {
        try {
            const kinaseData = JSON.parse(structDataScript.textContent);
            
            // Create visualizations for structure network kinase data
            createNetworkKinaseBarChart('struct-network-kinase-chart', kinaseData.top_kinases);
            createNetworkKinaseFamilyChart('struct-network-kinase-family-chart', kinaseData.top_kinases);
            // Add explicit 'structure' parameter
            createNetworkKinaseHeatmap('struct-network-kinase-heatmap', kinaseData.heatmap, 'structure');
            
        } catch (e) {
            console.error('Error initializing structural network kinase charts:', e);
        }
    }
    
    // Initialize sequence charts if data exists
    if (seqDataScript) {
        try {
            const kinaseData = JSON.parse(seqDataScript.textContent);
            
            // Create visualizations for sequence network kinase data
            createNetworkKinaseBarChart('seq-network-kinase-chart', kinaseData.top_kinases);
            createNetworkKinaseFamilyChart('seq-network-kinase-family-chart', kinaseData.top_kinases);
            // Add explicit 'sequence' parameter
            createNetworkKinaseHeatmap('seq-network-kinase-heatmap', kinaseData.heatmap, 'sequence');
            
        } catch (e) {
            console.error('Error initializing sequence network kinase charts:', e);
        }
    }
    
    // Initialize combined comparison if both data sources exist
    if (structDataScript && seqDataScript) {
        try {
            const structData = JSON.parse(structDataScript.textContent);
            const seqData = JSON.parse(seqDataScript.textContent);
            
            createNetworkKinaseComparisonChart(
                'network-kinase-comparison-chart',
                structData.top_kinases,
                seqData.top_kinases
            );
            
            fillNetworkKinaseTable(
                'network-kinase-table',
                structData.top_kinases,
                seqData.top_kinases
            );
            
        } catch (e) {
            console.error('Error initializing combined network kinase analysis:', e);
        }
    }
}

// Initialize threshold sliders and update buttons
function initThresholdSliders() {
    // RMSD threshold slider and update button for structural analysis
    const rmsdSlider = document.getElementById('rmsd-threshold');
    const rmsdValue = document.getElementById('rmsd-threshold-value');
    const updateRmsdBtn = document.getElementById('update-rmsd-threshold');
    
    if (rmsdSlider && rmsdValue) {
        rmsdSlider.addEventListener('input', function() {
            const threshold = parseFloat(this.value);
            rmsdValue.textContent = threshold.toFixed(1) + ' Å';
        });
    }
    
    if (updateRmsdBtn) {
        updateRmsdBtn.addEventListener('click', function() {
            const threshold = parseFloat(rmsdSlider.value);
            updateNetworkAnalysis('structure', threshold);
        });
    }
    
    // Sequence similarity slider and update button for sequence analysis
    const seqSlider = document.getElementById('seq-similarity-threshold');
    const seqValue = document.getElementById('seq-similarity-value');
    const updateSeqBtn = document.getElementById('update-seq-similarity-threshold');
    
    if (seqSlider && seqValue) {
        seqSlider.addEventListener('input', function() {
            const threshold = parseFloat(this.value);
            seqValue.textContent = threshold.toFixed(2);
        });
    }
    
    if (updateSeqBtn) {
        updateSeqBtn.addEventListener('click', function() {
            const threshold = parseFloat(seqSlider.value);
            updateNetworkAnalysis('sequence', threshold);
        });
    }
}

// Update network analysis based on new threshold
function updateNetworkAnalysis(analysisType, threshold) {
    // Show loading spinner
    const siteId = document.getElementById('site-id')?.value;
    if (!siteId) {
        console.error('Site ID not found');
        return;
    }
    
    // Create loading overlay
    const overlayId = analysisType === 'structure' ? 'struct-network-kinase' : 'seq-network-kinase';
    const tabContent = document.getElementById(overlayId);
    
    if (tabContent) {
        // Create overlay
        const loadingOverlay = document.createElement('div');
        loadingOverlay.id = 'loading-overlay';
        loadingOverlay.style.position = 'absolute';
        loadingOverlay.style.top = '0';
        loadingOverlay.style.left = '0';
        loadingOverlay.style.width = '100%';
        loadingOverlay.style.height = '100%';
        loadingOverlay.style.backgroundColor = 'rgba(255, 255, 255, 0.7)';
        loadingOverlay.style.display = 'flex';
        loadingOverlay.style.justifyContent = 'center';
        loadingOverlay.style.alignItems = 'center';
        loadingOverlay.style.zIndex = '1000';
        
        // Create spinner
        loadingOverlay.innerHTML = `
            <div class="spinner-border text-primary" role="status">
                <span class="visually-hidden">Loading...</span>
            </div>
        `;
        
        // Position relative for absolute overlay
        const originalPosition = tabContent.style.position;
        tabContent.style.position = 'relative';
        
        // Add overlay
        tabContent.appendChild(loadingOverlay);
        
        // Prepare request data
        const requestData = {
            score_type: analysisType
        };
        
        if (analysisType === 'structure') {
            requestData.rmsd_threshold = threshold;
        } else {
            requestData.similarity_threshold = threshold;
        }
        
        // Make AJAX request to update visualization
        fetch(`/api/update-network-kinases/${siteId}`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(requestData)
        })
        .then(response => {
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            return response.json();
        })
        .then(data => {
            
            // Update visualizations with new data
            const chartId = analysisType === 'structure' ? 'struct-network-kinase-chart' : 'seq-network-kinase-chart';
            const familyChartId = analysisType === 'structure' ? 'struct-network-kinase-family-chart' : 'seq-network-kinase-family-chart';
            const heatmapId = analysisType === 'structure' ? 'struct-network-kinase-heatmap' : 'seq-network-kinase-heatmap';
            
            createNetworkKinaseBarChart(chartId, data.top_kinases);
            createNetworkKinaseFamilyChart(familyChartId, data.top_kinases);
            createNetworkKinaseHeatmap(heatmapId, data.heatmap);
            
            // Update statistics table
            const tableBody = tabContent.querySelector('tbody');
            if (tableBody) {
                let tableHtml = '';
                data.top_kinases.forEach(kinase => {
                    tableHtml += `
                    <tr>
                        <td>${kinase.kinase}</td>
                        <td>${kinase.mean_score.toFixed(3)}</td>
                        <td>${kinase.median_score.toFixed(3)}</td>
                        <td>${kinase.min_score.toFixed(3)}</td>
                        <td>${kinase.max_score.toFixed(3)}</td>
                        <td>${kinase.sample_size}</td>
                        <td>${kinase.variability.toFixed(3)}</td>
                    </tr>
                    `;
                });
                tableBody.innerHTML = tableHtml;
            }
            
            // Update site count in info alert
            const infoAlert = tabContent.querySelector('.alert-info');
            if (infoAlert) {
                const countSpan = infoAlert.querySelector('strong');
                if (countSpan) {
                    countSpan.textContent = data.site_count;
                }
            }
            
            // Remove loading overlay
            tabContent.removeChild(loadingOverlay);
            tabContent.style.position = originalPosition;
        })
        .catch(error => {
            console.error('Error updating network analysis:', error);
            
            // Remove loading overlay
            tabContent.removeChild(loadingOverlay);
            tabContent.style.position = originalPosition;
            
        });
    }
}

// Create bar chart for network-based top kinases with error bars
function createNetworkKinaseBarChart(canvasId, kinaseData) {
    const canvas = document.getElementById(canvasId);
    if (!canvas || !kinaseData || kinaseData.length === 0) return;
    
    // Prepare chart data
    const labels = kinaseData.map(k => k.kinase);
    const meanScores = kinaseData.map(k => k.mean_score);
    const medianScores = kinaseData.map(k => k.median_score);
    const minScores = kinaseData.map(k => k.min_score);
    const maxScores = kinaseData.map(k => k.max_score);
    const sampleSizes = kinaseData.map(k => k.sample_size);
    
    // Define colors based on mean score
    const colors = meanScores.map(score => {
        if (score >= 0.8) return 'rgba(40, 167, 69, 0.7)';  // Green for high confidence
        if (score >= 0.5) return 'rgba(255, 193, 7, 0.7)';  // Yellow for medium confidence
        return 'rgba(220, 53, 69, 0.7)';  // Red for low confidence
    });
    
    // Load Chart.js annotation plugin if not already loaded
    if (!Chart.annotation && typeof chartjsPluginAnnotation !== 'undefined') {
        Chart.register(chartjsPluginAnnotation);
    }
    
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
            datasets: [
                {
                    label: 'Mean Score',
                    data: meanScores,
                    backgroundColor: colors,
                    borderColor: colors.map(c => c.replace('0.7', '1')),
                    borderWidth: 1
                },
                {
                    label: 'Median Score',
                    data: medianScores,
                    backgroundColor: 'rgba(0, 0, 0, 0)',
                    borderColor: 'rgba(0, 0, 0, 0.7)',
                    borderWidth: 2,
                    type: 'line',
                    pointStyle: 'rectRot',
                    pointRadius: 5,
                    pointHoverRadius: 7
                }
            ]
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
                        title: function(context) {
                            return `${context[0].label} Prediction`;
                        },
                        label: function(context) {
                            const index = context.dataIndex;
                            const dataset = context.dataset;
                            
                            if (dataset.label === 'Mean Score') {
                                return [
                                    `Mean Score: ${meanScores[index].toFixed(3)}`,
                                    `Median Score: ${medianScores[index].toFixed(3)}`,
                                    `Range: ${minScores[index].toFixed(3)} - ${maxScores[index].toFixed(3)}`,
                                    `Sample Size: ${sampleSizes[index]} sites`
                                ];
                            } else {
                                return `Median Score: ${medianScores[index].toFixed(3)}`;
                            }
                        }
                    }
                }
            }
        }
    });
    
    // Add error bars manually if annotation plugin isn't available
    try {
        // Add error bars
        for (let i = 0; i < kinaseData.length; i++) {
            const minY = ctx.height - (minScores[i] * ctx.height);
            const maxY = ctx.height - (maxScores[i] * ctx.height);
            const x = ((ctx.width - 40) / kinaseData.length) * (i + 0.5) + 20;
            
            // Draw vertical line
            ctx.beginPath();
            ctx.moveTo(x, minY);
            ctx.lineTo(x, maxY);
            ctx.strokeStyle = 'black';
            ctx.lineWidth = 1;
            ctx.stroke();
            
            // Draw horizontal caps
            ctx.beginPath();
            ctx.moveTo(x - 5, minY);
            ctx.lineTo(x + 5, minY);
            ctx.moveTo(x - 5, maxY);
            ctx.lineTo(x + 5, maxY);
            ctx.stroke();
        }
    } catch (e) {
        console.error('Error adding error bars manually:', e);
    }
}

// Enhanced createNetworkKinaseHeatmap function with better error handling and debugging
function createNetworkKinaseHeatmap(containerId, heatmapData, dataType) {
    const container = document.getElementById(containerId);
    
    if (!container) {
        return;
    }
    
    if (!heatmapData) {
        container.innerHTML = '<div class="alert alert-info">No heatmap data available.</div>';
        return;
    }
    
    // Verify scores array exists
    if (!heatmapData.scores) {
        container.innerHTML = '<div class="alert alert-info">Heatmap data format is invalid - missing scores array.</div>';
        return;
    }
    
    if (heatmapData.scores.length === 0) {
        container.innerHTML = '<div class="alert alert-info">No scores available for visualization.</div>';
        return;
    }
    
    // Clear container
    container.innerHTML = '';
    
    // Log some data to help with debugging
    
    try {
        // Prepare data for heatmap
        const sites = [...new Set(heatmapData.scores.map(d => d.site))];
        const kinases = [...new Set(heatmapData.scores.map(d => d.kinase))];
        
        
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
            .style('font-size', '10px')
            .text(d => {
                // Shorten site IDs for better display
                const parts = d.split('_');
                if (parts.length === 2) {
                    const uniprot = parts[0];
                    const site = parts[1];
                    return `${uniprot.substring(0, 6)}_${site}`;
                }
                return d;
            });
        
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
            .style('position', 'absolute')
            .style('z-index', '999');
        
        // Add squares
        svg.selectAll('rect')
            .data(heatmapData.scores)
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
            .text(`${dataType.charAt(0).toUpperCase() + dataType.slice(1)} Network Kinase Prediction Heatmap`);
        
        // Add color scale legend
        const legend = svg.append('g')
            .attr('transform', `translate(0,${height + 50})`);
        
        const legendWidth = width;
        const legendHeight = 20;
        
        // Create gradient for legend
        const defs = svg.append('defs');
        const gradientId = `${containerId}-gradient`;
        const gradient = defs.append('linearGradient')
            .attr('id', gradientId)
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
            .style('fill', `url(#${gradientId})`);
        
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
        
        
    } catch (error) {
        console.error(`Error creating heatmap for ${containerId}:`, error);
        container.innerHTML = `<div class="alert alert-danger">
            Error creating heatmap: ${error.message}<br>
            See console for details.
        </div>`;
    }
}


// Create a family distribution chart for network-based predictions
function createNetworkKinaseFamilyChart(canvasId, kinaseData) {
    // Use categorize_kinases_by_family logic from the server
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
    
    const canvas = document.getElementById(canvasId);
    if (!canvas || !kinaseData || kinaseData.length === 0) return;
    
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
                familyScores[family] += kinase.mean_score;
                assigned = true;
                break;
            }
        }
        
        // If not assigned to any family, put in 'Other'
        if (!assigned) {
            if (!familyScores['Other']) {
                familyScores['Other'] = 0;
            }
            familyScores['Other'] += kinase.mean_score;
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
        'rgba(128, 128, 255, 0.7)'
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

// Create a comparison chart for network-based structural vs sequence predictions
function createNetworkKinaseComparisonChart(canvasId, structuralData, sequenceData) {
    const canvas = document.getElementById(canvasId);
    if (!canvas || !structuralData || !sequenceData) return;
    
    // Get all kinases from both datasets
    const allKinases = new Set();
    structuralData.forEach(k => allKinases.add(k.kinase));
    sequenceData.forEach(k => allKinases.add(k.kinase));
    
    // Prepare data for chart
    const kinases = Array.from(allKinases);
    const structScores = [];
    const structVariability = [];
    const seqScores = [];
    const seqVariability = [];
    
    // Get scores for each kinase
    kinases.forEach(kinase => {
        const structKinase = structuralData.find(k => k.kinase === kinase);
        structScores.push(structKinase ? structKinase.mean_score : 0);
        structVariability.push(structKinase ? structKinase.variability : 0);
        
        const seqKinase = sequenceData.find(k => k.kinase === kinase);
        seqScores.push(seqKinase ? seqKinase.mean_score : 0);
        seqVariability.push(seqKinase ? seqKinase.variability : 0);
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
                    label: 'Structural Network Score',
                    data: structScores,
                    backgroundColor: 'rgba(54, 162, 235, 0.2)',
                    borderColor: 'rgb(54, 162, 235)',
                    pointBackgroundColor: 'rgb(54, 162, 235)',
                    pointBorderColor: '#fff',
                    pointHoverBackgroundColor: '#fff',
                    pointHoverBorderColor: 'rgb(54, 162, 235)'
                },
                {
                    label: 'Sequence Network Score',
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
            },
            plugins: {
                tooltip: {
                    callbacks: {
                        label: function(context) {
                            const index = context.dataIndex;
                            const datasetIndex = context.datasetIndex;
                            const dataset = context.dataset;
                            
                            // Get the score and variability
                            const score = dataset.data[index].toFixed(3);
                            const variability = (datasetIndex === 0 ? structVariability : seqVariability)[index].toFixed(3);
                            
                            return [
                                `${dataset.label}: ${score}`,
                                `Variability: ${variability}`
                            ];
                        }
                    }
                }
            }
        }
    });
}

// Fill in the network kinase table for the combined analysis
function fillNetworkKinaseTable(tableId, structData, seqData) {
    const tableBody = document.getElementById(tableId);
    if (!tableBody || !structData || !seqData) return;
    
    // Get all unique kinases
    const allKinases = new Set();
    structData.forEach(k => allKinases.add(k.kinase));
    seqData.forEach(k => allKinases.add(k.kinase));
    
    // Create map for easy lookups
    const structMap = new Map(structData.map(k => [k.kinase, k]));
    const seqMap = new Map(seqData.map(k => [k.kinase, k]));
    
    // Calculate combined scores
    const combinedScores = [];
    allKinases.forEach(kinase => {
        const structData = structMap.get(kinase);
        const seqData = seqMap.get(kinase);
        
        const structScore = structData ? structData.mean_score : 0;
        const seqScore = seqData ? seqData.mean_score : 0;
        const avgScore = (structScore + seqScore) / 2;
        
        const structSampleSize = structData ? structData.sample_size : 0;
        const seqSampleSize = seqData ? seqData.sample_size : 0;
        
        combinedScores.push({
            kinase: kinase,
            structScore: structScore,
            seqScore: seqScore,
            avgScore: avgScore,
            structSampleSize: structSampleSize,
            seqSampleSize: seqSampleSize
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
            <td>${score.structScore.toFixed(3)} <span class="badge bg-secondary">${score.structSampleSize}</span></td>
            <td>${score.seqScore.toFixed(3)} <span class="badge bg-secondary">${score.seqSampleSize}</span></td>
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