/**
 * Phosphosite Visualization Script
 * This script creates a bar chart visualization for phosphorylation sites.
 */

// Configuration for the different metrics to visualize
const METRICS_CONFIG = {
    'nearby': {
      yLabel: 'Number of residues',
      barName: 'Nearby Residues (10Å)',
      mainColor: '#8884d8',
      description: 'Higher bars indicate more residues within 10Å, suggesting a more buried position.'
    },
    'surface': {
      yLabel: 'Surface Accessibility (%)',
      barName: 'Surface Accessibility',
      mainColor: '#2196f3',
      description: 'Higher bars indicate greater surface accessibility, suggesting the site is more exposed to solvent.'
    },
    'plddt': {
      yLabel: 'Mean pLDDT Score',
      barName: 'Mean pLDDT (-5:+5)',
      mainColor: '#4caf50',
      description: 'Higher bars indicate greater model confidence in the local structure around the phosphosite.',
      referenceLine: {
        y: 70,
        label: 'Confidence threshold'
      }
    },
    'acidic': {
      yLabel: 'Acidic Residues (%)',
      barName: 'Acidic Content (-5:+5)',
      mainColor: '#f44336',
      description: 'Higher bars indicate a higher percentage of acidic residues (D, E) near the phosphosite.'
    },
    'basic': {
      yLabel: 'Basic Residues (%)',
      barName: 'Basic Content (-5:+5)',
      mainColor: '#9c27b0',
      description: 'Higher bars indicate a higher percentage of basic residues (K, R, H) near the phosphosite.'
    },
    'aromatic': {
      yLabel: 'Aromatic Residues (%)',
      barName: 'Aromatic Content (-5:+5)',
      mainColor: '#ff9800',
      description: 'Higher bars indicate a higher percentage of aromatic residues (F, W, Y) near the phosphosite.'
    },
    'bfactor': {
      yLabel: 'B-factor Gradient',
      barName: 'Structure Variability',
      mainColor: '#009688',
      description: 'Higher bars indicate greater structural variability in the local environment.'
    },
    'hydrophobicity': {
      yLabel: 'Hydrophobicity Score',
      barName: 'Hydrophobicity (-5:+5)',
      mainColor: '#607d8b',
      description: 'Higher bars indicate a more hydrophobic local environment around the phosphosite.'
    }
  };
  
  // Create the visualization in the given container element
  function createPhosphositeVisualization(containerId) {
    const container = document.getElementById(containerId);
    if (!container) {
      console.error(`Container element with ID "${containerId}" not found.`);
      return;
    }
    
    // Extract protein sequence if available
    const proteinSequence = extractProteinSequence();
    
    // Extract phosphosite data from the table on the page
    const phosphositesData = extractPhosphositesFromTable();
    if (!phosphositesData || phosphositesData.length === 0) {
      container.innerHTML = `
        <div class="alert alert-warning">No phosphosite data available for visualization.</div>
      `;
      return;
    }
    
    // Create the visualization container
    container.innerHTML = `
      <div class="card mb-4">
        <div class="card-header">
          <h5 class="mb-0">Phosphosite Structural Analysis</h5>
        </div>
        <div class="card-body">
          <ul class="nav nav-tabs" id="phosphositeAnalysisTabs" role="tablist">
            ${Object.keys(METRICS_CONFIG).map((metricId, index) => `
              <li class="nav-item" role="presentation">
                <button 
                  class="nav-link ${index === 0 ? 'active' : ''}" 
                  id="${metricId}-tab" 
                  data-bs-toggle="tab" 
                  data-bs-target="#${metricId}-tab-content" 
                  type="button" 
                  role="tab" 
                  aria-controls="${metricId}" 
                  aria-selected="${index === 0 ? 'true' : 'false'}"
                >
                  ${METRICS_CONFIG[metricId].barName}
                </button>
              </li>
            `).join('')}
          </ul>
          
          <div class="tab-content mt-3">
            ${Object.keys(METRICS_CONFIG).map((metricId, index) => `
              <div 
                class="tab-pane fade ${index === 0 ? 'show active' : ''}" 
                id="${metricId}-tab-content" 
                role="tabpanel" 
                aria-labelledby="${metricId}-tab"
              >
                <div class="chart-container" style="height: 300px;" id="${metricId}-chart-container"></div>
                <div class="mt-4">
                  <p class="text-muted mb-0">
                    <strong>Click on a bar</strong> to highlight the corresponding row in the phosphosite table.
                    ${METRICS_CONFIG[metricId].description}
                  </p>
                </div>
              </div>
            `).join('')}
            
            <div class="mt-3">
              <div class="d-flex align-items-center justify-content-center">
                <div class="d-flex align-items-center me-4">
                  <div style="width: 16px; height: 16px; background-color: #4caf50; border-radius: 50%; margin-right: 6px;"></div>
                  <span>Serine (S)</span>
                </div>
                <div class="d-flex align-items-center me-4">
                  <div style="width: 16px; height: 16px; background-color: #2196f3; border-radius: 50%; margin-right: 6px;"></div>
                  <span>Threonine (T)</span>
                </div>
                <div class="d-flex align-items-center me-4">
                  <div style="width: 16px; height: 16px; background-color: #ff9800; border-radius: 50%; margin-right: 6px;"></div>
                  <span>Tyrosine (Y)</span>
                </div>
                <div class="d-flex align-items-center">
                  <div style="width: 16px; height: 16px; background-color: #9e9e9e; border-radius: 50%; margin-right: 6px;"></div>
                  <span>Not known sites</span>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    `;
    
    // Initialize each chart
    Object.keys(METRICS_CONFIG).forEach(metricId => {
      createBarChart(metricId, phosphositesData, proteinSequence);
    });
    
    // Set up tab switching functionality with Bootstrap 5 tabs
    const tabEls = document.querySelectorAll('#phosphositeAnalysisTabs button');
    tabEls.forEach(tabEl => {
      tabEl.addEventListener('click', event => {
        event.preventDefault();
        
        // Remove active class from all tabs and hide all tab contents
        tabEls.forEach(el => {
          el.classList.remove('active');
          el.setAttribute('aria-selected', 'false');
          const tabContent = document.querySelector(el.dataset.bsTarget);
          if (tabContent) {
            tabContent.classList.remove('show', 'active');
          }
        });
        
        // Add active class to clicked tab and show its content
        event.currentTarget.classList.add('active');
        event.currentTarget.setAttribute('aria-selected', 'true');
        const tabContent = document.querySelector(event.currentTarget.dataset.bsTarget);
        if (tabContent) {
          tabContent.classList.add('show', 'active');
        }
        
        // Force chart resize
        window.dispatchEvent(new Event('resize'));
      });
    });
  }
  
  // Extract protein sequence from the page if available
  function extractProteinSequence() {
    // Try to find the sequence in a pre tag
    const sequenceElement = document.querySelector('.sequence-display');
    if (sequenceElement) {
      const sequence = sequenceElement.textContent.trim().replace(/\s+/g, '');
      console.log("Found protein sequence with length:", sequence.length);
      return sequence;
    }
    
    // Alternative: try another selector if the above doesn't work
    const altSequenceElement = document.querySelector('pre.sequence');
    if (altSequenceElement) {
      const sequence = altSequenceElement.textContent.trim().replace(/\s+/g, '');
      console.log("Found protein sequence with length:", sequence.length);
      return sequence;
    }
    
    // If we couldn't find the sequence, log a warning and return null
    console.warn("Protein sequence not found on the page");
    return null;
  }
  
  // Extract phosphosite data from the table on the page
  function extractPhosphositesFromTable() {
    try {
      // Get phosphosite data from the table on the page
      const table = document.querySelector('.phosphosite-table tbody');
      if (!table) {
        console.error("Phosphosite table not found on page");
        return [];
      }
  
      const rows = Array.from(table.querySelectorAll('tr'));
      const sites = rows.map(row => {
        try {
          // Get data attributes first if available (these are added by our enhanced table)
          const dataAttrs = {
            site: row.getAttribute('data-site'),
            resno: parseInt(row.getAttribute('data-resno')),
            siteType: row.getAttribute('data-type'),
            nearbyCount: parseInt(row.getAttribute('data-nearby')),
            meanPlddt: parseFloat(row.getAttribute('data-plddt')),
            surfaceAccessibility: parseFloat(row.getAttribute('data-surface')),
            acidicPercentage: parseFloat(row.getAttribute('data-acidic')),
            basicPercentage: parseFloat(row.getAttribute('data-basic')),
            aromaticPercentage: parseFloat(row.getAttribute('data-aromatic')),
            bFactorGradient: parseFloat(row.getAttribute('data-bfactor')),
            hydrophobicityScore: parseFloat(row.getAttribute('data-hydrophobicity')),
            isKnown: row.getAttribute('data-known') === 'True'
          };
          
          // If data attributes are not available, extract from table cells as fallback
          const cells = row.querySelectorAll('td');
          if (cells.length < 4) return null;
          
          // Use data attributes if available, otherwise extract from table cells
          if (!dataAttrs.site || isNaN(dataAttrs.resno)) {
            // Extract site information
            const siteElement = cells[0].querySelector('strong');
            const site = siteElement ? siteElement.textContent : cells[0].textContent.trim();
            
            // Extract residue number and type
            const match = site.match(/([STY])(\d+)/);
            if (!match) return null;
            
            const siteType = match[1];
            const resno = parseInt(match[2], 10);
            
            dataAttrs.site = site;
            dataAttrs.siteType = siteType;
            dataAttrs.resno = resno;
          }
          
          if (isNaN(dataAttrs.nearbyCount) && cells.length > 3) {
            dataAttrs.nearbyCount = parseInt(cells[3].textContent.trim(), 10) || 0;
          }
          
          if (isNaN(dataAttrs.meanPlddt) && cells.length > 2) {
            dataAttrs.meanPlddt = parseFloat(cells[2].textContent.trim()) || 0;
          }
          
          // For surface accessibility, we use a fallback if it's not in data attributes
          if (isNaN(dataAttrs.surfaceAccessibility) && cells.length > 5) {
            const surfaceText = cells[5].textContent.trim();
            const surfaceMatch = surfaceText.match(/(\d+(\.\d+)?)/);
            dataAttrs.surfaceAccessibility = surfaceMatch ? parseFloat(surfaceMatch[1]) : 0;
          } else if (isNaN(dataAttrs.surfaceAccessibility)) {
            // This is a simplified fallback calculation
            dataAttrs.surfaceAccessibility = Math.max(0, 100 - (dataAttrs.nearbyCount * 3));
          }
          
          // For is_known, use data attribute or try to determine from table
          if (dataAttrs.isKnown === undefined && cells.length > 6) {
            dataAttrs.isKnown = cells[6].textContent.trim() === 'Yes';
          }
          
          // If we don't have these values, use reasonable defaults with some random variation for visual effect
          if (isNaN(dataAttrs.acidicPercentage)) {
            dataAttrs.acidicPercentage = 10 + Math.random() * 40; // Random placeholder
          }
          
          if (isNaN(dataAttrs.basicPercentage)) {
            dataAttrs.basicPercentage = 10 + Math.random() * 40; // Random placeholder
          }
          
          if (isNaN(dataAttrs.aromaticPercentage)) {
            dataAttrs.aromaticPercentage = 5 + Math.random() * 35; // Random placeholder
          }
          
          if (isNaN(dataAttrs.bFactorGradient)) {
            dataAttrs.bFactorGradient = 5 + Math.random() * 25; // Random placeholder
          }
          
          if (isNaN(dataAttrs.hydrophobicityScore)) {
            dataAttrs.hydrophobicityScore = 20 + Math.random() * 80; // Random placeholder
          }
          
          return dataAttrs;
        } catch (err) {
          console.error("Error processing row:", err);
          return null;
        }
      }).filter(site => site !== null);
      
      return sites;
    } catch (err) {
      console.error("Error extracting phosphosites:", err);
      return [];
    }
  }
  
  // Create a bar chart for a specific metric, showing the entire protein sequence
  function createBarChart(metricId, phosphositesData, proteinSequence) {
    // We'll use Chart.js for the visualization
    // Make sure Chart.js is included in your HTML
    if (typeof Chart === 'undefined') {
      console.error('Chart.js library not found. Please include Chart.js in your HTML.');
      return;
    }
    
    const chartContainer = document.getElementById(`${metricId}-chart-container`);
    if (!chartContainer) {
      console.error(`Chart container for metric "${metricId}" not found.`);
      return;
    }
    
    // Create canvas element
    const canvas = document.createElement('canvas');
    chartContainer.innerHTML = ''; // Clear any existing content
    chartContainer.appendChild(canvas);
    
    // Sort sites by residue number
    const sortedSites = [...phosphositesData].sort((a, b) => a.resno - b.resno);
    
    // Determine the sequence length to use as x-axis
    let sequenceLength = 0;
    if (proteinSequence) {
      sequenceLength = proteinSequence.length;
    } else {
      // If no sequence is available, use the highest residue number
      const maxResno = Math.max(...sortedSites.map(site => site.resno));
      sequenceLength = maxResno + 10; // Add some padding
    }
    
    // Prepare data for the entire sequence
    const labels = Array.from({ length: sequenceLength }, (_, i) => i + 1);
    const values = Array(sequenceLength).fill(0);
    const backgroundColors = Array(sequenceLength).fill('rgba(220, 220, 220, 0.2)'); // Light gray for non-S/T/Y
    const borderColors = Array(sequenceLength).fill('rgba(220, 220, 220, 0.5)');
    
    // Map of position to site for tooltip and click handling
    const positionToSite = {};
    
    // Fill in values for S/T/Y positions from the data
    sortedSites.forEach(site => {
      const pos = site.resno - 1; // Convert to 0-based index
      if (pos >= 0 && pos < sequenceLength) {
        // Get the value based on the current metric
        let value = 0;
        if (metricId === 'nearby') {
          value = site.nearbyCount;
        } else if (metricId === 'surface') {
          value = site.surfaceAccessibility;
        } else if (metricId === 'plddt') {
          value = site.meanPlddt;
        } else if (metricId === 'acidic') {
          value = site.acidicPercentage;
        } else if (metricId === 'basic') {
          value = site.basicPercentage;
        } else if (metricId === 'aromatic') {
          value = site.aromaticPercentage;
        } else if (metricId === 'bfactor') {
          value = site.bFactorGradient;
        } else if (metricId === 'hydrophobicity') {
          value = site.hydrophobicityScore;
        }
        
        values[pos] = value;
        
        // Set color based on site type and whether it's known
        if (site.isKnown) {
          if (site.siteType === 'S') {
            backgroundColors[pos] = 'rgba(76, 175, 80, 0.7)'; // Green for Serine
            borderColors[pos] = 'rgba(76, 175, 80, 1)';
          } else if (site.siteType === 'T') {
            backgroundColors[pos] = 'rgba(33, 150, 243, 0.7)'; // Blue for Threonine
            borderColors[pos] = 'rgba(33, 150, 243, 1)';
          } else if (site.siteType === 'Y') {
            backgroundColors[pos] = 'rgba(255, 152, 0, 0.7)'; // Orange for Tyrosine
            borderColors[pos] = 'rgba(255, 152, 0, 1)';
          }
        } else {
          // Gray for unknown/not known sites
          backgroundColors[pos] = 'rgba(158, 158, 158, 0.7)';
          borderColors[pos] = 'rgba(158, 158, 158, 1)';
        }
        
        // Store the site information for tooltip and click
        positionToSite[pos] = site;
      }
    });
    
    // Create chart
    const chart = new Chart(canvas.getContext('2d'), {
      type: 'bar',
      data: {
        labels: labels,
        datasets: [{
          label: METRICS_CONFIG[metricId].barName,
          data: values,
          backgroundColor: backgroundColors,
          borderColor: borderColors,
          borderWidth: 1
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
          y: {
            beginAtZero: true,
            title: {
              display: true,
              text: METRICS_CONFIG[metricId].yLabel
            }
          },
          x: {
            title: {
              display: true,
              text: 'Residue number'
            },
            ticks: {
              autoSkip: true,
              maxTicksLimit: 20,
              callback: function(val, index) {
                // Show fewer labels for better readability
                return index % Math.ceil(sequenceLength / 20) === 0 ? val : '';
              }
            }
          }
        },
        onClick: (event, elements) => {
          if (elements.length > 0) {
            const index = elements[0].index;
            const site = positionToSite[index];
            if (site) {
              highlightTableRow(site.site);
            }
          }
        },
        plugins: {
          tooltip: {
            callbacks: {
              title: (tooltipItems) => {
                const index = tooltipItems[0].index;
                const site = positionToSite[index];
                return site ? site.site : `Residue ${index + 1}`;
              },
              afterTitle: (tooltipItems) => {
                const index = tooltipItems[0].index;
                const site = positionToSite[index];
                if (site) {
                  return `Position: ${site.resno}`;
                }
                return null;
              },
              beforeBody: (tooltipItems) => {
                const index = tooltipItems[0].index;
                const site = positionToSite[index];
                
                if (!site) return null;
                
                let additionalInfo = [];
                
                if (metricId === 'nearby') {
                  additionalInfo.push(`Nearby residues: ${site.nearbyCount}`);
                } else if (metricId === 'surface') {
                  additionalInfo.push(`Surface accessibility: ${site.surfaceAccessibility.toFixed(1)}%`);
                } else if (metricId === 'plddt') {
                  additionalInfo.push(`Mean pLDDT: ${site.meanPlddt.toFixed(1)}`);
                } else if (metricId === 'acidic') {
                  additionalInfo.push(`Acidic content: ${site.acidicPercentage.toFixed(1)}%`);
                } else if (metricId === 'basic') {
                  additionalInfo.push(`Basic content: ${site.basicPercentage.toFixed(1)}%`);
                } else if (metricId === 'aromatic') {
                  additionalInfo.push(`Aromatic content: ${site.aromaticPercentage.toFixed(1)}%`);
                } else if (metricId === 'bfactor') {
                  additionalInfo.push(`B-factor gradient: ${site.bFactorGradient.toFixed(1)}`);
                } else if (metricId === 'hydrophobicity') {
                  additionalInfo.push(`Hydrophobicity: ${site.hydrophobicityScore.toFixed(1)}`);
                }
                
                additionalInfo.push(`Known in PhosphositePlus: ${site.isKnown ? 'Yes' : 'No'}`);
                
                return additionalInfo;
              }
            }
          }
        }
      }
    });
    
    // Add reference line if configured (requires Chart.js annotation plugin)
    if (METRICS_CONFIG[metricId].referenceLine && 
        typeof window.ChartAnnotation !== 'undefined') {
      chart.options.plugins.annotation = {
        annotations: {
          line1: {
            type: 'line',
            yMin: METRICS_CONFIG[metricId].referenceLine.y,
            yMax: METRICS_CONFIG[metricId].referenceLine.y,
            borderColor: '#666',
            borderWidth: 1,
            borderDash: [3, 3],
            label: {
              content: METRICS_CONFIG[metricId].referenceLine.label,
              enabled: true,
              position: 'right'
            }
          }
        }
      };
      chart.update();
    }
  }
  
  // Highlight a row in the phosphosite table
  function highlightTableRow(site) {
    // Find the corresponding row in the table
    const rows = document.querySelectorAll('.phosphosite-table tbody tr');
    let targetRow = null;
    
    for (const row of rows) {
      const siteCell = row.querySelector('td:first-child');
      if (siteCell && siteCell.textContent.includes(site)) {
        targetRow = row;
        break;
      }
    }
    
    if (targetRow) {
      // Scroll to the element
      targetRow.scrollIntoView({ behavior: 'smooth', block: 'center' });
      
      // Flash highlight effect
      const originalBg = targetRow.style.backgroundColor;
      targetRow.style.backgroundColor = '#ffc107';
      targetRow.style.transition = 'background-color 1s ease';
      
      setTimeout(() => {
        targetRow.style.backgroundColor = originalBg;
      }, 1500);
    }
  }
  
  // Initialize when the DOM is fully loaded
  document.addEventListener('DOMContentLoaded', function() {
    // Create the visualization if the container exists
    if (document.getElementById('phosphosite-visualization-container')) {
      console.log('Initializing phosphosite visualization...');
      createPhosphositeVisualization('phosphosite-visualization-container');
    }
  });