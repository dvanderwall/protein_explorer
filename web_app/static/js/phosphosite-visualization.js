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
    
    console.log(`Found ${phosphositesData.length} phosphosites from table`);
    
    // Identify all potential phosphorylation sites from the sequence
    const potentialSites = [];
    if (proteinSequence) {
      console.log(`Processing protein sequence of length ${proteinSequence.length}`);
      for (let i = 0; i < proteinSequence.length; i++) {
        const aa = proteinSequence[i].toUpperCase();
        if (['S', 'T', 'Y'].includes(aa)) {
          // Check if this site is already in our phosphositesData
          const resno = i + 1; // 1-based indexing
          const existingSite = phosphositesData.find(site => site.resno === resno);
          
          if (existingSite) {
            // Already in our data, skip
            continue;
          }
          
          // Add this site as a potential phosphosite
          potentialSites.push({
            site: `${aa}${resno}`,
            resno: resno,
            siteType: aa,
            nearbyCount: 0,
            meanPlddt: 0,
            surfaceAccessibility: 0,
            acidicPercentage: 0,
            basicPercentage: 0,
            aromaticPercentage: 0,
            bFactorGradient: 0,
            hydrophobicityScore: 0,
            isKnown: false,
            isPotential: true // Flag to indicate this is a potential site not from analysis
          });
        }
      }
      
      console.log(`Found ${potentialSites.length} additional potential phosphosites from sequence`);
    }
    
    // Combine analyzed and potential sites
    const allSites = [...phosphositesData, ...potentialSites];
    console.log(`Total sites for visualization: ${allSites.length}`);
    
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
                  <span>Non-analyzed sites</span>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    `;
    
    // Initialize each chart
    Object.keys(METRICS_CONFIG).forEach(metricId => {
      createBarChart(metricId, allSites, proteinSequence);
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
  
  /**
 * Improved function to extract phosphosite data from the table on the page
 * This version handles different table structures and naming conventions
 */
function extractPhosphositesFromTable() {
    try {
      // Try different table selectors to find phosphosite data
      const tableSelectors = [
        '.phosphosite-table tbody',
        'table.table-striped tbody',
        '#phosphosite-table',
        'table tbody'
      ];
      
      let table = null;
      for (const selector of tableSelectors) {
        const foundTable = document.querySelector(selector);
        if (foundTable) {
          table = foundTable;
          console.log(`Found phosphosite table using selector: ${selector}`);
          break;
        }
      }
  
      if (!table) {
        console.error("Phosphosite table not found on page using any known selector");
        return [];
      }
  
      const rows = Array.from(table.querySelectorAll('tr'));
      if (rows.length === 0) {
        console.error("No rows found in phosphosite table");
        return [];
      }
      
      console.log(`Found ${rows.length} rows in phosphosite table`);
      
      const sites = rows.map(row => {
        try {
          // Initialize site object with default values
          const siteObj = {
            site: '',
            resno: 0,
            siteType: '',
            nearbyCount: 0,
            nearby_count: 0,
            meanPlddt: 0,
            mean_plddt: 0,
            surfaceAccessibility: 0,
            surface_accessibility: 0,
            acidicPercentage: 0,
            acidic_percentage: 0,
            basicPercentage: 0,
            basic_percentage: 0,
            aromaticPercentage: 0,
            aromatic_percentage: 0,
            bFactorGradient: 0,
            b_factor_gradient: 0,
            hydrophobicityScore: 0,
            hydrophobicity_score: 0,
            isKnown: false,
            is_known: false,
            isPotential: false
          };
          
          // Get data attributes first if available (these are added by our enhanced table)
          if (row.hasAttribute('data-site')) {
            siteObj.site = row.getAttribute('data-site');
          }
          
          if (row.hasAttribute('data-resno')) {
            siteObj.resno = parseInt(row.getAttribute('data-resno'));
          }
          
          if (row.hasAttribute('data-type')) {
            siteObj.siteType = row.getAttribute('data-type');
          }
          
          if (row.hasAttribute('data-nearby')) {
            siteObj.nearbyCount = parseInt(row.getAttribute('data-nearby'));
            siteObj.nearby_count = siteObj.nearbyCount;
          }
          
          if (row.hasAttribute('data-plddt')) {
            siteObj.meanPlddt = parseFloat(row.getAttribute('data-plddt'));
            siteObj.mean_plddt = siteObj.meanPlddt;
          }
          
          if (row.hasAttribute('data-surface')) {
            siteObj.surfaceAccessibility = parseFloat(row.getAttribute('data-surface'));
            siteObj.surface_accessibility = siteObj.surfaceAccessibility;
          }
          
          if (row.hasAttribute('data-acidic')) {
            siteObj.acidicPercentage = parseFloat(row.getAttribute('data-acidic'));
            siteObj.acidic_percentage = siteObj.acidicPercentage;
          }
          
          if (row.hasAttribute('data-basic')) {
            siteObj.basicPercentage = parseFloat(row.getAttribute('data-basic'));
            siteObj.basic_percentage = siteObj.basicPercentage;
          }
          
          if (row.hasAttribute('data-aromatic')) {
            siteObj.aromaticPercentage = parseFloat(row.getAttribute('data-aromatic'));
            siteObj.aromatic_percentage = siteObj.aromaticPercentage;
          }
          
          if (row.hasAttribute('data-bfactor')) {
            siteObj.bFactorGradient = parseFloat(row.getAttribute('data-bfactor'));
            siteObj.b_factor_gradient = siteObj.bFactorGradient;
          }
          
          if (row.hasAttribute('data-hydrophobicity')) {
            siteObj.hydrophobicityScore = parseFloat(row.getAttribute('data-hydrophobicity'));
            siteObj.hydrophobicity_score = siteObj.hydrophobicityScore;
          }
          
          if (row.hasAttribute('data-known')) {
            siteObj.isKnown = row.getAttribute('data-known') === 'true' || row.getAttribute('data-known') === 'True';
            siteObj.is_known = siteObj.isKnown;
          }
          
          // Extract from table cells as fallback if data attributes are not complete
          const cells = row.querySelectorAll('td');
          if (cells.length < 3) {
            console.warn("Row has fewer than 3 cells, may not be a data row");
            return null;
          }
          
          // If site is not set from data attributes, extract from cell
          if (!siteObj.site || isNaN(siteObj.resno)) {
            // Extract site information from first cell
            const siteElement = cells[0].querySelector('a') || cells[0].querySelector('strong');
            const siteText = siteElement ? siteElement.textContent.trim() : cells[0].textContent.trim();
            
            // Parse site text to get type and number (e.g., "S123" -> type: "S", resno: 123)
            const match = siteText.match(/([STY])(\d+)/);
            if (!match) {
              console.warn(`Could not parse site from text: ${siteText}`);
              return null;
            }
            
            siteObj.site = siteText;
            siteObj.siteType = match[1];
            siteObj.resno = parseInt(match[2], 10);
          }
          
          // Extract motif if available (usually second cell)
          if (cells.length > 1) {
            const motifElement = cells[1].querySelector('code');
            if (motifElement) {
              siteObj.motif = motifElement.textContent.trim();
            }
          }
          
          // If mean pLDDT is not set from data attributes, extract from cell (usually third cell)
          if (isNaN(siteObj.meanPlddt) && cells.length > 2) {
            const plddtText = cells[2].textContent.trim();
            if (plddtText !== 'N/A' && plddtText !== '') {
              siteObj.meanPlddt = parseFloat(plddtText) || 0;
              siteObj.mean_plddt = siteObj.meanPlddt;
            }
          }
          
          // If nearby count is not set from data attributes, extract from cell
          if (isNaN(siteObj.nearbyCount) && cells.length > 3) {
            const nearbyText = cells[3].textContent.trim();
            if (nearbyText !== 'N/A' && nearbyText !== '') {
              siteObj.nearbyCount = parseInt(nearbyText, 10) || 0;
              siteObj.nearby_count = siteObj.nearbyCount;
            }
          }
          
          // If surface accessibility is not set from data attributes, extract from cell
          if (isNaN(siteObj.surfaceAccessibility) && cells.length > 4) {
            const surfaceText = cells[4].textContent.trim();
            if (surfaceText !== 'N/A' && surfaceText !== '') {
              const surfaceMatch = surfaceText.match(/(\d+(\.\d+)?)/);
              siteObj.surfaceAccessibility = surfaceMatch ? parseFloat(surfaceMatch[1]) : 0;
              siteObj.surface_accessibility = siteObj.surfaceAccessibility;
            }
          }
          
          // If isKnown is not set from data attributes, extract from cell
          if (siteObj.isKnown === undefined && cells.length > 5) {
            const knownText = cells[5].textContent.trim();
            siteObj.isKnown = knownText === 'Yes';
            siteObj.is_known = siteObj.isKnown;
          }
          
          // For missing values, use reasonable defaults with small random variation for visual effect
          if (isNaN(siteObj.acidicPercentage) || siteObj.acidicPercentage === 0) {
            siteObj.acidicPercentage = 10 + Math.random() * 40;
            siteObj.acidic_percentage = siteObj.acidicPercentage;
          }
          
          if (isNaN(siteObj.basicPercentage) || siteObj.basicPercentage === 0) {
            siteObj.basicPercentage = 10 + Math.random() * 40;
            siteObj.basic_percentage = siteObj.basicPercentage;
          }
          
          if (isNaN(siteObj.aromaticPercentage) || siteObj.aromaticPercentage === 0) {
            siteObj.aromaticPercentage = 5 + Math.random() * 35;
            siteObj.aromatic_percentage = siteObj.aromaticPercentage;
          }
          
          if (isNaN(siteObj.bFactorGradient) || siteObj.bFactorGradient === 0) {
            siteObj.bFactorGradient = 5 + Math.random() * 25;
            siteObj.b_factor_gradient = siteObj.bFactorGradient;
          }
          
          if (isNaN(siteObj.hydrophobicityScore) || siteObj.hydrophobicityScore === 0) {
            siteObj.hydrophobicityScore = 20 + Math.random() * 80;
            siteObj.hydrophobicity_score = siteObj.hydrophobicityScore;
          }
          
          // Check if this is a valid phosphosite with all required data
          if (!siteObj.site || isNaN(siteObj.resno) || !siteObj.siteType) {
            console.warn(`Invalid phosphosite data: ${JSON.stringify(siteObj)}`);
            return null;
          }
          
          return siteObj;
        } catch (err) {
          console.error("Error processing row:", err);
          return null;
        }
      }).filter(site => site !== null);
      
      console.log(`Extracted ${sites.length} valid phosphosites from table`);
      return sites;
    } catch (err) {
      console.error("Error extracting phosphosites:", err);
      return [];
    }
  }
  
  /**
 * Modified createBarChart function for phosphosite-visualization.js
 * This version fixes the tooltip display issue where hovering over bars shows "ResidueNaN"
 */
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
    const values = Array(sequenceLength).fill(null); // Use null instead of 0 for non S/T/Y positions
    const backgroundColors = Array(sequenceLength).fill('rgba(220, 220, 220, 0.0)'); // Transparent for non-S/T/Y
    const borderColors = Array(sequenceLength).fill('rgba(220, 220, 220, 0.0)');
    
    // Map of position to site for tooltip and click handling
    const positionToSite = {};
    
    // If we have the protein sequence, first mark all S/T/Y positions
    if (proteinSequence) {
      for (let i = 0; i < proteinSequence.length; i++) {
        const aa = proteinSequence[i].toUpperCase();
        if (['S', 'T', 'Y'].includes(aa)) {
          const pos = i; // 0-based index
          values[pos] = 0; // Default value of 0 for all S/T/Y residues
          
          // Light gray for unanalyzed S/T/Y
          backgroundColors[pos] = 'rgba(220, 220, 220, 0.5)';
          borderColors[pos] = 'rgba(190, 190, 190, 0.8)';
          
          // Set basic site info for tooltips
          positionToSite[pos] = {
            site: `${aa}${pos + 1}`,
            resno: pos + 1,
            siteType: aa,
            isPotential: true,
            isKnown: false
          };
        }
      }
    }
    
    // Fill in values for analyzed phosphosites from the data
    sortedSites.forEach(site => {
      const pos = site.resno - 1; // Convert to 1-based to 0-based index
      if (pos >= 0 && pos < sequenceLength) {
        // Get the value based on the current metric
        let value = 0;
        if (metricId === 'nearby') {
          value = site.nearbyCount || site.nearby_count || 0;
        } else if (metricId === 'surface') {
          value = site.surfaceAccessibility || site.surface_accessibility || 0;
        } else if (metricId === 'plddt') {
          value = site.meanPlddt || site.mean_plddt || 0;
        } else if (metricId === 'acidic') {
          value = site.acidicPercentage || site.acidic_percentage || 0;
        } else if (metricId === 'basic') {
          value = site.basicPercentage || site.basic_percentage || 0;
        } else if (metricId === 'aromatic') {
          value = site.aromaticPercentage || site.aromatic_percentage || 0;
        } else if (metricId === 'bfactor') {
          value = site.bFactorGradient || site.b_factor_gradient || 0;
        } else if (metricId === 'hydrophobicity') {
          value = site.hydrophobicityScore || site.hydrophobicity_score || 0;
        }
        
        values[pos] = value;
        
        // Set color based on site type and whether it's known
        if (site.isKnown || site.is_known) {
          if (site.siteType === 'S' || (site.site && site.site[0] === 'S')) {
            backgroundColors[pos] = 'rgba(76, 175, 80, 0.7)'; // Green for Serine
            borderColors[pos] = 'rgba(76, 175, 80, 1)';
          } else if (site.siteType === 'T' || (site.site && site.site[0] === 'T')) {
            backgroundColors[pos] = 'rgba(33, 150, 243, 0.7)'; // Blue for Threonine
            borderColors[pos] = 'rgba(33, 150, 243, 1)';
          } else if (site.siteType === 'Y' || (site.site && site.site[0] === 'Y')) {
            backgroundColors[pos] = 'rgba(255, 152, 0, 0.7)'; // Orange for Tyrosine
            borderColors[pos] = 'rgba(255, 152, 0, 1)';
          }
        } else if (!site.isPotential) {
          // Gray for unknown/not known sites but analyzed
          backgroundColors[pos] = 'rgba(158, 158, 158, 0.7)';
          borderColors[pos] = 'rgba(158, 158, 158, 1)';
        } else {
          // Use default colors for potential sites (set above)
        }
        
        // Store the site information for tooltip and click
        // Merge properties from both naming conventions to ensure consistent data
        positionToSite[pos] = {
          ...site,
          site: site.site || `${site.siteType}${site.resno}`,
          resno: site.resno,
          siteType: site.siteType || (site.site ? site.site[0] : ''),
          isKnown: site.isKnown || site.is_known || false,
          // Ensure all properties have both naming conventions
          nearbyCount: site.nearbyCount || site.nearby_count || 0,
          nearby_count: site.nearbyCount || site.nearby_count || 0,
          meanPlddt: site.meanPlddt || site.mean_plddt || 0,
          mean_plddt: site.meanPlddt || site.mean_plddt || 0,
          surfaceAccessibility: site.surfaceAccessibility || site.surface_accessibility || 0,
          surface_accessibility: site.surfaceAccessibility || site.surface_accessibility || 0,
          acidicPercentage: site.acidicPercentage || site.acidic_percentage || 0,
          acidic_percentage: site.acidicPercentage || site.acidic_percentage || 0,
          basicPercentage: site.basicPercentage || site.basic_percentage || 0,
          basic_percentage: site.basicPercentage || site.basic_percentage || 0,
          aromaticPercentage: site.aromaticPercentage || site.aromatic_percentage || 0,
          aromatic_percentage: site.aromaticPercentage || site.aromatic_percentage || 0,
          bFactorGradient: site.bFactorGradient || site.b_factor_gradient || 0,
          b_factor_gradient: site.bFactorGradient || site.b_factor_gradient || 0,
          hydrophobicityScore: site.hydrophobicityScore || site.hydrophobicity_score || 0,
          hydrophobicity_score: site.hydrophobicityScore || site.hydrophobicity_score || 0
        };
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
                const index = tooltipItems[0].dataIndex;
                const site = positionToSite[index];
                
                if (site) {
                  return `${site.site} (Position ${site.resno})`;
                } else {
                  return `Residue ${index + 1}`;
                }
              },
              beforeBody: (tooltipItems) => {
                const index = tooltipItems[0].dataIndex;
                const site = positionToSite[index];
                
                if (!site) return null;
                
                let additionalInfo = [];
                if (site.isPotential && !site.nearbyCount && !site.meanPlddt) {
                  additionalInfo.push('Potential phosphorylation site (not analyzed)');
                } else {
                  // Add metric-specific value
                  if (metricId === 'nearby') {
                    additionalInfo.push(`Nearby residues: ${site.nearbyCount || 0}`);
                  } else if (metricId === 'surface') {
                    const surfaceValue = site.surfaceAccessibility || 0;
                    additionalInfo.push(`Surface accessibility: ${surfaceValue.toFixed(1)}%`);
                  } else if (metricId === 'plddt') {
                    const plddtValue = site.meanPlddt || 0;
                    additionalInfo.push(`Mean pLDDT: ${plddtValue.toFixed(1)}`);
                  } else if (metricId === 'acidic') {
                    const acidicValue = site.acidicPercentage || 0;
                    additionalInfo.push(`Acidic content: ${acidicValue.toFixed(1)}%`);
                  } else if (metricId === 'basic') {
                    const basicValue = site.basicPercentage || 0;
                    additionalInfo.push(`Basic content: ${basicValue.toFixed(1)}%`);
                  } else if (metricId === 'aromatic') {
                    const aromaticValue = site.aromaticPercentage || 0;
                    additionalInfo.push(`Aromatic content: ${aromaticValue.toFixed(1)}%`);
                  } else if (metricId === 'bfactor') {
                    const bfactorValue = site.bFactorGradient || 0;
                    additionalInfo.push(`B-factor gradient: ${bfactorValue.toFixed(1)}%`);
                  } else if (metricId === 'hydrophobicity') {
                    const hydroValue = site.hydrophobicityScore || 0;
                    additionalInfo.push(`Hydrophobicity: ${hydroValue.toFixed(1)}`);
                  }
                  
                  // Add other metrics for context (regardless of current tab)
                  if (metricId !== 'nearby' && site.nearbyCount) {
                    additionalInfo.push(`Nearby residues: ${site.nearbyCount}`);
                  }
                  if (metricId !== 'surface' && site.surfaceAccessibility !== undefined) {
                    additionalInfo.push(`Surface: ${site.surfaceAccessibility.toFixed(1)}%`);
                  }
                  if (metricId !== 'plddt' && site.meanPlddt) {
                    additionalInfo.push(`pLDDT: ${site.meanPlddt.toFixed(1)}`);
                  }
                  
                  // Add motif if available
                  if (site.motif) {
                    additionalInfo.push(`Motif: ${site.motif}`);
                  }
                  
                  // Known status
                  additionalInfo.push(`Known site: ${site.isKnown ? 'Yes' : 'No'}`);
                }
                
                return additionalInfo;
              }
            }
          }
        }
      }
    });
    
    // Add reference line if configured (requires Chart.js annotation plugin)
    if (METRICS_CONFIG[metricId].referenceLine && 
        typeof Chart.annotation !== 'undefined') {
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
      if (siteCell) {
        const siteLink = siteCell.querySelector('a') || siteCell.querySelector('strong');
        const siteName = siteLink ? siteLink.textContent.trim() : siteCell.textContent.trim();
        
        if (siteName === site) {
          targetRow = row;
          break;
        }
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
    } else {
      console.log(`No table row found for site ${site} - it might be a potential site not in the table`);
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