# Protein Explorer

A comprehensive web application for visualizing and analyzing protein phosphorylation data, with a focus on structural and sequence similarity analysis of phosphosites.

## Overview

Protein Explorer is a sophisticated platform designed to help researchers explore protein structures, analyze phosphorylation sites, and understand kinase-substrate relationships. The application integrates multiple data sources including AlphaFold protein structures, UniProt data, and custom phosphosite structural similarity databases to provide deep insights into protein phosphorylation mechanisms.

### Key Features

- **Protein Structure Visualization**: Interactive 3D visualization of AlphaFold protein structures
- **Phosphosite Analysis**: Identification and analysis of potential phosphorylation sites
- **Structural Similarity Networks**: Visualization of structurally similar phosphosites
- **Sequence Similarity Analysis**: Identification of sequence-similar phosphosites
- **Kinase Prediction**: Structure-based and sequence-based kinase predictions
- **Motif Conservation**: Analysis of conserved sequence patterns around phosphosites
- **Interactive Data Exploration**: Rich interactive visualizations and filtering capabilities

## Architecture

The application follows a classic web application architecture with:

### Core Components

1. **Flask Web Application** (`web_app/app.py`): The main server that handles HTTP requests and serves the application
2. **Analysis Modules** (`protein_explorer/analysis/`): Core analytical functionality for proteins and phosphorylation sites
3. **Data Modules** (`protein_explorer/data/`): Handles data loading, caching, and retrieval from external sources
4. **Visualization Modules** (`protein_explorer/visualization/`): Generates visualizations for protein structures and networks
5. **Templates** (`web_app/templates/`): HTML templates for the user interface
6. **Static Files** (`web_app/static/`): JavaScript, CSS, and other static assets

### Data Flow

1. User submits a query (protein identifier or phosphosite)
2. Flask app routes the request to the appropriate handler
3. Data is retrieved from external sources (UniProt, AlphaFold) or local databases
4. Analysis modules process the data
5. Results are visualized through templates and client-side JavaScript
6. Interactive visualizations are rendered in the browser

## Installation and Setup

### Prerequisites

- Python 3.8+
- Flask and required packages (see `requirements.txt`)
- Node.js and npm (for frontend development)

### Installation Steps

1. Clone the repository
   ```bash
   git clone https://github.com/yourusername/protein-explorer.git
   cd protein-explorer
   ```

2. Set up a virtual environment
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install Python dependencies
   ```bash
   pip install -r requirements.txt
   ```

4. Download the required data files
   - `Combined_Kinome_10A_Master_Filtered_2.feather`: Structural similarity database
   - `Sequence_Similarity_Edges.parquet`: Sequence similarity database
   - `PhosphositeSuppData.feather`: Phosphosite supplementary data
   - `Structure_Kinase_Scores.feather`: Structure-based kinase scores
   - `Sequence_Kinase_Scores.feather`: Sequence-based kinase scores

5. Place the data files in the project root directory

6. Start the Flask development server
   ```bash
   python web_app/app.py
   ```

7. Open your browser and navigate to `http://localhost:5000`

## Core Modules

### Analysis Modules

- **enhanced_table.py**: Generates enhanced HTML tables for phosphosite visualization
- **kinase_predictor.py**: Predicts kinases for phosphorylation sites based on similarity
- **networks.py**: Analyzes protein interaction networks using linear algebra
- **phospho.py**: Analyzes potential phosphorylation sites in proteins
- **phospho_analyzer.py**: Handles phosphosite structural analysis and comparisons
- **sequence_analyzer.py**: Analyzes sequence similarity between phosphorylation sites
- **structure.py**: Functions for analyzing protein structures using linear algebra

### Web Application (app.py)

The main application file (`web_app/app.py`) contains:

- Flask route definitions for all pages
- API endpoints for data retrieval
- Data preprocessing and integration logic
- Template rendering with appropriate context data

### Key Routes

- `/`: Home page
- `/search`: Protein search page
- `/protein/<identifier>`: Detailed protein information page
- `/site/<uniprot_id>/<site>`: Detailed phosphorylation site page
- `/phosphosite`: Phosphosite structural analysis page
- `/site-search`: Search page for specific phosphorylation sites
- `/analyze`: Tool for analyzing multiple proteins
- `/faq`: Frequently asked questions

### API Endpoints

- `/api/phosphosites/<uniprot_id>`: Get phosphorylation site data for a protein
- `/api/structure/<uniprot_id>`: Get structure information for a protein
- `/api/network/<uniprot_id>`: Get protein interaction network data
- `/api/sequence_matches/<site_id>`: Get sequence similarity matches for a phosphosite
- `/api/sequence_conservation/<site_id>`: Get sequence conservation analysis for a phosphosite
- `/api/kinases/<site_id>`: Get kinase prediction scores for a phosphosite
- `/api/kinases/compare`: Compare kinase predictions across multiple sites

## Frontend Components

### Templates

The application uses a series of HTML templates located in `web_app/templates/`:

- **index.html**: Home page with search functionality and feature overview
- **protein.html**: Detailed protein information page with structure visualization
- **site.html**: Detailed phosphorylation site analysis page
- **site_structural_section.html**: Tab-based structural analysis section
- **site_sequence_section.html**: Tab-based sequence analysis section
- **combined_kinase_tab.html**: Combined kinase analysis tab content
- **structural_network_script.html**: JavaScript for structural network visualization
- **sequence_network_script.html**: JavaScript for sequence network visualization

### JavaScript Components

- **kinase_prediction.js**: Visualizes kinase prediction results with charts
- **phosphosite-visualization.js**: Handles visualization of phosphosites in protein context

### Visualization Techniques

- **3D Protein Structures**: Uses NGL Viewer for 3D molecular visualization
- **Network Graphs**: Uses D3.js for interactive network visualizations
- **Charts and Plots**: Uses Chart.js for kinase prediction visualizations
- **Heatmaps**: D3.js-based heatmaps for multi-site comparisons
- **Interactive Tables**: Enhanced tables with filtering and sorting capabilities

## Data Sources and Integration

### External Data Sources

- **AlphaFold**: Provides protein structure predictions
- **UniProt**: Provides protein sequence and annotation data
- **PhosphositePlus**: Reference data for known phosphorylation sites and kinases

### Internal Databases

- **Structural Similarity Database**: Precomputed RMSD values between phosphosites
- **Sequence Similarity Database**: Precomputed sequence similarity scores
- **Kinase Prediction Scores**: Precomputed structure and sequence-based kinase scores

### Data Loading and Caching

The application implements efficient data loading and caching strategies:

1. **Preloading**: Critical data is preloaded at application startup
2. **Local Caching**: Downloaded structures and data are cached locally
3. **In-Memory Caching**: Frequently accessed data is cached in memory
4. **Progressive Loading**: Data is loaded progressively as needed

## Key Analysis Features

### Phosphosite Identification

Identifies potential phosphorylation sites (S, T, Y residues) in protein sequences and analyzes their structural contexts. For each phosphosite, the application calculates:

- Mean pLDDT score (confidence metric from AlphaFold)
- Number of nearby residues within 10Å
- Surface accessibility
- Sequence motif (-7 to +7 positions)
- Comparison with known phosphorylation sites

### Structural Similarity Analysis

Compares the 3D structure of phosphorylation sites to identify similar binding patterns:

- RMSD-based structural similarity calculations
- Interactive network visualization of similar sites
- Filtering by RMSD threshold
- Detailed comparison of structural features

### Sequence Similarity Analysis

Analyzes sequence patterns around phosphorylation sites:

- Sequence similarity scoring
- Motif conservation analysis
- Position-specific amino acid distributions
- N-terminal and C-terminal region analysis

### Kinase Prediction

Predicts potential kinases for phosphorylation sites based on:

1. **Structure-based prediction**: Using structural similarity to known kinase substrates
2. **Sequence-based prediction**: Using sequence patterns recognized by specific kinases
3. **Combined analysis**: Integration of both approaches for more robust predictions

## Contributing

Contributions to Protein Explorer are welcome! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

[MIT License](LICENSE)

## Acknowledgements

- AlphaFold team for providing protein structure predictions
- UniProt for comprehensive protein data
- NGL Viewer for molecular visualization
- D3.js and Chart.js for data visualization
- Bootstrap for UI components

