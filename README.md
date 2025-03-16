# ProteinExplorer

A Python package for exploring protein structures, interactions, and functional networks.

## Features

- Retrieve protein structures from AlphaFold DB by UniProt ID or gene symbol
- Navigate protein-protein interaction networks
- Analyze protein structures using PCA and other techniques
- Identify functional clusters in protein networks
- Visualize protein structures and interaction networks interactively

## Installation

Install from PyPI:

```bash
pip install protein-explorer
```

Or install directly from GitHub:

```bash
git clone https://github.com/class-account/protein-explorer.git
cd protein-explorer
pip install -e .
```

## Usage

### Basic Usage

```python
import protein_explorer as pe

# Get protein by UniProt ID
protein = pe.data.get_protein_by_id(uniprot_id="P53_HUMAN")

# Get protein by gene symbol
protein = pe.data.get_protein_by_id(gene_symbol="TP53")

# Get AlphaFold structure
structure = pe.data.get_alphafold_structure("P53_HUMAN")

# Visualize the structure
viewer_html = pe.visualization.visualize_structure(structure)
```

### Network Analysis

```python
# Build a protein interaction network
network = pe.navigation.build_interaction_network(["P53_HUMAN"], max_depth=2)

# Find path between two proteins
path = pe.navigation.find_path(network, "P53_HUMAN", "MDM2_HUMAN")

# Analyze the network
centrality = pe.analysis.compute_eigenvector_centrality(network)
clusters = pe.analysis.perform_spectral_clustering(network, n_clusters=5)

# Visualize the network
network_viz = pe.visualization.visualize_network(network, highlight_nodes=["P53_HUMAN"])
```

## Web Application

An interactive web application is available at: https://protein-explorer-app.herokuapp.com

The web application allows you to:
- Search for proteins by UniProt ID or gene symbol
- View protein structures in 3D using MolStar
- Explore protein-protein interaction networks
- Perform various analyses on selected proteins

## Documentation

Full documentation is available at: https://protein-explorer.readthedocs.io

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.