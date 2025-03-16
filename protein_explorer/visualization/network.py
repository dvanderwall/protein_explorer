"""
Functions for visualizing protein interaction networks.
"""

import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import logging
from typing import Dict, List, Tuple, Optional, Set, Union

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def create_network_figure(network: nx.Graph, 
                        highlighted_nodes: Optional[List[str]] = None,
                        node_color_map: Optional[Dict[str, str]] = None,
                        node_size_map: Optional[Dict[str, float]] = None,
                        title: str = "Protein Interaction Network") -> go.Figure:
    """
    Create a Plotly figure for a protein interaction network.
    
    Args:
        network: NetworkX graph object
        highlighted_nodes: List of nodes to highlight
        node_color_map: Dictionary mapping node IDs to colors
        node_size_map: Dictionary mapping node IDs to sizes
        title: Plot title
        
    Returns:
        Plotly figure object
    """
    # Create positions using a layout algorithm
    try:
        pos = nx.spring_layout(network, seed=42)
    except:
        pos = nx.random_layout(network)
    
    # Default colors and sizes
    default_node_color = '#6175c1'  # Medium blue
    highlight_color = '#FF4500'     # Orange-red
    default_node_size = 10
    highlight_size = 15
    
    # Create node trace
    node_x = []
    node_y = []
    node_colors = []
    node_sizes = []
    node_texts = []
    
    for node in network.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        
        # Node text (for hover information)
        node_info = f"ID: {node}"
        if 'name' in network.nodes[node]:
            node_info += f"<br>Name: {network.nodes[node]['name']}"
        node_texts.append(node_info)
        
        # Node color
        if node_color_map and node in node_color_map:
            node_colors.append(node_color_map[node])
        elif highlighted_nodes and node in highlighted_nodes:
            node_colors.append(highlight_color)
        else:
            node_colors.append(default_node_color)
        
        # Node size
        if node_size_map and node in node_size_map:
            node_sizes.append(node_size_map[node])
        elif highlighted_nodes and node in highlighted_nodes:
            node_sizes.append(highlight_size)
        else:
            node_sizes.append(default_node_size)
    
    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode='markers',
        hoverinfo='text',
        text=node_texts,
        marker=dict(
            color=node_colors,
            size=node_sizes,
            line=dict(width=1, color='#888')
        )
    )
    
    # Create edge trace
    edge_x = []
    edge_y = []
    edge_texts = []
    
    for edge in network.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        
        # Edge text (for hover information)
        if 'confidence' in network.edges[edge]:
            confidence = network.edges[edge]['confidence']
            edge_text = f"Confidence: {confidence:.2f}"
            edge_texts.extend([edge_text, edge_text, None])
    
    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode='lines',
        hoverinfo='text' if edge_texts else 'none',
        text=edge_texts if edge_texts else None,
        line=dict(width=0.5, color='#888'),
        opacity=0.7
    )
    
    # Create the figure
    fig = go.Figure(
        data=[edge_trace, node_trace],
        layout=go.Layout(
            title=title,
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='rgba(247,247,247,1)'
        )
    )
    
    return fig

def visualize_network(network: nx.Graph,
                     highlight_nodes: Optional[List[str]] = None,
                     node_colors: Optional[Dict[str, str]] = None,
                     title: str = "Protein Interaction Network") -> str:
    """
    Create an interactive network visualization.
    
    Args:
        network: NetworkX graph object
        highlight_nodes: List of nodes to highlight
        node_colors: Dictionary mapping node IDs to colors
        title: Plot title
        
    Returns:
        HTML string with the interactive visualization
    """
    # Create the network figure
    fig = create_network_figure(
        network,
        highlighted_nodes=highlight_nodes,
        node_color_map=node_colors,
        title=title
    )
    
    # Convert to HTML
    html = fig.to_html(full_html=False, include_plotlyjs='cdn')
    
    return html

def visualize_path(network: nx.Graph,
                  path: List[str],
                  title: str = "Protein Interaction Path") -> str:
    """
    Visualize a path through the protein interaction network.
    
    Args:
        network: NetworkX graph object
        path: List of proteins in the path
        title: Plot title
        
    Returns:
        HTML string with the path visualization
    """
    if not path or len(path) < 2:
        logger.error("Path must contain at least two proteins")
        return "<p>Error: Path must contain at least two proteins</p>"
    
    # Create a subgraph with only the nodes and edges in the path
    path_nodes = set(path)
    path_edges = [(path[i], path[i+1]) for i in range(len(path)-1)]
    
    # Create color map for path nodes
    node_colors = {}
    node_colors[path[0]] = '#4CAF50'  # Start node (green)
    node_colors[path[-1]] = '#F44336'  # End node (red)
    for node in path[1:-1]:
        node_colors[node] = '#FFC107'  # Intermediate nodes (amber)
    
    # Create the subgraph
    subgraph = network.subgraph(path_nodes).copy()
    
    # Create the network figure
    fig = create_network_figure(
        subgraph,
        node_color_map=node_colors,
        title=title
    )
    
    # Add annotations for start and end nodes
    try:
        pos = nx.spring_layout(subgraph, seed=42)
    except:
        pos = nx.random_layout(subgraph)
    
    start_x, start_y = pos[path[0]]
    end_x, end_y = pos[path[-1]]
    
    fig.add_annotation(
        x=start_x,
        y=start_y,
        text="Start",
        showarrow=True,
        arrowhead=1,
        ax=0,
        ay=-40
    )
    
    fig.add_annotation(
        x=end_x,
        y=end_y,
        text="End",
        showarrow=True,
        arrowhead=1,
        ax=0,
        ay=-40
    )
    
    # Convert to HTML
    html = fig.to_html(full_html=False, include_plotlyjs='cdn')
    
    return html

def visualize_clusters(network: nx.Graph,
                      clusters: Dict[str, int],
                      title: str = "Protein Interaction Network Clusters") -> str:
    """
    Visualize protein clusters in the network.
    
    Args:
        network: NetworkX graph object
        clusters: Dictionary mapping protein IDs to cluster IDs
        title: Plot title
        
    Returns:
        HTML string with the cluster visualization
    """
    # Get unique cluster IDs
    unique_clusters = sorted(set(clusters.values()))
    num_clusters = len(unique_clusters)
    
    # Generate a color palette for clusters
    if num_clusters <= 10:
        colors = px.colors.qualitative.Set1[:num_clusters]
    else:
        colors = px.colors.qualitative.Alphabet[:num_clusters]
    
    # Map clusters to colors
    cluster_colors = {cluster_id: colors[i % len(colors)] 
                     for i, cluster_id in enumerate(unique_clusters)}
    
    # Create node color map
    node_colors = {node: cluster_colors[clusters[node]] 
                  for node in clusters if node in network}
    
    # Create the network figure
    fig = create_network_figure(
        network,
        node_color_map=node_colors,
        title=title
    )
    
    # Add a legend for clusters
    for i, cluster_id in enumerate(unique_clusters):
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(size=10, color=cluster_colors[cluster_id]),
            name=f"Cluster {cluster_id}",
            showlegend=True
        ))
    
    fig.update_layout(showlegend=True)
    
    # Convert to HTML
    html = fig.to_html(full_html=False, include_plotlyjs='cdn')
    
    return html