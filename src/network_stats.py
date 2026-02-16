from __future__ import annotations

import os
from typing import Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

from src.config import NODES_PATH, EDGES_PATH, FIGURES_DIR


def _ensure_dirs() -> None:
    os.makedirs(FIGURES_DIR, exist_ok=True)


def load_nodes_edges() -> tuple[pd.DataFrame, pd.DataFrame]:
    nodes = pd.read_csv(NODES_PATH)
    edges = pd.read_csv(EDGES_PATH)
    return nodes, edges


def node_counts(nodes: pd.DataFrame) -> pd.DataFrame:
    out = nodes["node_type"].value_counts().reset_index()
    out.columns = ["node_type", "count"]
    return out


def edge_counts(edges: pd.DataFrame) -> pd.DataFrame:
    out = edges["relation"].value_counts().reset_index()
    out.columns = ["relation", "count"]
    return out


def build_nx_graph(edges: pd.DataFrame) -> nx.Graph:
    # For basic structure stats we use an undirected projection
    # Poster-friendly: components, degree distribution, density
    g = nx.Graph()
    g.add_edges_from(edges[["source", "target"]].itertuples(index=False, name=None))
    return g


def compute_graph_stats(g: nx.Graph) -> Dict[str, float]:
    n = g.number_of_nodes()
    m = g.number_of_edges()

    density = nx.density(g) if n > 1 else 0.0
    num_components = nx.number_connected_components(g) if n > 0 else 0

    largest_cc = 0
    if n > 0 and num_components > 0:
        largest_cc = max(len(c) for c in nx.connected_components(g))

    return {
        "num_nodes": float(n),
        "num_edges": float(m),
        "density": float(density),
        "num_components": float(num_components),
        "largest_component_nodes": float(largest_cc),
        "largest_component_frac": float(largest_cc / n) if n > 0 else 0.0,
    }


def plot_degree_distribution(g: nx.Graph, out_png: str) -> None:
    degrees = np.array([d for _, d in g.degree()])
    degrees = degrees[degrees >= 0]

    plt.figure()
    plt.hist(degrees, bins=50)
    plt.xlabel("Degree")
    plt.ylabel("Count")
    plt.title("Degree distribution (undirected projection)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def top_hubs(g: nx.Graph, top_k: int = 20) -> pd.DataFrame:
    deg = sorted(g.degree(), key=lambda x: x[1], reverse=True)[:top_k]
    return pd.DataFrame(deg, columns=["node_id", "degree"])


def main() -> None:
    _ensure_dirs()

    nodes, edges = load_nodes_edges()

    # Quick correctness checks (poster and sanity)
    kinase_nodes = (nodes["node_id"].astype(str).str.startswith("KINASE:")).sum()
    phosph_nodes = (nodes["node_id"].astype(str).str.startswith("PHOSPHATASE:")).sum()

    print("=== Basic node and edge counts ===")
    print(node_counts(nodes).to_string(index=False))
    print()
    print(edge_counts(edges).to_string(index=False))
    print()

    print("=== Identity sanity checks ===")
    print("KINASE nodes:", int(kinase_nodes))
    print("PHOSPHATASE nodes:", int(phosph_nodes))
    print()

    print("=== Building NetworkX graph (undirected projection) ===")
    g = build_nx_graph(edges)

    stats = compute_graph_stats(g)
    for k, v in stats.items():
        if "num_" in k or "largest_" in k:
            print(f"{k}: {int(v)}" if float(v).is_integer() else f"{k}: {v:.4f}")
        else:
            print(f"{k}: {v:.6f}")
    print()

    hubs = top_hubs(g, top_k=20)
    print("=== Top 20 hubs by degree (undirected) ===")
    print(hubs.to_string(index=False))
    print()

    # Save artifacts for poster
    stats_df = pd.DataFrame([stats])
    stats_df.to_csv(os.path.join(FIGURES_DIR, "network_summary_stats.csv"), index=False)
    node_counts(nodes).to_csv(os.path.join(FIGURES_DIR, "node_counts.csv"), index=False)
    edge_counts(edges).to_csv(os.path.join(FIGURES_DIR, "edge_counts.csv"), index=False)
    hubs.to_csv(os.path.join(FIGURES_DIR, "top20_hubs.csv"), index=False)

    plot_degree_distribution(g, os.path.join(FIGURES_DIR, "degree_hist.png"))

    print("Wrote poster-ready outputs to /figures:")
    print("  network_summary_stats.csv")
    print("  node_counts.csv")
    print("  edge_counts.csv")
    print("  top20_hubs.csv")
    print("  degree_hist.png")


if __name__ == "__main__":
    main()
