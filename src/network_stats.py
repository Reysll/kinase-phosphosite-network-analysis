# src/network_stats.py
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import networkx as nx

from src.config import NODES_PATH, EDGES_PATH, FIGURES_DIR


# ----------------------------
# Plot styling for posters
# ----------------------------
def _set_poster_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "font.size": 12,
            "axes.titlesize": 18,
            "axes.labelsize": 14,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 12,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linestyle": "-",
        }
    )


def _comma_axis(ax) -> None:
    ax.xaxis.set_major_formatter(mtick.StrMethodFormatter("{x:,.0f}"))
    ax.yaxis.set_major_formatter(mtick.StrMethodFormatter("{x:,.0f}"))


def _savefig(path: str) -> None:
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight")
    plt.close()


# ----------------------------
# Core loading
# ----------------------------
@dataclass
class DataBundle:
    nodes: pd.DataFrame
    edges: pd.DataFrame


def load_data(nodes_path: str, edges_path: str) -> DataBundle:
    nodes = pd.read_csv(nodes_path)
    edges = pd.read_csv(edges_path)
    return DataBundle(nodes=nodes, edges=edges)


# ----------------------------
# NetworkX graph builder
# ----------------------------
def build_undirected_projection(edges: pd.DataFrame) -> nx.Graph:
    """
    Poster analytics: treat everything as undirected for global structure stats.
    Deduplicate undirected edges to avoid double counting.
    """
    g = nx.Graph()
    # Add nodes later when we load nodes.csv.gz
    # Here we only add edges.
    for _, r in edges.iterrows():
        s = r["source"]
        t = r["target"]
        if s == t:
            continue
        g.add_edge(s, t)
    return g


# ----------------------------
# Counts tables + CSV exports
# ----------------------------
def compute_counts(nodes: pd.DataFrame, edges: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    node_counts = (
        nodes.groupby("node_type")["node_id"]
        .count()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )

    edge_counts = (
        edges.groupby("relation")["source"]
        .count()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )
    return node_counts, edge_counts


def write_csvs(node_counts: pd.DataFrame, edge_counts: pd.DataFrame, hubs: pd.DataFrame, summary: pd.DataFrame) -> None:
    os.makedirs(FIGURES_DIR, exist_ok=True)
    node_counts.to_csv(os.path.join(FIGURES_DIR, "node_counts.csv"), index=False)
    edge_counts.to_csv(os.path.join(FIGURES_DIR, "edge_counts.csv"), index=False)
    hubs.to_csv(os.path.join(FIGURES_DIR, "top20_hubs.csv"), index=False)
    summary.to_csv(os.path.join(FIGURES_DIR, "network_summary_stats.csv"), index=False)


# ----------------------------
# Poster plots
# ----------------------------
def plot_node_counts(node_counts: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.barh(node_counts["node_type"][::-1], node_counts["count"][::-1])
    ax.set_title("Node counts by type")
    ax.set_xlabel("Count")
    ax.set_ylabel("Node type")
    _comma_axis(ax)

    # labels
    for y, v in enumerate(node_counts["count"][::-1].to_list()):
        ax.text(v, y, f"  {v:,}", va="center", fontsize=12)

    _savefig(os.path.join(FIGURES_DIR, "node_type_counts.png"))


def plot_edge_counts(edge_counts: pd.DataFrame) -> None:
    # Full (linear) with comma formatting
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(edge_counts["relation"][::-1], edge_counts["count"][::-1])
    ax.set_title("Edge counts by relation (linear)")
    ax.set_xlabel("Count")
    ax.set_ylabel("Relation")
    _comma_axis(ax)

    # labels
    for y, v in enumerate(edge_counts["count"][::-1].to_list()):
        ax.text(v, y, f"  {v:,}", va="center", fontsize=11)

    _savefig(os.path.join(FIGURES_DIR, "edge_type_counts_linear.png"))

    # Full (log x) so small categories are visible
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(edge_counts["relation"][::-1], edge_counts["count"][::-1])
    ax.set_xscale("log")
    ax.set_title("Edge counts by relation (log scale)")
    ax.set_xlabel("Count (log scale)")
    ax.set_ylabel("Relation")
    ax.xaxis.set_major_formatter(mtick.LogFormatterSciNotation(base=10))

    _savefig(os.path.join(FIGURES_DIR, "edge_type_counts_log.png"))

    # Excluding dominant category (usually protein_same_pathway)
    dominant = "protein_same_pathway"
    subset = edge_counts[edge_counts["relation"] != dominant].copy()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(subset["relation"][::-1], subset["count"][::-1])
    ax.set_title(f"Edge counts by relation (excluding {dominant})")
    ax.set_xlabel("Count")
    ax.set_ylabel("Relation")
    _comma_axis(ax)

    for y, v in enumerate(subset["count"][::-1].to_list()):
        ax.text(v, y, f"  {v:,}", va="center", fontsize=11)

    _savefig(os.path.join(FIGURES_DIR, "edge_type_counts_excluding_dominant.png"))


def plot_degree_ccdf(degrees: np.ndarray) -> None:
    # CCDF P(D >= k)
    deg = np.asarray(degrees, dtype=int)
    deg = deg[deg > 0]
    if len(deg) == 0:
        return

    k_vals, counts = np.unique(deg, return_counts=True)
    # CCDF computed from tail cumulative
    tail = np.cumsum(counts[::-1])[::-1]
    ccdf = tail / tail[0]

    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.plot(k_vals, ccdf, linewidth=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Degree CCDF (log-log)")
    ax.set_xlabel("k (degree)")
    ax.set_ylabel("P(D ≥ k)")
    ax.grid(True, which="both", alpha=0.25)

    _savefig(os.path.join(FIGURES_DIR, "degree_ccdf_loglog.png"))


def plot_degree_log_binned_pdf(degrees: np.ndarray, bins: int = 30) -> None:
    """
    Log-binned degree distribution (PDF-like) looks much better than raw histogram
    for heavy-tailed networks.
    """
    deg = np.asarray(degrees, dtype=float)
    deg = deg[deg > 0]
    if len(deg) == 0:
        return

    dmin = max(1.0, deg.min())
    dmax = deg.max()

    edges = np.logspace(np.log10(dmin), np.log10(dmax), num=bins + 1)
    hist, bin_edges = np.histogram(deg, bins=edges, density=True)

    # Use geometric midpoints for x values
    mids = np.sqrt(bin_edges[:-1] * bin_edges[1:])

    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.plot(mids, hist, marker="o", linewidth=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Log-binned degree distribution (log-log)")
    ax.set_xlabel("Degree (log scale)")
    ax.set_ylabel("Density (log scale)")
    ax.grid(True, which="both", alpha=0.25)

    _savefig(os.path.join(FIGURES_DIR, "degree_logbinned_pdf_loglog.png"))


def plot_degree_hist_linear_p99(degrees: np.ndarray) -> None:
    deg = np.asarray(degrees, dtype=int)
    if len(deg) == 0:
        return

    p99 = np.percentile(deg, 99)
    capped = deg[deg <= p99]

    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.hist(capped, bins=60)
    ax.set_title("Degree distribution (linear, capped at 99th percentile)")
    ax.set_xlabel("Degree")
    ax.set_ylabel("Node count")
    _comma_axis(ax)

    _savefig(os.path.join(FIGURES_DIR, "degree_hist_linear_p99.png"))


# ----------------------------
# Main
# ----------------------------
def main() -> None:
    _set_poster_style()
    os.makedirs(FIGURES_DIR, exist_ok=True)

    bundle = load_data(NODES_PATH, EDGES_PATH)
    nodes = bundle.nodes
    edges = bundle.edges

    print("=== Basic node and edge counts ===")
    node_counts, edge_counts = compute_counts(nodes, edges)
    print(node_counts.to_string(index=False))
    print()
    print(edge_counts.to_string(index=False))
    print()

    print("=== Identity sanity checks ===")
    kinase_nodes = (nodes["node_id"].astype(str).str.startswith("KINASE:")).sum()
    ppase_nodes = (nodes["node_id"].astype(str).str.startswith("PHOSPHATASE:")).sum()
    print(f"KINASE nodes: {kinase_nodes}")
    print(f"PHOSPHATASE nodes: {ppase_nodes}")
    print()

    print("=== Building NetworkX graph (undirected projection) ===")
    g = build_undirected_projection(edges)

    # Ensure all nodes included (even isolated, though you currently have 1 component)
    for nid in nodes["node_id"].astype(str).tolist():
        if nid not in g:
            g.add_node(nid)

    num_nodes = g.number_of_nodes()
    num_edges = g.number_of_edges()
    density = nx.density(g)

    comps = list(nx.connected_components(g))
    num_components = len(comps)
    largest = max(comps, key=len) if comps else set()
    largest_n = len(largest)
    largest_frac = (largest_n / num_nodes) if num_nodes else 0.0

    print(f"num_nodes: {num_nodes}")
    print(f"num_edges: {num_edges}")
    print(f"density: {density:.6f}")
    print(f"num_components: {num_components}")
    print(f"largest_component_nodes: {largest_n}")
    print(f"largest_component_frac: {largest_frac:.3f}")
    print()

    print("=== Top 20 hubs by degree (undirected) ===")
    deg_series = pd.Series(dict(g.degree()), name="degree")
    hubs = (
        deg_series.sort_values(ascending=False)
        .head(20)
        .reset_index()
        .rename(columns={"index": "node_id"})
    )
    print(hubs.to_string(index=False))
    print()

    summary = pd.DataFrame(
        [
            {"metric": "num_nodes", "value": num_nodes},
            {"metric": "num_edges", "value": num_edges},
            {"metric": "density", "value": density},
            {"metric": "num_components", "value": num_components},
            {"metric": "largest_component_nodes", "value": largest_n},
            {"metric": "largest_component_frac", "value": largest_frac},
        ]
    )

    # CSV outputs
    write_csvs(node_counts, edge_counts, hubs, summary)

    # Plots
    plot_node_counts(node_counts)
    plot_edge_counts(edge_counts)

    degrees = np.array([d for _, d in g.degree()], dtype=int)
    plot_degree_hist_linear_p99(degrees)
    plot_degree_ccdf(degrees)
    plot_degree_log_binned_pdf(degrees)

    print(f"Wrote poster-ready outputs to: {FIGURES_DIR}")
    print("CSV files:")
    print("  node_counts.csv")
    print("  edge_counts.csv")
    print("  top20_hubs.csv")
    print("  network_summary_stats.csv")
    print("Plots:")
    print("  node_type_counts.png")
    print("  edge_type_counts_linear.png")
    print("  edge_type_counts_log.png")
    print("  edge_type_counts_excluding_dominant.png")
    print("  degree_hist_linear_p99.png")
    print("  degree_ccdf_loglog.png")
    print("  degree_logbinned_pdf_loglog.png")


if __name__ == "__main__":
    main()
