# src/subgraph_viz.py
from __future__ import annotations

import os
import math
import random
from collections import defaultdict

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from src.config import NODES_PATH, EDGES_PATH, FIGURES_DIR


# Relations present in your current pipeline
RELATION_ORDER = [
    "protein_same_pathway",
    "ppi_high_confidence",
    "site_coevolution",
    "site_same_pathway",
    "phosphorylates",
    "dephosphorylates",
    "has_site",
]

DEFAULT_HUBS = [
    "PROTEIN:CDK1",
    "PROTEIN:MAPK1",
    "PROTEIN:AKT1",
    "PROTEIN:EGFR",
    "PROTEIN:PRKCA",
]


def _safe_mkdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def load_nodes_edges() -> tuple[pd.DataFrame, pd.DataFrame]:
    nodes = pd.read_csv(NODES_PATH)
    edges = pd.read_csv(EDGES_PATH)

    # Ensure columns exist
    if "node_id" not in nodes.columns or "node_type" not in nodes.columns:
        raise ValueError("nodes file must have columns: node_id, node_type")

    for col in ["source", "target", "relation"]:
        if col not in edges.columns:
            raise ValueError(f"edges file must have column: {col}")

    # Some edges have blank confidence_score, keep it but fill missing
    if "confidence_score" not in edges.columns:
        edges["confidence_score"] = None

    return nodes, edges


def build_graph(edges: pd.DataFrame, undirected: bool = True) -> nx.Graph:
    if undirected:
        G = nx.Graph()
    else:
        G = nx.DiGraph()

    for row in edges.itertuples(index=False):
        s = getattr(row, "source")
        t = getattr(row, "target")
        rel = getattr(row, "relation")
        conf = getattr(row, "confidence_score", None)

        # Add edge with attributes
        G.add_edge(s, t, relation=rel, confidence_score=conf)

    return G


def ego_subgraph(
    G: nx.Graph,
    center: str,
    radius: int = 1,
    max_neighbors: int = 250,
    neighbor_priority_relations: list[str] | None = None,
) -> nx.Graph:
    """
    Build an ego subgraph around `center`, but cap the number of neighbors
    so the visualization stays readable.

    Strategy:
    - Collect all nodes within `radius` using ego_graph
    - If too many neighbors at radius=1, keep the top `max_neighbors` neighbors
      based on a score that prefers edges from selected relations.
    """
    if center not in G:
        raise ValueError(f"Center node not found in graph: {center}")

    H = nx.ego_graph(G, center, radius=radius, undirected=True)

    if radius != 1:
        return H

    # Cap neighbors for readability
    neighbors = list(H.neighbors(center))
    if len(neighbors) <= max_neighbors:
        return H

    rel_priority = neighbor_priority_relations or [
        "phosphorylates",
        "dephosphorylates",
        "has_site",
        "ppi_high_confidence",
        "protein_same_pathway",
        "site_same_pathway",
        "site_coevolution",
    ]
    rel_rank = {r: i for i, r in enumerate(rel_priority)}

    # Score neighbors by relation priority first, then by neighbor degree
    scores = []
    for nb in neighbors:
        edata = H.get_edge_data(center, nb) or {}
        rel = edata.get("relation", "")
        rank = rel_rank.get(rel, len(rel_rank))
        nb_deg = H.degree(nb)
        scores.append((rank, -nb_deg, nb))

    scores.sort()
    keep = set([center] + [nb for _, _, nb in scores[:max_neighbors]])

    # Induce a subgraph containing center + selected neighbors + edges among them
    return H.subgraph(keep).copy()


def filter_edges_by_relation(edges: pd.DataFrame, keep_relations: set[str]) -> pd.DataFrame:
    return edges[edges["relation"].isin(keep_relations)].copy()


def draw_subgraph(
    H: nx.Graph,
    node_types: dict[str, str],
    outpath: str,
    title: str,
    seed: int = 7,
    label_center: str | None = None,
) -> None:
    _safe_mkdir(os.path.dirname(outpath))

    # Node categories
    proteins = [n for n in H.nodes if node_types.get(n) == "protein"]
    sites = [n for n in H.nodes if node_types.get(n) == "site"]
    unknown = [n for n in H.nodes if node_types.get(n) not in ("protein", "site")]

    # Edge categories
    rel_to_edges = defaultdict(list)
    for u, v, d in H.edges(data=True):
        rel = d.get("relation", "unknown")
        rel_to_edges[rel].append((u, v))

    # Layout
    pos = nx.spring_layout(H, seed=seed, k=None)

    # Poster-friendly sizing
    plt.figure(figsize=(12, 9), dpi=300)
    plt.title(title, fontsize=18, pad=12)

    # Draw edges by relation (separate calls so legend works)
    # Use a stable order so the plot is consistent run-to-run
    relations_in_H = [r for r in RELATION_ORDER if r in rel_to_edges] + [
        r for r in rel_to_edges.keys() if r not in RELATION_ORDER
    ]

    for rel in relations_in_H:
        edgelist = rel_to_edges[rel]
        # Slightly different alpha depending on edge type
        alpha = 0.25 if rel in ("protein_same_pathway", "site_same_pathway") else 0.55
        width = 0.6 if rel in ("protein_same_pathway", "site_same_pathway") else 1.1
        nx.draw_networkx_edges(H, pos, edgelist=edgelist, alpha=alpha, width=width, label=rel)

    # Draw nodes
    nx.draw_networkx_nodes(H, pos, nodelist=proteins, node_size=55, alpha=0.90, label="protein")
    nx.draw_networkx_nodes(H, pos, nodelist=sites, node_size=22, alpha=0.85, label="site")
    if unknown:
        nx.draw_networkx_nodes(H, pos, nodelist=unknown, node_size=30, alpha=0.85, label="other")

    # Label only the center to avoid clutter
    if label_center and label_center in H:
        nx.draw_networkx_labels(H, pos, labels={label_center: label_center}, font_size=10)

    # Legend
    plt.legend(loc="upper right", fontsize=8, frameon=True)

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()


def main() -> None:
    nodes, edges = load_nodes_edges()
    node_types = dict(zip(nodes["node_id"], nodes["node_type"]))

    # Build full undirected graph for subgraph extraction
    G = build_graph(edges, undirected=True)

    # Option 1: Hub ego networks (recommended)
    for center in DEFAULT_HUBS[:3]:
        H = ego_subgraph(G, center=center, radius=1, max_neighbors=250)
        out = os.path.join(FIGURES_DIR, f"subgraph_ego_{center.replace(':', '_')}.png")
        title = f"Ego subgraph (radius=1) around {center}\nNodes={H.number_of_nodes():,} Edges={H.number_of_edges():,}"
        draw_subgraph(H, node_types, out, title, label_center=center)

    # Option 2: Phosphorylation-only subgraph (whole graph, then sample one hub)
    phospho_rels = {"phosphorylates", "dephosphorylates", "has_site"}
    edges_ph = filter_edges_by_relation(edges, phospho_rels)
    Gph = build_graph(edges_ph, undirected=True)

    # Pick a center that exists in phospho layer
    center = "PROTEIN:CDK1" if "PROTEIN:CDK1" in Gph else next(iter(Gph.nodes))
    Hph = ego_subgraph(Gph, center=center, radius=1, max_neighbors=350)
    out = os.path.join(FIGURES_DIR, f"subgraph_phospho_ego_{center.replace(':', '_')}.png")
    title = f"Phosphorylation layer ego (radius=1) around {center}\nNodes={Hph.number_of_nodes():,} Edges={Hph.number_of_edges():,}"
    draw_subgraph(Hph, node_types, out, title, label_center=center)

    # Option 3: Multi-layer ego (phospho + PPI + pathway) around CDK1
    multi_rels = {"phosphorylates", "dephosphorylates", "has_site", "ppi_high_confidence", "protein_same_pathway"}
    edges_ml = filter_edges_by_relation(edges, multi_rels)
    Gml = build_graph(edges_ml, undirected=True)

    center = "PROTEIN:CDK1" if "PROTEIN:CDK1" in Gml else next(iter(Gml.nodes))
    Hml = ego_subgraph(Gml, center=center, radius=1, max_neighbors=300)
    out = os.path.join(FIGURES_DIR, f"subgraph_multilayer_ego_{center.replace(':', '_')}.png")
    title = f"Multi-layer ego (phospho + PPI + pathway) around {center}\nNodes={Hml.number_of_nodes():,} Edges={Hml.number_of_edges():,}"
    draw_subgraph(Hml, node_types, out, title, label_center=center)

    print("Wrote subgraph figures to figures/:")
    print("  subgraph_ego_PROTEIN_CDK1.png (and 2 more hubs)")
    print("  subgraph_phospho_ego_*.png")
    print("  subgraph_multilayer_ego_*.png")


if __name__ == "__main__":
    main()
