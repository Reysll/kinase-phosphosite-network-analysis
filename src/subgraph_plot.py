# src/subgraph_plot.py
from __future__ import annotations

import argparse
import os
import random
from typing import Optional, Set, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


print("[subgraph_plot] module imported")


def _default_paths(liver: bool) -> Tuple[str, str, str]:
    """
    Returns (nodes_path, edges_path, out_dir)
    Uses local repo conventions shown in your screenshot:
      - generic: data/nodes.csv.gz, data/edges.csv.gz, figures/
      - liver:   data/liver/Liver_network_nodes.csv.gz, data/liver/Liver_network_edges.csv.gz, figures/figures_liver/
    """
    if liver:
        nodes_path = os.path.join("data", "liver", "Liver_network_nodes.csv.gz")
        edges_path = os.path.join("data", "liver", "Liver_network_edges.csv.gz")
        out_dir = os.path.join("figures", "figures_liver")
    else:
        nodes_path = os.path.join("data", "nodes.csv.gz")
        edges_path = os.path.join("data", "edges.csv.gz")
        out_dir = os.path.join("figures")

    return nodes_path, edges_path, out_dir


def load_network(nodes_path: str, edges_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    nodes = pd.read_csv(nodes_path)
    edges = pd.read_csv(edges_path)
    # Normalize expected column names
    # nodes: node_id,node_type
    # edges: source,target,relation,(optional)confidence_score
    required_nodes = {"node_id", "node_type"}
    required_edges = {"source", "target", "relation"}
    if not required_nodes.issubset(set(nodes.columns)):
        raise ValueError(f"nodes file missing columns: {required_nodes - set(nodes.columns)}")
    if not required_edges.issubset(set(edges.columns)):
        raise ValueError(f"edges file missing columns: {required_edges - set(edges.columns)}")
    return nodes, edges


def build_graph_undirected(edges: pd.DataFrame, allowed_relations: Optional[Set[str]] = None) -> nx.Graph:
    if allowed_relations is not None:
        edges = edges[edges["relation"].isin(sorted(allowed_relations))].copy()

    G = nx.Graph()
    # Add edges
    for row in edges.itertuples(index=False):
        # row.source, row.target, row.relation, maybe row.confidence_score
        src = getattr(row, "source")
        tgt = getattr(row, "target")
        rel = getattr(row, "relation")
        # Keep relation as attribute so we can optionally style later
        G.add_edge(src, tgt, relation=rel)
    return G


def cap_ego_graph(
    H: nx.Graph,
    seed_node: str,
    max_nodes: int,
    rng_seed: int = 7,
) -> nx.Graph:
    """
    Downsample nodes in an ego subgraph to keep the plot poster-friendly.
    Strategy:
      - Always keep seed
      - Preferentially keep highest-degree neighbors first
      - If still too many, randomly sample remaining
    """
    if max_nodes <= 0 or H.number_of_nodes() <= max_nodes:
        return H

    rng = random.Random(rng_seed)

    if seed_node not in H:
        return H

    # Sort nodes by degree (desc), keep seed at top
    degrees = dict(H.degree())
    others = [n for n in H.nodes() if n != seed_node]
    others_sorted = sorted(others, key=lambda n: degrees.get(n, 0), reverse=True)

    keep = [seed_node]

    # Keep top-degree nodes first
    for n in others_sorted:
        if len(keep) >= max_nodes:
            break
        keep.append(n)

    # If we somehow still need random (should not), do it
    if len(keep) < max_nodes and len(others) > len(keep) - 1:
        remaining = [n for n in others if n not in keep]
        need = max_nodes - len(keep)
        if need > 0 and len(remaining) > 0:
            keep.extend(rng.sample(remaining, k=min(need, len(remaining))))

    return H.subgraph(keep).copy()


def cap_edges(H: nx.Graph, max_edges: int, rng_seed: int = 7) -> nx.Graph:
    """
    If the subgraph is still too dense, downsample edges while keeping connectivity as much as possible.
    Approach:
      - Keep a spanning tree first (preserves connectivity if graph is connected)
      - Then add a random sample of remaining edges up to max_edges
    """
    if max_edges <= 0 or H.number_of_edges() <= max_edges:
        return H

    rng = random.Random(rng_seed)

    if H.number_of_nodes() <= 1:
        return H

    # Start with spanning tree on each component
    edges_keep = set()
    for comp_nodes in nx.connected_components(H):
        comp = H.subgraph(comp_nodes)
        if comp.number_of_nodes() <= 1:
            continue
        T = nx.minimum_spanning_tree(comp)  # unweighted spanning tree
        edges_keep.update(tuple(sorted(e)) for e in T.edges())

    # Add additional edges randomly
    all_edges = [tuple(sorted(e)) for e in H.edges()]
    remaining = [e for e in all_edges if e not in edges_keep]

    target = max_edges
    edges_keep_list = list(edges_keep)
    if len(edges_keep_list) >= target:
        edges_keep_list = edges_keep_list[:target]
    else:
        need = target - len(edges_keep_list)
        if need > 0 and len(remaining) > 0:
            edges_keep_list.extend(rng.sample(remaining, k=min(need, len(remaining))))

    H2 = nx.Graph()
    H2.add_nodes_from(H.nodes())
    for u, v in edges_keep_list:
        # preserve relation attribute if exists
        rel = H.get_edge_data(u, v).get("relation") if H.has_edge(u, v) else None
        if rel is None:
            H2.add_edge(u, v)
        else:
            H2.add_edge(u, v, relation=rel)

    # Remove isolated nodes except seed will remain if present
    return H2


def node_label_map(H: nx.Graph, top_labels: int, seed_node: str) -> dict:
    """
    Label only the top-degree nodes + seed node for readability.
    """
    if top_labels <= 0:
        return {seed_node: seed_node} if seed_node in H else {}

    deg = dict(H.degree())
    top = sorted(H.nodes(), key=lambda n: deg.get(n, 0), reverse=True)[:top_labels]
    labels = {n: n for n in top}
    if seed_node in H:
        labels[seed_node] = seed_node
    return labels


def draw_subgraph(
    H: nx.Graph,
    title: str,
    outpath: str,
    seed_node: str,
    top_labels: int = 20,
) -> None:
    # Layout that does NOT require SciPy
    # Works well for 50-300 nodes
    pos = nx.kamada_kawai_layout(H)

    plt.figure(figsize=(12, 9), dpi=200)

    # Draw nodes and edges
    nx.draw_networkx_edges(H, pos, alpha=0.25, width=0.6)

    # Highlight seed node
    node_sizes = []
    for n in H.nodes():
        if n == seed_node:
            node_sizes.append(220)
        else:
            node_sizes.append(40)

    nx.draw_networkx_nodes(H, pos, node_size=node_sizes)

    labels = node_label_map(H, top_labels=top_labels, seed_node=seed_node)
    if labels:
        nx.draw_networkx_labels(H, pos, labels=labels, font_size=7)

    plt.title(title, fontsize=16)
    plt.axis("off")
    plt.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()


def main() -> None:
    print("[subgraph_plot] starting...")
    ap = argparse.ArgumentParser(description="Plot poster-friendly ego subgraphs from generic or liver network.")
    ap.add_argument("--seed", required=True, help="Seed node id, e.g. PROTEIN:EGFR or PROTEIN:CDK1")
    ap.add_argument("--k", type=int, default=1, help="Ego radius (default 1)")
    ap.add_argument("--liver", action="store_true", help="Use liver network files under data/liver/")
    ap.add_argument("--max-nodes", type=int, default=200, help="Cap subgraph nodes for readability (default 200)")
    ap.add_argument("--max-edges", type=int, default=8000, help="Cap subgraph edges for readability (default 8000)")
    ap.add_argument("--top-labels", type=int, default=25, help="Label only top-degree nodes (default 25)")
    ap.add_argument(
        "--relations",
        default="",
        help=(
            "Optional comma-separated list of relations to include. "
            "Example: phosphorylates,dephosphorylates,has_site,ppi_high_confidence"
        ),
    )
    args = ap.parse_args()

    nodes_path, edges_path, out_dir = _default_paths(liver=args.liver)

    print(f"[subgraph_plot] network={'liver' if args.liver else 'generic'}")
    print(f"[subgraph_plot] nodes={os.path.abspath(nodes_path)}")
    print(f"[subgraph_plot] edges={os.path.abspath(edges_path)}")

    nodes_df, edges_df = load_network(nodes_path, edges_path)

    allowed_relations = None
    if args.relations.strip():
        allowed_relations = set(r.strip() for r in args.relations.split(",") if r.strip())

    G = build_graph_undirected(edges_df, allowed_relations=allowed_relations)
    print(f"[subgraph_plot] loaded: |V|={G.number_of_nodes():,} |E|={G.number_of_edges():,}")

    if args.seed not in G:
        # Give a helpful hint: common issue is wrong prefix/case
        sample = [n for n in G.nodes() if n.endswith(args.seed.split(":")[-1])]
        hint = f" Similar nodes: {sample[:10]}" if sample else ""
        raise SystemExit(f"[subgraph_plot] seed node not found: {args.seed}.{hint}")

    # Ego subgraph
    H = nx.ego_graph(G, args.seed, radius=args.k)
    print(f"[subgraph_plot] ego radius k={args.k}: |V|={H.number_of_nodes():,} |E|={H.number_of_edges():,}")

    # Cap nodes and edges for readability
    H = cap_ego_graph(H, seed_node=args.seed, max_nodes=args.max_nodes)
    if args.max_edges and H.number_of_edges() > args.max_edges:
        H = cap_edges(H, max_edges=args.max_edges)

    print(f"[subgraph_plot] after caps: |V|={H.number_of_nodes():,} |E|={H.number_of_edges():,}")

    safe_seed = args.seed.replace(":", "_").replace("/", "_")
    rel_tag = ""
    if allowed_relations:
        rel_tag = "_rels-" + "-".join(sorted(allowed_relations))[:80]

    outname = f"subgraph_ego_{safe_seed}_k{args.k}_n{H.number_of_nodes()}{rel_tag}.png"
    outpath = os.path.join(out_dir, outname)

    title = f"Ego subgraph (k={args.k}) for {args.seed} ({'liver' if args.liver else 'generic'})"
    draw_subgraph(H, title=title, outpath=outpath, seed_node=args.seed, top_labels=args.top_labels)

    print(f"[subgraph_plot] saved: {os.path.abspath(outpath)}")


if __name__ == "__main__":
    main()
