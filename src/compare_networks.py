# src/compare_networks.py
from __future__ import annotations

from pathlib import Path
import sys
import pandas as pd

from src.config import DATA_DIR, LIVER_NODES_PATH, LIVER_EDGES_PATH, FIGURES_LIVER_DIR


def _read_csv_any(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path)


def _summarize(nodes: pd.DataFrame, edges: pd.DataFrame, label: str) -> dict:
    if "node_type" not in nodes.columns:
        raise ValueError(f"{label} nodes missing 'node_type' column. Columns: {list(nodes.columns)}")
    if not {"source", "target", "relation"}.issubset(edges.columns):
        raise ValueError(f"{label} edges missing required cols. Columns: {list(edges.columns)}")

    return {
        "network": label,
        "num_nodes": int(len(nodes)),
        "num_edges": int(len(edges)),
        "proteins": int((nodes["node_type"] == "protein").sum()),
        "sites": int((nodes["node_type"] == "site").sum()),
        "edge_types": int(edges["relation"].nunique()),
    }


def main() -> None:
    print("[compare_networks] starting...", flush=True)

    generic_nodes_path = DATA_DIR / "nodes.csv.gz"
    generic_edges_path = DATA_DIR / "edges.csv.gz"

    print("[compare_networks] paths:", flush=True)
    print("  generic_nodes:", generic_nodes_path, flush=True)
    print("  generic_edges:", generic_edges_path, flush=True)
    print("  liver_nodes:  ", LIVER_NODES_PATH, flush=True)
    print("  liver_edges:  ", LIVER_EDGES_PATH, flush=True)

    g_nodes = _read_csv_any(generic_nodes_path)
    g_edges = _read_csv_any(generic_edges_path)
    l_nodes = _read_csv_any(LIVER_NODES_PATH)
    l_edges = _read_csv_any(LIVER_EDGES_PATH)

    print(f"[compare_networks] loaded rows: generic nodes={len(g_nodes):,} edges={len(g_edges):,}", flush=True)
    print(f"[compare_networks] loaded rows: liver   nodes={len(l_nodes):,} edges={len(l_edges):,}", flush=True)

    s1 = _summarize(g_nodes, g_edges, "generic")
    s2 = _summarize(l_nodes, l_edges, "liver")

    out = pd.DataFrame([s1, s2])
    delta = pd.DataFrame([{
        "network": "delta(liver-generic)",
        "num_nodes": s2["num_nodes"] - s1["num_nodes"],
        "num_edges": s2["num_edges"] - s1["num_edges"],
        "proteins": s2["proteins"] - s1["proteins"],
        "sites": s2["sites"] - s1["sites"],
        "edge_types": s2["edge_types"] - s1["edge_types"],
    }])

    out = pd.concat([out, delta], ignore_index=True)

    print("\n[compare_networks] summary:", flush=True)
    print(out.to_string(index=False), flush=True)

    FIGURES_LIVER_DIR.mkdir(parents=True, exist_ok=True)
    out_path = FIGURES_LIVER_DIR / "generic_vs_liver_summary.csv"
    out.to_csv(out_path, index=False)
    print(f"\n[compare_networks] wrote: {out_path}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("[compare_networks] ERROR:", repr(e), file=sys.stderr, flush=True)
        raise
