# src/config.py
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import List
# src/config.py
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_ROOT / "data"
FIGURES_DIR = PROJECT_ROOT / "figures"

# Liver-specific
LIVER_DIR = DATA_DIR / "liver"
LIVER_NODES_PATH = LIVER_DIR / "Liver_network_nodes.csv.gz"
LIVER_EDGES_PATH = LIVER_DIR / "Liver_network_edges.csv.gz"

FIGURES_LIVER_DIR = FIGURES_DIR / "figures_liver"

# Liver Excel
LIVER_XLSX = DATA_DIR / "LiverCancer_ProtExp_Phospho_casecntrl.xlsx"
SHEET_PROTEIN_EXPR = "ProteinExpression"
SHEET_PHOSPHO = "Phosphorylation"


def _first_existing(paths: List[str]) -> str:
    for p in paths:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        "None of these paths exist:\n" + "\n".join(paths)
    )


@dataclass(frozen=True)
class NetworkPaths:
    nodes: str
    edges: str


def get_generic_paths() -> NetworkPaths:
    # Edit these if you store generic somewhere else
    candidates_nodes = [
        os.path.join(DATA_DIR, "generic", "nodes.csv.gz"),
        os.path.join(DATA_DIR, "nodes.csv.gz"),
        os.path.join(PROJECT_ROOT, "outputs", "nodes.csv.gz"),
    ]
    candidates_edges = [
        os.path.join(DATA_DIR, "generic", "edges.csv.gz"),
        os.path.join(DATA_DIR, "edges.csv.gz"),
        os.path.join(PROJECT_ROOT, "outputs", "edges.csv.gz"),
    ]
    return NetworkPaths(
        nodes=_first_existing(candidates_nodes),
        edges=_first_existing(candidates_edges),
    )


def get_liver_paths() -> NetworkPaths:
    # Edit these if you store liver somewhere else
    candidates_nodes = [
        os.path.join(DATA_DIR, "liver", "nodes.csv.gz"),
        os.path.join(DATA_DIR, "nodes_liver.csv.gz"),
        os.path.join(PROJECT_ROOT, "outputs_liver", "nodes.csv.gz"),
        os.path.join(PROJECT_ROOT, "figures_liver", "nodes.csv.gz"),
    ]
    candidates_edges = [
        os.path.join(DATA_DIR, "liver", "edges.csv.gz"),
        os.path.join(DATA_DIR, "edges_liver.csv.gz"),
        os.path.join(PROJECT_ROOT, "outputs_liver", "edges.csv.gz"),
        os.path.join(PROJECT_ROOT, "figures_liver", "edges.csv.gz"),
    ]
    return NetworkPaths(
        nodes=_first_existing(candidates_nodes),
        edges=_first_existing(candidates_edges),
    )
