# src/config.py
from __future__ import annotations
import os

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

DATA_DIR = os.path.join(PROJECT_ROOT, "data")
FIGURES_DIR = os.path.join(PROJECT_ROOT, "figures")

NODES_PATH = os.path.join(DATA_DIR, "nodes.csv.gz")
EDGES_PATH = os.path.join(DATA_DIR, "edges.csv.gz")

LIVER_XLSX_PATH = os.path.join(DATA_DIR, "LiverCancer_ProtExp_Phospho_casecntrl.xlsx")
