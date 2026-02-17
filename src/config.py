# src/config.py
from __future__ import annotations

import os

# Project root = one folder above /src
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

DATA_DIR = os.path.join(PROJECT_ROOT, "data")
FIGURES_DIR = os.path.join(PROJECT_ROOT, "figures")

# Generic network outputs (from kinase-phosphosite-graph repo)
# Put your nodes/edges .csv.gz inside data/ or point these to the right path.
NODES_PATH = os.path.join(DATA_DIR, "nodes.csv.gz")
EDGES_PATH = os.path.join(DATA_DIR, "edges.csv.gz")

# Liver excel
LIVER_XLSX = os.path.join(DATA_DIR, "LiverCancer_ProtExp_Phospho_casecntrl.xlsx")

# Standardized sheet names (as you requested)
SHEET_PROTEIN_EXPR = "ProteinExpression"
SHEET_PHOSPHO = "Phosphorylation"
