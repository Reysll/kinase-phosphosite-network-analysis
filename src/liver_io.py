# src/liver_io.py
from __future__ import annotations

import pandas as pd

from src.config import LIVER_XLSX, SHEET_PROTEIN_EXPR, SHEET_PHOSPHO


def load_liver_protein_expression() -> pd.DataFrame:
    # first row is a description, real headers start on row 2
    return pd.read_excel(LIVER_XLSX, sheet_name=SHEET_PROTEIN_EXPR, header=1)


def load_liver_phosphorylation() -> pd.DataFrame:
    return pd.read_excel(LIVER_XLSX, sheet_name=SHEET_PHOSPHO, header=1)
