from __future__ import annotations

import os
import re
from typing import Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.config import LIVER_XLSX_PATH, FIGURES_DIR


def _ensure_dirs() -> None:
    os.makedirs(FIGURES_DIR, exist_ok=True)


def _load_sheet_with_row0_header(xlsx_path: str, sheet_name: str) -> pd.DataFrame:
    """
    Your liver Excel sheets have the true header in row 1 of the sheet content.
    We read normally, then promote the first row to header.
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet_name)
    header = df.iloc[0].astype(str).tolist()
    df = df.iloc[1:].copy()
    df.columns = header
    df = df.reset_index(drop=True)
    return df


def _split_control_case_columns(cols: List[str]) -> Tuple[List[str], List[str]]:
    """
    The sheet columns include strings like:
      "Abundances (Normalized): F26: Control"
      "Abundance: F1: Control"
      "Abundance: F18: Sample"  (Sample == case in tutor email)
    """
    control_cols = []
    case_cols = []

    for c in cols:
        c_str = str(c)
        if re.search(r"\bControl\b", c_str, flags=re.IGNORECASE):
            control_cols.append(c)
        elif re.search(r"\bSample\b", c_str, flags=re.IGNORECASE):
            case_cols.append(c)

    return control_cols, case_cols


def _pairwise_correlations(
    mat: pd.DataFrame,
    min_common: int = 6,
) -> np.ndarray:
    """
    Compute correlations between all pairs of rows, using:
      - only common non-NaN samples per pair
      - if common samples < min_common: skip
    Returns an array of correlation values (r).
    """
    values = mat.to_numpy(dtype=float)
    n = values.shape[0]
    corrs = []

    for i in range(n):
        xi = values[i]
        for j in range(i + 1, n):
            xj = values[j]
            mask = np.isfinite(xi) & np.isfinite(xj)
            if mask.sum() < min_common:
                continue
            r = np.corrcoef(xi[mask], xj[mask])[0, 1]
            if np.isfinite(r):
                corrs.append(r)

    return np.array(corrs, dtype=float)


def analyze_liver_protein_expression() -> None:
    _ensure_dirs()

    df = _load_sheet_with_row0_header(LIVER_XLSX_PATH, "ProteinExpression")

    # Identify columns
    all_cols = df.columns.tolist()
    control_cols, case_cols = _split_control_case_columns(all_cols)

    # Choose a stable ID column for row identity
    # In your file, these exist: "Gene Symbol" is a good label.
    gene_col = "Gene Symbol" if "Gene Symbol" in df.columns else None
    if gene_col is None:
        raise ValueError("Could not find 'Gene Symbol' column in ProteinExpression sheet.")

    # Convert measurement columns to numeric
    df_control = df[[gene_col] + control_cols].copy()
    df_case = df[[gene_col] + case_cols].copy()

    for c in control_cols:
        df_control[c] = pd.to_numeric(df_control[c], errors="coerce")
    for c in case_cols:
        df_case[c] = pd.to_numeric(df_case[c], errors="coerce")

    # Filter: protein must be identified in >=6 samples within the category
    control_ok = df_control[control_cols].notna().sum(axis=1) >= 6
    case_ok = df_case[case_cols].notna().sum(axis=1) >= 6

    df_control = df_control.loc[control_ok].reset_index(drop=True)
    df_case = df_case.loc[case_ok].reset_index(drop=True)

    # Matrix (rows=proteins, cols=samples)
    mat_control = df_control.set_index(gene_col)[control_cols]
    mat_case = df_case.set_index(gene_col)[case_cols]

    print("ProteinExpression:")
    print("  proteins_control_kept:", mat_control.shape[0])
    print("  proteins_case_kept:", mat_case.shape[0])
    print("  control_samples:", mat_control.shape[1])
    print("  case_samples:", mat_case.shape[1])

    # Correlations with pairwise common-sample >=6 rule
    corr_control = _pairwise_correlations(mat_control, min_common=6)
    corr_case = _pairwise_correlations(mat_case, min_common=6)

    pd.DataFrame({"r": corr_control}).to_csv(os.path.join(FIGURES_DIR, "liver_protein_corr_control.csv"), index=False)
    pd.DataFrame({"r": corr_case}).to_csv(os.path.join(FIGURES_DIR, "liver_protein_corr_case.csv"), index=False)

    # Histograms
    plt.figure()
    plt.hist(corr_control, bins=60)
    plt.xlabel("Pearson r")
    plt.ylabel("Count")
    plt.title("Protein expression correlations (Control)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "liver_protein_corr_control_hist.png"), dpi=200)
    plt.close()

    plt.figure()
    plt.hist(corr_case, bins=60)
    plt.xlabel("Pearson r")
    plt.ylabel("Count")
    plt.title("Protein expression correlations (Case)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "liver_protein_corr_case_hist.png"), dpi=200)
    plt.close()

    print("Wrote liver protein correlation outputs to /figures.")


def main() -> None:
    analyze_liver_protein_expression()
    print("Done.")


if __name__ == "__main__":
    main()
