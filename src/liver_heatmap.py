# src/liver_heatmap.py
from __future__ import annotations

import re
from pathlib import Path
from typing import Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.config import (
    LIVER_XLSX,
    SHEET_PROTEIN_EXPR,
    SHEET_PHOSPHO,
    FIGURES_LIVER_DIR,
)


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _split_cols(df: pd.DataFrame) -> Tuple[List[str], List[str]]:
    # Works for both tabs based on column text containing "Control" or "Sample"
    control_cols = [c for c in df.columns if isinstance(c, str) and "Control" in c]
    sample_cols = [c for c in df.columns if isinstance(c, str) and "Sample" in c]
    return control_cols, sample_cols


def _filter_rows_min_nonnull(df: pd.DataFrame, cols: List[str], min_n: int = 6) -> pd.DataFrame:
    nonnull = df[cols].notna().sum(axis=1)
    return df.loc[nonnull >= min_n].copy()


def _to_matrix(df: pd.DataFrame, id_col: str, value_cols: List[str]) -> pd.DataFrame:
    # Build a matrix: rows = features, columns = samples
    mat = df[[id_col] + value_cols].copy()
    mat = mat.set_index(id_col)
    # Ensure numeric
    mat = mat.apply(pd.to_numeric, errors="coerce")
    return mat


def _top_variable_rows(mat: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    # Use variance across samples (ignoring NaN) to pick features
    var = mat.var(axis=1, skipna=True)
    var = var.replace([np.inf, -np.inf], np.nan).dropna()
    if len(var) <= top_n:
        return mat.loc[var.index]
    keep = var.sort_values(ascending=False).head(top_n).index
    return mat.loc[keep]


def _pairwise_corr_with_min_overlap(mat: pd.DataFrame, min_overlap: int = 6) -> pd.DataFrame:
    """
    Compute correlation matrix with pairwise deletion and a minimum overlap count.
    If overlap < min_overlap, set corr to NaN for that pair.
    """
    X = mat.to_numpy(dtype=float)
    n = X.shape[0]
    corr = np.full((n, n), np.nan, dtype=float)

    for i in range(n):
        corr[i, i] = 1.0
        xi = X[i, :]
        for j in range(i + 1, n):
            xj = X[j, :]
            mask = np.isfinite(xi) & np.isfinite(xj)
            m = int(mask.sum())
            if m < min_overlap:
                continue
            a = xi[mask]
            b = xj[mask]
            # If constant vectors, corr is undefined
            if np.std(a) == 0 or np.std(b) == 0:
                continue
            r = np.corrcoef(a, b)[0, 1]
            corr[i, j] = r
            corr[j, i] = r

    return pd.DataFrame(corr, index=mat.index, columns=mat.index)


def _plot_heatmap(corr_df: pd.DataFrame, title: str, outpath: Path) -> None:
    fig_w = 12
    fig_h = 10
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=200)
    ax = plt.gca()

    data = corr_df.to_numpy(dtype=float)

    im = ax.imshow(data, aspect="auto", interpolation="nearest")
    ax.set_title(title)

    ax.set_xticks(range(len(corr_df.columns)))
    ax.set_yticks(range(len(corr_df.index)))
    ax.set_xticklabels(corr_df.columns, rotation=90, fontsize=6)
    ax.set_yticklabels(corr_df.index, fontsize=6)

    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Pearson r")

    plt.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


def _load_protein_expression() -> pd.DataFrame:
    # header=1 because row0 is a description line
    df = pd.read_excel(LIVER_XLSX, sheet_name=SHEET_PROTEIN_EXPR, header=1)
    return df


def _load_phospho() -> pd.DataFrame:
    df = pd.read_excel(LIVER_XLSX, sheet_name=SHEET_PHOSPHO, header=1)
    return df


def main() -> None:
    _ensure_dir(FIGURES_LIVER_DIR)

    # Protein expression
    pe = _load_protein_expression()
    pe_control_cols, pe_sample_cols = _split_cols(pe)

    # Use Gene Symbol if present
    pe_id_col = "Gene Symbol" if "Gene Symbol" in pe.columns else pe.columns[0]

    pe_control = _filter_rows_min_nonnull(pe, pe_control_cols, min_n=6)
    pe_sample = _filter_rows_min_nonnull(pe, pe_sample_cols, min_n=6)

    pe_mat_control = _to_matrix(pe_control, pe_id_col, pe_control_cols)
    pe_mat_sample = _to_matrix(pe_sample, pe_id_col, pe_sample_cols)

    pe_mat_control = _top_variable_rows(pe_mat_control, top_n=50)
    pe_mat_sample = _top_variable_rows(pe_mat_sample, top_n=50)

    pe_corr_control = _pairwise_corr_with_min_overlap(pe_mat_control, min_overlap=6)
    pe_corr_sample = _pairwise_corr_with_min_overlap(pe_mat_sample, min_overlap=6)

    _plot_heatmap(
        pe_corr_control,
        f"Protein expression correlation heatmap (Control) top50 | rows={pe_mat_control.shape[0]}",
        FIGURES_LIVER_DIR / "liver_protein_corr_heatmap_control_top50.png",
    )
    _plot_heatmap(
        pe_corr_sample,
        f"Protein expression correlation heatmap (Sample) top50 | rows={pe_mat_sample.shape[0]}",
        FIGURES_LIVER_DIR / "liver_protein_corr_heatmap_sample_top50.png",
    )

    # Phosphorylation
    ph = _load_phospho()
    ph_control_cols, ph_sample_cols = _split_cols(ph)

    # Use "Modifications in Master Proteins" if present (acts as a phosphosite id)
    if "Modifications in Master Proteins" in ph.columns:
        ph_id_col = "Modifications in Master Proteins"
    elif "Gene Symbol" in ph.columns:
        ph_id_col = "Gene Symbol"
    else:
        ph_id_col = ph.columns[0]

    ph_control = _filter_rows_min_nonnull(ph, ph_control_cols, min_n=6)
    ph_sample = _filter_rows_min_nonnull(ph, ph_sample_cols, min_n=6)

    ph_mat_control = _to_matrix(ph_control, ph_id_col, ph_control_cols)
    ph_mat_sample = _to_matrix(ph_sample, ph_id_col, ph_sample_cols)

    ph_mat_control = _top_variable_rows(ph_mat_control, top_n=50)
    ph_mat_sample = _top_variable_rows(ph_mat_sample, top_n=50)

    ph_corr_control = _pairwise_corr_with_min_overlap(ph_mat_control, min_overlap=6)
    ph_corr_sample = _pairwise_corr_with_min_overlap(ph_mat_sample, min_overlap=6)

    _plot_heatmap(
        ph_corr_control,
        f"Phosphosite correlation heatmap (Control) top50 | rows={ph_mat_control.shape[0]}",
        FIGURES_LIVER_DIR / "liver_phospho_corr_heatmap_control_top50.png",
    )
    _plot_heatmap(
        ph_corr_sample,
        f"Phosphosite correlation heatmap (Sample) top50 | rows={ph_mat_sample.shape[0]}",
        FIGURES_LIVER_DIR / "liver_phospho_corr_heatmap_sample_top50.png",
    )

    print("Saved heatmaps to:", str(FIGURES_LIVER_DIR))


if __name__ == "__main__":
    main()
