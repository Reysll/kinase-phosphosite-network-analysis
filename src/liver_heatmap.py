# src/liver_heatmap.py
from __future__ import annotations

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.config import FIGURES_DIR, LIVER_XLSX, SHEET_PROTEIN_EXPR, SHEET_PHOSPHO


def _safe_mkdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _pick_group_columns(df: pd.DataFrame, group: str) -> list[str]:
    """
    group: "Control" or "Sample"
    Works for both tabs based on column naming.
    """
    cols = []
    for c in df.columns:
        if isinstance(c, str) and ("Abundance" in c or "Abundances" in c):
            if f": {group}" in c:
                cols.append(c)
    return cols


def _to_numeric_matrix(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    X = df[cols].copy()
    for c in cols:
        X[c] = pd.to_numeric(X[c], errors="coerce")
    return X


def _filter_rows_min_samples(X: pd.DataFrame, min_samples: int) -> pd.DataFrame:
    keep = X.notna().sum(axis=1) >= min_samples
    return X.loc[keep].copy()


def _select_top_variable_rows(
    X: pd.DataFrame,
    feature_names: pd.Series,
    top_n: int,
) -> pd.DataFrame:
    """
    Select rows with highest variance (ignoring NaNs).
    Returns matrix indexed by feature name.
    """
    var = X.var(axis=1, skipna=True)
    idx = var.sort_values(ascending=False).head(top_n).index
    Xs = X.loc[idx].copy()

    names = feature_names.loc[idx].astype(str)
    # Ensure uniqueness for heatmap labels
    names = names.where(~names.duplicated(), names + "_" + names.groupby(names).cumcount().astype(str))
    Xs.index = names.values
    return Xs


def _corr_heatmap(corr: pd.DataFrame, title: str, outpath: str) -> None:
    _safe_mkdir(os.path.dirname(outpath))

    arr = corr.values
    n = arr.shape[0]

    plt.figure(figsize=(10, 9), dpi=300)
    im = plt.imshow(arr, vmin=-1, vmax=1, aspect="auto", interpolation="nearest")
    plt.title(title, fontsize=16, pad=12)

    # For poster readability, show ticks only if not too many
    if n <= 60:
        plt.xticks(range(n), corr.columns, rotation=90, fontsize=6)
        plt.yticks(range(n), corr.index, fontsize=6)
    else:
        plt.xticks([])
        plt.yticks([])

    cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
    cbar.set_label("Pearson r", rotation=90)

    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()


def _protein_expression_feature_name(df: pd.DataFrame) -> pd.Series:
    # Preferred: Gene Symbol, fallback: Accession
    if "Gene Symbol" in df.columns:
        return df["Gene Symbol"].fillna(df.get("Accession", "NA"))
    if "GeneSymbol" in df.columns:
        return df["GeneSymbol"].fillna(df.get("Accession", "NA"))
    return df.get("Gene Symbol", df.iloc[:, 0]).astype(str)


def _phospho_feature_name(df: pd.DataFrame) -> pd.Series:
    """
    Build a readable ID like: GENE_S311 or ALB_T539 using:
      - Gene Symbol
      - Modifications in Master Proteins (contains [S311(100)] etc.)
    """
    gene = df.get("Gene Symbol", pd.Series(["NA"] * len(df))).astype(str)

    mod_col = None
    for c in df.columns:
        if isinstance(c, str) and "Modifications" in c:
            mod_col = c
            break

    if mod_col is None:
        # fallback to first column
        return gene

    mods = df[mod_col].astype(str)

    # Extract first site token like S311, T539, Y476
    site = mods.str.extract(r"\[.*?([STY][0-9]+)\(", expand=False)
    site = site.fillna("UNK")

    return (gene + "_" + site).astype(str)


def run_heatmaps(
    sheet: str,
    feature_mode: str,
    top_n: int = 50,
    min_samples: int = 6,
) -> None:
    """
    sheet: SHEET_PROTEIN_EXPR or SHEET_PHOSPHO
    feature_mode: "protein" or "phospho"
    """
    df = pd.read_excel(LIVER_XLSX, sheet_name=sheet, header=1)

    if feature_mode == "protein":
        feat = _protein_expression_feature_name(df)
        prefix = "liver_protein"
    else:
        feat = _phospho_feature_name(df)
        prefix = "liver_phospho"

    for group in ["Control", "Sample"]:
        cols = _pick_group_columns(df, group)
        if not cols:
            raise ValueError(f"No abundance columns found for group={group} in sheet={sheet}")

        X = _to_numeric_matrix(df, cols)
        X = _filter_rows_min_samples(X, min_samples=min_samples)

        # Select top N most variable rows within this group
        Xtop = _select_top_variable_rows(X, feat, top_n=top_n)

        # Correlation across samples: corr between features (rows)
        corr = Xtop.T.corr(method="pearson", min_periods=min_samples)

        out = os.path.join(FIGURES_DIR, f"{prefix}_corr_heatmap_{group.lower()}_top{top_n}.png")
        title = f"{prefix.replace('_', ' ').title()} correlation heatmap ({group})\nTop {top_n} variable, min_samples={min_samples}"
        _corr_heatmap(corr, title, out)

    print("Saved heatmaps to figures/:")
    print(f"  {prefix}_corr_heatmap_control_top{top_n}.png")
    print(f"  {prefix}_corr_heatmap_sample_top{top_n}.png")


def main() -> None:
    # Protein expression heatmaps
    run_heatmaps(sheet=SHEET_PROTEIN_EXPR, feature_mode="protein", top_n=50, min_samples=6)

    # Phosphorylation heatmaps
    run_heatmaps(sheet=SHEET_PHOSPHO, feature_mode="phospho", top_n=50, min_samples=6)


if __name__ == "__main__":
    main()
