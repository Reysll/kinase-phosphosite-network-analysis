# src/liver_correlations.py
from __future__ import annotations

import os
import math
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.config import FIGURES_DIR
from src.liver_io import load_liver_protein_expression, load_liver_phosphorylation


def _split_cols_by_group(columns: List[str]) -> Tuple[List[str], List[str]]:
    """
    Split columns into Control vs Sample based on suffix in the column name.
    Works for both tabs.
    """
    control_cols: List[str] = []
    sample_cols: List[str] = []
    for c in columns:
        c_str = str(c)
        if ": Control" in c_str:
            control_cols.append(c)
        elif ": Sample" in c_str:
            sample_cols.append(c)
    return control_cols, sample_cols


def _filter_rows_min_obs(mat: pd.DataFrame, min_obs: int = 6) -> pd.DataFrame:
    keep = mat.notna().sum(axis=1) >= min_obs
    return mat.loc[keep].copy()


def _corr_pairwise_with_missing(x: np.ndarray, y: np.ndarray, min_common: int = 6) -> Optional[float]:
    mask = np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) < min_common:
        return None
    xv = x[mask]
    yv = y[mask]
    if np.nanstd(xv) == 0 or np.nanstd(yv) == 0:
        return None
    return float(np.corrcoef(xv, yv)[0, 1])


def _sample_pairs(n: int, max_pairs: int, seed: int = 7) -> List[Tuple[int, int]]:
    """
    Randomly sample index pairs (i, j) with i < j.
    If all pairs <= max_pairs, return all pairs.
    """
    rng = np.random.default_rng(seed)
    total_pairs = n * (n - 1) // 2

    if total_pairs <= max_pairs:
        return [(i, j) for i in range(n) for j in range(i + 1, n)]

    idxs = rng.choice(total_pairs, size=max_pairs, replace=False)

    def unrank(k: int) -> Tuple[int, int]:
        # Unrank k in lex order of pairs to (i,j)
        i = int(n - 2 - math.floor(math.sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5))
        before = i * (2 * n - i - 1) // 2
        j = int(k - before + i + 1)
        return i, j

    return [unrank(int(k)) for k in idxs]


def _compute_corr_distribution(
    mat: pd.DataFrame,
    max_pairs: int = 250_000,
    min_common: int = 6,
    seed: int = 7,
) -> np.ndarray:
    mat_np = mat.to_numpy(dtype=float)
    n = mat_np.shape[0]
    if n < 2:
        return np.array([], dtype=float)

    pairs = _sample_pairs(n, max_pairs=max_pairs, seed=seed)
    out: List[float] = []

    for i, j in pairs:
        r = _corr_pairwise_with_missing(mat_np[i, :], mat_np[j, :], min_common=min_common)
        if r is not None and np.isfinite(r):
            out.append(r)

    return np.array(out, dtype=float)


def _save_hist(values: np.ndarray, title: str, out_path: str, bins: int = 60) -> None:
    plt.figure()
    if values.size == 0:
        plt.text(0.5, 0.5, "No correlations computed", ha="center", va="center")
        plt.title(title)
        plt.axis("off")
    else:
        plt.hist(values, bins=bins)
        plt.title(title)
        plt.xlabel("Pearson correlation (r)")
        plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def run_liver_correlation_figures(
    min_obs: int = 6,
    min_common: int = 6,
    max_pairs: int = 250_000,
    seed: int = 7,
) -> None:
    """
    Main entrypoint used by src/liver_plots.py

    Produces:
      figures/liver_protein_corr_control_hist.png
      figures/liver_protein_corr_sample_hist.png
      figures/liver_phospho_corr_control_hist.png
      figures/liver_phospho_corr_sample_hist.png
      figures/liver_correlation_summary.csv
    """
    os.makedirs(FIGURES_DIR, exist_ok=True)

    # Protein expression
    df_prot = load_liver_protein_expression()
    control_cols, sample_cols = _split_cols_by_group(list(df_prot.columns))
    if not control_cols or not sample_cols:
        raise ValueError("Could not find Control/Sample columns in ProteinExpression sheet.")

    prot_control = _filter_rows_min_obs(df_prot[control_cols].copy(), min_obs=min_obs)
    prot_sample = _filter_rows_min_obs(df_prot[sample_cols].copy(), min_obs=min_obs)

    r_prot_ctrl = _compute_corr_distribution(prot_control, max_pairs=max_pairs, min_common=min_common, seed=seed)
    r_prot_samp = _compute_corr_distribution(prot_sample, max_pairs=max_pairs, min_common=min_common, seed=seed)

    _save_hist(
        r_prot_ctrl,
        f"Protein expression correlations (Control)\nrows_kept={len(prot_control):,}, pairs_sampled<= {max_pairs:,}",
        os.path.join(FIGURES_DIR, "liver_protein_corr_control_hist.png"),
    )
    _save_hist(
        r_prot_samp,
        f"Protein expression correlations (Sample)\nrows_kept={len(prot_sample):,}, pairs_sampled<= {max_pairs:,}",
        os.path.join(FIGURES_DIR, "liver_protein_corr_sample_hist.png"),
    )

    # Phosphorylation
    df_ph = load_liver_phosphorylation()
    control_cols2, sample_cols2 = _split_cols_by_group(list(df_ph.columns))
    if not control_cols2 or not sample_cols2:
        raise ValueError("Could not find Control/Sample columns in Phosphorylation sheet.")

    ph_control = _filter_rows_min_obs(df_ph[control_cols2].copy(), min_obs=min_obs)
    ph_sample = _filter_rows_min_obs(df_ph[sample_cols2].copy(), min_obs=min_obs)

    r_ph_ctrl = _compute_corr_distribution(ph_control, max_pairs=max_pairs, min_common=min_common, seed=seed)
    r_ph_samp = _compute_corr_distribution(ph_sample, max_pairs=max_pairs, min_common=min_common, seed=seed)

    _save_hist(
        r_ph_ctrl,
        f"Phosphosite correlations (Control)\nrows_kept={len(ph_control):,}, pairs_sampled<= {max_pairs:,}",
        os.path.join(FIGURES_DIR, "liver_phospho_corr_control_hist.png"),
    )
    _save_hist(
        r_ph_samp,
        f"Phosphosite correlations (Sample)\nrows_kept={len(ph_sample):,}, pairs_sampled<= {max_pairs:,}",
        os.path.join(FIGURES_DIR, "liver_phospho_corr_sample_hist.png"),
    )

    # Summary CSV for poster table
    summary = pd.DataFrame(
        [
            {
                "dataset": "ProteinExpression",
                "group": "Control",
                "rows_total": len(df_prot),
                "rows_kept_minobs": len(prot_control),
                "corr_n": int(r_prot_ctrl.size),
                "corr_mean": float(np.nanmean(r_prot_ctrl)) if r_prot_ctrl.size else np.nan,
                "corr_median": float(np.nanmedian(r_prot_ctrl)) if r_prot_ctrl.size else np.nan,
            },
            {
                "dataset": "ProteinExpression",
                "group": "Sample",
                "rows_total": len(df_prot),
                "rows_kept_minobs": len(prot_sample),
                "corr_n": int(r_prot_samp.size),
                "corr_mean": float(np.nanmean(r_prot_samp)) if r_prot_samp.size else np.nan,
                "corr_median": float(np.nanmedian(r_prot_samp)) if r_prot_samp.size else np.nan,
            },
            {
                "dataset": "Phosphorylation",
                "group": "Control",
                "rows_total": len(df_ph),
                "rows_kept_minobs": len(ph_control),
                "corr_n": int(r_ph_ctrl.size),
                "corr_mean": float(np.nanmean(r_ph_ctrl)) if r_ph_ctrl.size else np.nan,
                "corr_median": float(np.nanmedian(r_ph_ctrl)) if r_ph_ctrl.size else np.nan,
            },
            {
                "dataset": "Phosphorylation",
                "group": "Sample",
                "rows_total": len(df_ph),
                "rows_kept_minobs": len(ph_sample),
                "corr_n": int(r_ph_samp.size),
                "corr_mean": float(np.nanmean(r_ph_samp)) if r_ph_samp.size else np.nan,
                "corr_median": float(np.nanmedian(r_ph_samp)) if r_ph_samp.size else np.nan,
            },
        ]
    )

    summary_path = os.path.join(FIGURES_DIR, "liver_correlation_summary.csv")
    summary.to_csv(summary_path, index=False)

    print("Wrote liver correlation outputs to figures/:")
    print("  liver_protein_corr_control_hist.png")
    print("  liver_protein_corr_sample_hist.png")
    print("  liver_phospho_corr_control_hist.png")
    print("  liver_phospho_corr_sample_hist.png")
    print("  liver_correlation_summary.csv")
