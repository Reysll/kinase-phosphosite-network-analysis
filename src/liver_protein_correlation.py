import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from itertools import combinations

DATA_PATH = Path("data/LiverCancer_ProtExp_Phospho_casecntrl.xlsx")
FIGURES_DIR = Path("figures")
FIGURES_DIR.mkdir(exist_ok=True)

MIN_SAMPLES = 6


def load_protein_expression():
    df = pd.read_excel(DATA_PATH, sheet_name="Protein Expression", header=1)

    control_cols = [c for c in df.columns if "Control" in str(c)]
    case_cols = [c for c in df.columns if "Sample" in str(c)]

    return df, control_cols, case_cols


def filter_proteins(df, cols):
    counts = df[cols].notna().sum(axis=1)
    return df.loc[counts >= MIN_SAMPLES].copy()


def compute_correlations(df, cols):
    correlations = []
    values = df[cols]

    for (_, row_i), (_, row_j) in combinations(values.iterrows(), 2):
        x = row_i.values.astype(float)
        y = row_j.values.astype(float)

        mask = ~np.isnan(x) & ~np.isnan(y)
        if mask.sum() >= MIN_SAMPLES:
            r = np.corrcoef(x[mask], y[mask])[0, 1]
            if not np.isnan(r):
                correlations.append(r)

    return np.array(correlations)


def plot_hist(control_corr, case_corr):
    plt.figure(figsize=(8,6))
    plt.hist(control_corr, bins=50, alpha=0.6, density=True, label="Control")
    plt.hist(case_corr, bins=50, alpha=0.6, density=True, label="Cancer")
    plt.xlabel("Pearson correlation")
    plt.ylabel("Density")
    plt.title("Protein Expression Correlation Distribution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "protein_corr_control_vs_case.png", dpi=300)
    plt.close()


def main():
    df, control_cols, case_cols = load_protein_expression()

    df_control = filter_proteins(df, control_cols)
    df_case = filter_proteins(df, case_cols)

    print("Proteins in control:", len(df_control))
    print("Proteins in case:", len(df_case))

    control_corr = compute_correlations(df_control, control_cols)
    case_corr = compute_correlations(df_case, case_cols)

    print("Control mean:", np.mean(control_corr))
    print("Case mean:", np.mean(case_corr))

    plot_hist(control_corr, case_corr)

    print("Figure saved.")


if __name__ == "__main__":
    main()
