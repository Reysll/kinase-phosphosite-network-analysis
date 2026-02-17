# src/liver_plots.py
from __future__ import annotations

from src.liver_correlations import run_liver_correlation_figures


def main() -> None:
    run_liver_correlation_figures(
        min_obs=6,
        min_common=6,
        max_pairs=250_000,
        seed=7,
    )


if __name__ == "__main__":
    main()
