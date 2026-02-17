# src/check_liver.py
from __future__ import annotations

from src.liver_io import load_liver_protein_expression, load_liver_phosphorylation

#Sanity check: can we load the data without error? Are there any obvious issues with the data?

def main() -> None:
    df1 = load_liver_protein_expression()
    df2 = load_liver_phosphorylation()
    print(f"{len(df1)} proteins")
    print(f"{len(df2)} phosphosites")

if __name__ == "__main__":
    main()
