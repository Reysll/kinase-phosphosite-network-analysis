# Kinase Phosphosite Network Analysis (Poster)

This repository contains poster-ready analytics for the generic kinase–phosphosite network.

## Inputs (not stored in git)
Place these files in `./data/`:
- `nodes.csv.gz`
- `edges.csv.gz`
- (optional) `LiverCancer_ProtExp_Phospho_casecntrl.xlsx`

Excel sheet names expected: ProteinExpression and Phosphorylation (headers start at row 2)





## Setup
```powershell
python -m venv .venv
.\.venv\Scripts\activate
pip install -r requirements.txt
