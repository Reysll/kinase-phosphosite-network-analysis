Kinase–Phosphosite Network Analysis

Comparative Systems-Level Modeling of Signaling in Generic and Liver Cancer Contexts

Overview

This repository contains the analysis and visualization framework built on top of a heterogeneous kinase–phosphosite signaling network.
The project consists of two major network instances:

1.-Generic Network
Constructed from curated signaling databases:
Kinase–substrate relationships
Phosphatase–substrate relationships
Protein–protein interactions
Pathway-based co-membership
Site-level co-evolution relationships

2.- Liver Network
Built by integrating:
The Generic Network
Liver cancer proteomics data (protein expression)
Liver phosphoproteomics data (phosphosite abundance)
Correlation-derived edges (control vs cancer)

The goal of this repository is to:

Quantify structural properties of signaling networks
Compare global topology between generic and disease-specific contexts
Visualize hub structure and local signaling neighborhoods
Characterize correlation shifts in protein and phosphosite behavior between control and cancer samples
Produce publication-ready figures for poster or manuscript use

Scientific Motivation

Cellular signaling networks are complex, highly interconnected systems.

Kinases and phosphatases regulate protein function via phosphorylation at specific sites. These regulatory events:

Control cell cycle progression
Drive oncogenic signaling
Alter transcriptional programs
Modify protein–protein interactions

Traditional pathway analysis treats signaling as modular and discrete.

In contrast, this project models signaling as a multi-layer network, integrating:

Protein-level relationships
Site-level regulation
Physical interaction data
Pathway co-membership
Disease-specific correlation structure

This systems-level view allows us to:

Identify hub proteins
Quantify network density and connectivity
Compare topology between healthy and cancer states
Visualize local rewiring around critical oncogenic proteins



Network Architecture
Node Types

All protein-level entities are unified as:

PROTEIN:<gene>

There are no separate KINASE or PHOSPHATASE node types.

Site-level nodes are represented as:

SITE:<gene>-<residue>

Example:
SITE:MAPK1-T185


This unification ensures:

No duplication of biological identity
All connections to a gene are centralized
Graph semantics remain clean and interpretable

Edge Types

The network includes the following relations:

Relation	                Meaning
has_site	                Protein contains phosphosite
phosphorylates	            Protein phosphorylates site
dephosphorylates	        Protein removes phosphorylation
ppi_high_confidence	        Protein–protein interaction
protein_same_pathway	    Proteins share pathway membership
site_same_pathway	        Sites share pathway membership
site_coevolution	        Co-evolution evidence between sites


The Liver Network additionally includes:

Relation	                Meaning
protein_corr_control	    Protein correlation (control)
protein_corr_sample	        Protein correlation (cancer)
phospho_corr_control	    Phosphosite correlation (control)
phospho_corr_sample	           Phosphosite correlation (cancer)


Repository Structure
data/
    nodes.csv.gz
    edges.csv.gz
    liver/
        Liver_network_nodes.csv.gz
        Liver_network_edges.csv.gz
        LiverCancer_ProtExp_Phospho_casecntrl.xlsx

figures/
    network plots
    degree distributions
    hub summaries
    subgraphs

figures/figures_liver/
    liver-specific plots
    correlation histograms
    heatmaps
    comparison summaries

src/
    network_stats.py
    liver_io.py
    liver_correlations.py
    liver_plots.py
    compare_networks.py
    subgraph_plot.py
    liver_heatmap.py

Installation
1. Clone repository
git clone https://github.com/Reysll/kinase-phosphosite-network-analysis.git
cd kinase-phosphosite-network-analysis

2. Create virtual environment

Windows:

python -m venv .venv
.\.venv\Scripts\activate


Mac/Linux:

python3 -m venv .venv
source .venv/bin/activate

3. Install dependencies
pip install -r requirements.txt

Running Analyses
1. Basic Network Statistics
python -m src.network_stats


Outputs:

node counts
edge counts
density
connected components
top hubs
degree distribution plots

Saved to:
figures/

2. Liver Data Verification
python -m src.check_liver

Confirms:
Number of proteins
Number of phosphosites

3. Liver Correlation Analysis
python -m src.liver_plots


Generates:
Protein correlation histogram (control)
Protein correlation histogram (cancer)
Phosphosite correlation histogram (control)
Phosphosite correlation histogram (cancer)

Saved to:
figures/

Filtering rules applied:

Proteins/sites must be detected in ≥6 samples per group
Correlation only computed if ≥6 shared non-missing values
Ensures statistical reliability
Reduces noise

4. Generic vs Liver Network Comparison
python -m src.compare_networks


Outputs:

Node count differences
Edge count differences
Edge type expansion
Summary CSV

Saved to:
figures/figures_liver/

5. Subgraph Visualization (Ego Networks)

Example:
python -m src.subgraph_plot --seed PROTEIN:EGFR --k 1

Liver version:
python -m src.subgraph_plot --seed PROTEIN:EGFR --k 1 --liver

Filtered mechanistic edges:
python -m src.subgraph_plot --seed PROTEIN:EGFR --k 1 --relations phosphorylates,dephosphorylates,has_site,ppi_high_confidence --max-nodes 250 --max-edges 6000


This creates readable signaling neighborhoods.

6. Liver Heatmap
python -m src.liver_heatmap


Generates:
Correlation heatmaps
Control vs cancer visual comparison

Why Visualize?
Network visualization is meaningful because:
Degree distribution reveals hub dominance
Heavy-tailed distributions suggest scale-free topology
Subgraphs expose local signaling architecture
Correlation histograms show distribution shifts between conditions
Heatmaps reveal structured co-regulation patterns
Visualization translates raw network statistics into biological insight.

Key Results (Current Version)
Generic Network

~13,318 nodes
~1.09M edges
Single giant connected component
Heavy-tailed degree distribution
Hub proteins include:
CDK1
MAPK1
AKT1
PRKCA
MYC

Liver Network

~16,605 nodes
~1.29M edges
+3,287 nodes vs generic
+201,704 edges vs generic
+3 edge types from correlation integration

Correlation Insights

Protein expression:

Control distribution differs from cancer
Broader spread observed in tumor samples

Phosphosite abundance:

Increased variability in cancer
Suggests signaling rewiring

Biological Interpretation

The liver network demonstrates:

Increased site-level complexity
Expanded connectivity
Altered correlation structure
Potential rewiring around hub oncogenes

This supports the hypothesis that:

Cancer alters not just expression levels,
but the topology of signaling relationships.

Reproducibility Notes

All datasets included
All plots are generated from raw data
No manual intervention required
Outputs are deterministic
All filtering thresholds explicitly documented

Future Extensions

Node2Vec embeddings
Kinase–substrate prediction via Random Forest
5-fold cross validation
AUC, sensitivity, specificity reporting

Citation
If using this repository:
Please cite the associated abstract and network methodology.

Maintainer

Salvador
GitHub: https://github.com/Reysll