# Multi-Omics Analysis of Syndromic and Idiopathic Autism

## Overview

This repository contains data and scripts used for an integrative multi-omics analysis and deep learning–based drug repurposing study in autism. The workflow is designed to systematically compare syndromic and idiopathic autism at the molecular, regulatory, and therapeutic levels.

---

## Data

The repository includes curated and processed datasets used in the study:

* Syndromic autism gene lists
* Idiopathic autism gene lists
* Hub gene datasets
* miRNA–gene interaction data
* Drug repurposing input and output data

---

## Methods

The analysis pipeline integrates multiple computational and bioinformatics approaches:

* **Gene Curation:** High-confidence autism-associated genes retrieved from SFARI and AutismKB
* **Functional Enrichment:** GO and KEGG pathway analysis
* **PPI Network Analysis:** Network construction using STRING and visualization in Cytoscape
* **Hub Gene Identification:** Centrality-based prioritization of key genes
* **Transcriptomic Validation:** Expression profiling using GTEx, Human Protein Atlas, and Allen Brain Atlas
* **miRNA Regulatory Analysis:** Construction of brain-specific miRNA–hub gene interaction networks
* **Drug Repurposing:** Deep learning–based drug–target interaction prediction using pretrained DeepPurpose models
* **BBB Filtering:** In silico prediction of blood–brain barrier permeability

---

## Repository Structure

```
transcriptomic_validation/    # GTEx, HPA, and scRNA-seq analysis scripts
mirna_analysis/               # miRNA–gene regulatory network scripts
functional_enrichment/        # GO and KEGG enrichment analysis
network_analysis/             # PPI network and hub gene identification
drug_repurposing/             # DeepPurpose-based drug prediction
data/                         # Input datasets
results/                      # Output files and processed results
figures/                      # Figures and visualizations
```

---

## Applicability to Autism Subtypes

All scripts in this repository are designed to be applied to both syndromic and idiopathic autism gene sets. The analytical workflow remains consistent across subtypes, with differences arising only from input data. This enables direct comparison of molecular mechanisms, regulatory networks, and therapeutic targets between autism subtypes.

---

## Requirements

Python 3.13

Required packages:

* numpy
* pandas
* scipy
* statsmodels
* matplotlib
* seaborn
* networkx

(Additional dependencies may be required for specific analyses such as DeepPurpose.)

---

## Reproducibility

All analyses were performed using open-source tools. Scripts are organized by analysis module and can be executed independently. The repository provides sufficient information to reproduce the results presented in the study.

---

## Author

Rudhra Ondippili

---

## Contact

For queries or collaboration:
[ashauthra@gmail.com]

---

## License

This repository is intended for academic and research use. A suitable open-source license (e.g., MIT) can be added if required.

---

## Citation

If you use this repository, please cite the associated manuscript (details to be added upon publication).
