# Physiological architecture and evolutionary origins of cellular adaptability

**Dea et al., 2026** | Pincus & Raman Labs, University of Chicago

This repository contains analysis notebooks, processed data, and metadata 
supporting the findings in our paper. We performed single-cell transcriptional 
profiling of budding yeast (*S. cerevisiae*) across 20 complex environments 
before and after long-term laboratory evolution for increased osmotolerance, 
revealing a hierarchical architecture of cellular adaptability shaped by 
evolutionary history.

## Repository structure

your-paper-repo/
├── data/
│   ├── raw/              # CellRanger output (h5 files, barcodes) — not tracked
│   ├── processed/        # Processed AnnData objects (.h5ad)
│   └── metadata/         # Sample and condition metadata
└── notebooks/
    ├── 01_cellranger.md          # CellRanger pipeline documentation
    ├── 02_cell_calling.ipynb     # Adaptive cell calling
    ├── 03_merge_metadata.ipynb   # Merging h5 files, QC filtering, metadata
    ├── 04_preprocess_anc.ipynb   # Ancestral strain: QC, PCA, UMAP
    ├── 05_preprocess_evo.ipynb   # Evolved strain: QC, PCA, UMAP
    ├── 06_combined_preprocess.ipynb  # Combined SVD analysis
    └── 07_micdf_analysis.ipynb   # SCALES/MICDF: environmental hierarchy

## Data availability

Raw scRNA-seq and bulk RNA-seq data are deposited in NCBI GEO under 
accession [GSEXXXXXX]. Processed count matrices and metadata are available 
in this repository.

## Requirements

See `environment.yml` for full dependencies. Key packages:
- Python 3.x
- scanpy, anndata, numpy, pandas, matplotlib, seaborn, spipy (from Raman Lab repo)

## Installation

conda env create -f environment.yml
conda activate <env-name>

## Contact

Annisa Dea — annisa@uchicago.edu 
David Pincus — pincus@uchicago.edu  
Arjun Raman — araman@bsd.uchicago.edu