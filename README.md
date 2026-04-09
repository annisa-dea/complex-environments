# SCALES Analysis

Reproducible analysis pipeline for the SCALES scRNA-seq paper.  
Goes from a directory of CellRanger `.h5` files to all figures and CSVs in the paper.

## Repository layout

```
scales-analysis/
├── scales/                    # importable package
│   ├── merge.py               # merge_h5()
│   ├── preprocess.py          # preprocess()
│   ├── umap.py                # run_umap()
│   └── micdf.py               # run_svd_folds(), compute_micdf(),
│                              # plot_micdf(), top_genes_per_pc()
├── exploratory/               # original development notebooks (read-only reference)
│   ├── 02_merge_h5.ipynb
│   ├── 03_preprocess_anc.ipynb
│   └── 06_generate_micdf_evo.ipynb
├── data/
│   └── ConditionsMatrixCombined.csv
├── run_pipeline.py            # top-level runner
└── environment.yml
```

## Quickstart

```bash
conda env create -f environment.yml
conda activate scales

# Ancestral dataset
python run_pipeline.py \
    --ancestry anc \
    --h5_dir /path/to/anc_h5 \
    --out_dir results/anc

# Evolved dataset
python run_pipeline.py \
    --ancestry evo \
    --h5_dir /path/to/evo_h5 \
    --out_dir results/evo
```

### Re-running only MICDF (SVD already done)

```bash
python run_pipeline.py \
    --ancestry evo \
    --h5_dir /path/to/evo_h5 \
    --out_dir results/evo \
    --skip_svd
```

## Output layout

```
results/anc/
├── merged.h5ad
├── processed.h5ad
├── qc/
│   └── qc_prefilter.pdf
├── umap/
│   ├── umap_condition.pdf
│   └── umap_coords.csv
├── svd/
│   ├── svd_fold0.npz  ..  svd_fold5.npz
│   ├── hvg_gene_names.csv
│   └── fold_assignments.csv
└── micdf/
    ├── micdf_mean.csv
    ├── micdf_ste.csv
    ├── micdf_plot.pdf
    └── top_genes_pc0-4.csv
```

## Customizing ENV_LIST or colors

Open `scales/micdf.py` and edit the block at the top of the file:

```python
ENV_LIST = ["dglucose", "rich_media", "ros", "30c", "osmolytes"]

ENV_COLORS = {
    "dglucose"  : "#2ca02c",
    "rich_media": "#ffdd00",
    ...
}
```

## Dependencies

See `environment.yml`. Key packages: `scanpy`, `anndata`, `numpy`, `scipy`,
`scikit-learn`, `pandas`, `matplotlib`, `spipy` (aramanlab/spipy).


## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

[![CC BY 4.0](https://licensebuttons.net/l/by/4.0/88x31.png)](https://creativecommons.org/licenses/by/4.0/)
