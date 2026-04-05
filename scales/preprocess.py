"""
preprocess.py
-------------
QC filtering and normalization of a merged AnnData.

Applies:
  1. QC metric computation (n_genes, total_counts — direct, fast)
  2. Global cell filters (min/max genes, min UMI)
  3. Gene filter (expressed in at least 1 cell)
  4. Log-normalization
  5. Saves log-norm counts to adata.raw before HVG subsetting
  6. HVG selection
  7. Scaling (z-score, max_value=10)
  8. PCA (for UMAP only; MICDF SVD runs on adata.raw in micdf.py)

Note: mitochondrial filtering is intentionally omitted. Yeast mito genes
(Q0xx) are a small fraction and are not used as a filter criterion in this
pipeline. adata.raw is frozen post-log-norm, pre-scale — matching the
ground-truth pipeline in load_data_raw.ipynb.

Usage
-----
    from scales.preprocess import preprocess

    adata = preprocess(
        adata,
        output_path = "/path/to/processed.h5ad",  # optional
    )
"""

import numpy as np
import scipy.sparse
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path


# ── FILTER THRESHOLDS ─────────────────────────────────────────────────────────
# Adjust here if needed; these apply to both ancestral and evolved datasets.
MIN_GENES  = 200    # hard floor: removes near-empty barcodes CR may have included
MIN_COUNTS = 300    # UMI floor consistent with MIN_GENES
MAX_GENES  = 4000   # doublet guard (yeast transcriptome ~6,000 genes total)

# HVG selection parameters
HVG_MIN_MEAN  = 0.0125
HVG_MAX_MEAN  = 3.0
HVG_MIN_DISP  = 0.5

# PCA components (for UMAP only; not used by MICDF)
N_PCA_COMPS = 50


def preprocess(
    adata: ad.AnnData,
    output_path: str | Path | None = None,
    plot: bool = True,
    plot_dir: str | Path | None = None,
) -> ad.AnnData:
    """
    Run QC filtering, normalization, HVG selection, scaling, and PCA.

    adata.raw is frozen after log-normalization but before HVG subsetting
    and scaling — this is the matrix used by MICDF SVD in micdf.py.

    Parameters
    ----------
    adata : AnnData
        Merged AnnData from merge.merge_h5().
    output_path : str or Path or None
        If provided, saves processed AnnData to this path.
    plot : bool
        If True, generates QC distribution plots.
    plot_dir : str or Path or None
        Directory to save plots. If None and plot=True, saves to cwd.

    Returns
    -------
    AnnData
        Filtered, normalized, HVG-subset, scaled AnnData with PCA embedding.
        adata.raw contains log-normalized counts across all genes (pre-HVG).
    """
    print(f"[preprocess] Input: {adata.n_obs} cells x {adata.n_vars} genes")

    # ── Fast QC metrics ───────────────────────────────────────────────────────
    # Compute directly from the sparse matrix — avoids the overhead of
    # sc.pp.calculate_qc_metrics (percent_top passes, mito annotation, log1p).
    X = adata.X
    if scipy.sparse.issparse(X):
        adata.obs["n_genes"]      = np.asarray((X > 0).sum(axis=1)).ravel()
        adata.obs["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
    else:
        adata.obs["n_genes"]      = (X > 0).sum(axis=1)
        adata.obs["total_counts"] = X.sum(axis=1)

    print(f"[preprocess] n_genes     — min: {adata.obs['n_genes'].min():.0f}  "
          f"median: {adata.obs['n_genes'].median():.0f}  "
          f"max: {adata.obs['n_genes'].max():.0f}")
    print(f"[preprocess] total_counts — min: {adata.obs['total_counts'].min():.0f}  "
          f"median: {adata.obs['total_counts'].median():.0f}  "
          f"max: {adata.obs['total_counts'].max():.0f}")

    # ── Pre-filter plots ─────────────────────────────────────────────────────
    if plot:
        _plot_qc(adata, stage="prefilter", plot_dir=plot_dir)

    # ── Cell filters ─────────────────────────────────────────────────────────
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
    sc.pp.filter_cells(adata, max_genes=MAX_GENES)
    sc.pp.filter_genes(adata, min_cells=1)
    n_after = adata.n_obs

    print(f"[preprocess] Cells before filter : {n_before}")
    print(f"[preprocess] Cells after filter  : {n_after}  "
          f"({100 * (n_before - n_after) / n_before:.1f}% removed)")
    print(f"[preprocess] Genes retained: {adata.n_vars}")

    # ── Per-condition summary ─────────────────────────────────────────────────
    summary = (
        adata.obs
        .groupby("condition")[["n_genes", "total_counts"]]
        .agg(n_cells=("n_genes", "count"),
             median_genes=("n_genes", "median"),
             median_umi=("total_counts", "median"))
        .round(1)
    )
    print("\n[preprocess] Post-filter per-condition summary:")
    print(summary.to_string())

    # ── Normalization ────────────────────────────────────────────────────────
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Freeze log-norm full-gene matrix BEFORE HVG subsetting and scaling.
    # micdf.py reads adata.raw.X for SVD — log-norm, unscaled, all genes.
    adata.raw = adata

    # ── HVG selection ────────────────────────────────────────────────────────
    sc.pp.highly_variable_genes(
        adata,
        min_mean=HVG_MIN_MEAN,
        max_mean=HVG_MAX_MEAN,
        min_disp=HVG_MIN_DISP,
    )
    n_hvg = adata.var["highly_variable"].sum()
    print(f"\n[preprocess] HVGs selected: {n_hvg}")
    adata = adata[:, adata.var.highly_variable].copy()

    # ── Scaling ──────────────────────────────────────────────────────────────
    sc.pp.scale(adata, max_value=10)

    # ── PCA (for UMAP; not used by MICDF) ───────────────────────────────────
    sc.tl.pca(adata, n_comps=N_PCA_COMPS, svd_solver="arpack")

    print(f"\n[preprocess] Final AnnData: {adata.n_obs} cells x {adata.n_vars} HVGs")

    # ── Save ─────────────────────────────────────────────────────────────────
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(output_path)
        print(f"[preprocess] Saved -> {output_path}")

    return adata


# ── INTERNAL: QC plot helper ──────────────────────────────────────────────────

def _plot_qc(
    adata: ad.AnnData,
    stage: str,
    plot_dir: str | Path | None,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    axes[0].hist(adata.obs["n_genes"], bins=100,
                 color="#4C72B0", edgecolor="none")
    axes[0].axvline(MIN_GENES, color="red",    linestyle="--", linewidth=1,
                    label=f"min={MIN_GENES}")
    axes[0].axvline(MAX_GENES, color="orange", linestyle="--", linewidth=1,
                    label=f"max={MAX_GENES}")
    axes[0].set_xlabel("Genes per cell")
    axes[0].set_ylabel("Cells")
    axes[0].legend(fontsize=8)
    axes[0].spines[["top", "right"]].set_visible(False)

    axes[1].hist(adata.obs["total_counts"], bins=100,
                 color="#4C72B0", edgecolor="none")
    axes[1].axvline(MIN_COUNTS, color="red", linestyle="--", linewidth=1,
                    label=f"min={MIN_COUNTS}")
    axes[1].set_xlabel("UMI counts per cell")
    axes[1].set_ylabel("Cells")
    axes[1].legend(fontsize=8)
    axes[1].spines[["top", "right"]].set_visible(False)

    axes[2].scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes"],
        s=1, alpha=0.2, color="#4C72B0", rasterized=True,
    )
    axes[2].axvline(MIN_COUNTS, color="red",    linestyle="--", linewidth=0.8)
    axes[2].axhline(MIN_GENES,  color="red",    linestyle="--", linewidth=0.8)
    axes[2].axhline(MAX_GENES,  color="orange", linestyle="--", linewidth=0.8)
    axes[2].set_xlabel("UMI counts per cell")
    axes[2].set_ylabel("Genes per cell")
    axes[2].spines[["top", "right"]].set_visible(False)

    plt.suptitle(f"QC distributions - {stage}", fontsize=11)
    plt.tight_layout()

    if plot_dir is not None:
        plot_dir = Path(plot_dir)
        plot_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_dir / f"qc_{stage}.pdf", bbox_inches="tight", dpi=150)

    #plt.show()
    plt.close(fig)
