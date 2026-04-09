"""
umap.py
-------
Computes UMAP embedding and scree plot from the standard scanpy PCA.

Both run_umap() and plot_scree() operate on adata.obsm['X_pca'] /
adata.uns['pca']['variance_ratio'] — the scanpy PCA run in preprocess.py
(log-norm → HVG → scale → PCA). This is intentionally decoupled from the
SVD run inside micdf.py, which operates on unscaled log-norm HVGs for the
MI computation. The scree plot is a property of the data matrix and should
match standard scRNA-seq convention for comparability with the literature.

Usage
-----
    from scales.umap import run_umap, plot_scree

    adata = run_umap(adata, output_dir="/path/to/figures")
    fig   = plot_scree(adata, output_dir="/path/to/figures")
"""

import anndata as ad
import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


# ── UMAP / SCREE PARAMETERS ───────────────────────────────────────────────────
N_NEIGHBORS    = 15
N_PCS_FOR_UMAP = 50    # must be <= N_PCA_COMPS set in preprocess.py
N_PCS_SCREE    = 50    # number of PCs shown on the scree plot
RANDOM_SEED    = 12345


def run_umap(
    adata: ad.AnnData,
    output_dir: str | Path | None = None,
    color_by: list[str] | None = None,
) -> ad.AnnData:
    """
    Compute neighbors → UMAP, save colored UMAP plots and a wide-format CSV.

    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData with adata.obsm['X_pca'] present.
    output_dir : str or Path or None
        Directory to write plots and CSV. If None, does not save.
    color_by : list of str or None
        obs columns to use for UMAP coloring. Defaults to ['condition'].

    Returns
    -------
    AnnData
        With adata.obsm['X_umap'] populated.
    """
    if color_by is None:
        color_by = ["condition"]

    print(f"[umap] Computing neighbors (n_neighbors={N_NEIGHBORS}, "
          f"n_pcs={N_PCS_FOR_UMAP})")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS_FOR_UMAP, random_state=RANDOM_SEED)
    sc.tl.umap(adata, random_state=RANDOM_SEED)
    print("[umap] UMAP complete")

    # ── Plot ──────────────────────────────────────────────────────────────────
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        for col in color_by:
            if col not in adata.obs.columns:
                print(f"[umap] WARNING: '{col}' not in adata.obs, skipping plot")
                continue
            fig = sc.pl.umap(adata, color=col, show=False, return_fig=True)
            fig.savefig(output_dir / f"umap_{col}.pdf",
                        bbox_inches="tight", dpi=150)
            plt.close(fig)
            print(f"[umap] Plot saved → umap_{col}.pdf")

           # Full CSV (one row per cell)
        umap_df = pd.DataFrame(
            adata.obsm["X_umap"],
            index=adata.obs_names,
            columns=["UMAP_1", "UMAP_2"],
        )
        umap_df = pd.concat([umap_df, adata.obs.reset_index(drop=True)
                             .set_index(umap_df.index)], axis=1)
        csv_path = output_dir / "umap_coords.csv"
        umap_df.to_csv(csv_path)
        print(f"[umap] CSV saved -> {csv_path}  ({len(umap_df)} rows)")

        # Prism-friendly CSV — diagonal block layout
        # Each condition occupies its own contiguous row block.
        # Within a block: X = UMAP_1, Y_condN = UMAP_2 for that condition,
        # all other Y columns = NaN. Conditions are sorted alphabetically.
        #
        # X        Y_cond1  Y_cond2  Y_cond3
        # x1       y1       NaN      NaN
        # x2       y2       NaN      NaN
        # xN1+1    NaN      y1       NaN
        # ...
        conditions = sorted(umap_df["condition"].unique())
        col_names = [f"Y_{cond}" for cond in conditions]

        blocks = []
        for cond in conditions:
            cond_df = umap_df.loc[umap_df["condition"] == cond, ["UMAP_1", "UMAP_2"]].copy()
            block = pd.DataFrame(float("nan"), index=range(len(cond_df)), columns=["X"] + col_names)
            block["X"] = cond_df["UMAP_1"].values
            block[f"Y_{cond}"] = cond_df["UMAP_2"].values
            blocks.append(block)

        prism_df = pd.concat(blocks, ignore_index=True)
        prism_path = output_dir / "umap_prism.csv"
        prism_df.to_csv(prism_path, index=False)
        print(f"[umap] Prism CSV saved -> {prism_path}  "
              f"({len(conditions)} conditions, {len(prism_df)} rows)")

    return adata


# ── SCREE PLOT ────────────────────────────────────────────────────────────────

def plot_scree(
    adata: ad.AnnData,
    output_dir: str | Path | None = None,
    n_pcs: int = N_PCS_SCREE,
) -> plt.Figure:
    """
    Plot the PCA scree (variance explained per PC) from the standard scanpy
    PCA stored in adata.uns['pca']['variance_ratio'].

    This is intentionally independent of the SVD run in micdf.py. The scree
    characterizes the data matrix and should match standard scRNA-seq
    convention (log-norm → HVG → scale → PCA via sc.tl.pca).

    Parameters
    ----------
    adata : AnnData
        Must have adata.uns['pca']['variance_ratio'] populated by sc.tl.pca().
    output_dir : str or Path or None
        If provided, saves scree_plot.pdf and scree_values.csv here.
    n_pcs : int
        Number of PCs to display. Default N_PCS_SCREE (50).

    Returns
    -------
    matplotlib Figure
    """
    if "pca" not in adata.uns or "variance_ratio" not in adata.uns["pca"]:
        raise ValueError(
            "adata.uns['pca']['variance_ratio'] not found. "
            "Run sc.tl.pca(adata) before calling plot_scree()."
        )

    var_ratio = adata.uns["pca"]["variance_ratio"]
    n_pcs     = min(n_pcs, len(var_ratio))
    pcs       = np.arange(1, n_pcs + 1)
    vr        = var_ratio[:n_pcs]

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))

    # Linear scale
    axes[0].plot(pcs, vr, color="#333333", linewidth=1.2)
    axes[0].set_xlabel("PC")
    axes[0].set_ylabel("Fraction variance explained")
    axes[0].set_xlim(1, n_pcs)
    axes[0].set_ylim(0, None)
    axes[0].spines[["top", "right"]].set_visible(False)

    # Log scale — better resolves the long tail
    axes[1].plot(pcs, vr, color="#333333", linewidth=1.2)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("PC")
    axes[1].set_ylabel("Fraction variance explained (log)")
    axes[1].set_xlim(1, n_pcs)
    axes[1].spines[["top", "right"]].set_visible(False)

    fig.suptitle("Scree plot (scanpy PCA)", fontsize=10)
    plt.tight_layout()

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        fig.savefig(output_dir / "scree_plot.pdf", bbox_inches="tight", dpi=150)
        print(f"[scree] Figure saved → {output_dir / 'scree_plot.pdf'}")

        pd.DataFrame({
            "PC"              : pcs,
            "variance_ratio"  : vr,
            "cumulative_ratio": np.cumsum(vr),
        }).to_csv(output_dir / "scree_values.csv", index=False)
        print(f"[scree] Values saved → {output_dir / 'scree_values.csv'}")

    return fig
