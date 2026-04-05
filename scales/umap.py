"""
umap.py
-------
Computes UMAP embedding and saves a wide-format CSV of coordinates + metadata.

Usage
-----
    from scales.umap import run_umap

    adata = run_umap(
        adata,
        output_dir = "/path/to/figures",
    )
"""

import anndata as ad
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


# ── UMAP PARAMETERS ──────────────────────────────────────────────────────────
N_NEIGHBORS    = 15
N_PCS_FOR_UMAP = 50    # must be <= N_PCA_COMPS set in preprocess.py
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
