"""
merge.py
--------
Merges per-condition CellRanger .h5 files into a single AnnData object,
annotating each cell with its per-condition metadata from ConditionsMatrixCombined.csv.

Usage
-----
    from scales.merge import merge_h5

    adata = merge_h5(
        h5_dir   = "/path/to/h5_files",
        metadata_path = "/path/to/ConditionsMatrixCombined.csv",
        ancestry = "anc",   # "anc", "evo", or "combined"
        output_path = "/path/to/merged.h5ad",   # optional
    )

Ancestry filter
---------------
    "anc"      : ancestry == 0  (20 ancestral conditions)
    "evo"      : ancestry == 1  (18 evolved conditions)
    "combined" : all 38 conditions
"""

import os
import numpy as np
import anndata as ad
import scanpy as sc
import pandas as pd
from pathlib import Path


# ── RANDOM SEED ──────────────────────────────────────────────────────────────
RANDOM_SEED = 12345


# ── ANCESTRY FILTER MAP ───────────────────────────────────────────────────────
# Maps the ancestry argument to the integer value used in ConditionsMatrixCombined.csv
_ANCESTRY_VALUES = {
    "anc": 0,
    "evo": 1,
    "combined": None,   # no filter applied
}


def merge_h5(
    h5_dir: str | Path,
    metadata_path: str | Path,
    ancestry: str = "combined",
    output_path: str | Path | None = None,
) -> ad.AnnData:
    """
    Merge per-condition CellRanger filtered_feature_bc_matrix.h5 files into
    a single AnnData, with per-cell metadata columns from the conditions matrix.

    Parameters
    ----------
    h5_dir : str or Path
        Directory containing one .h5 file per condition.
        File names must match the condition index in the metadata CSV
        (e.g. "anc_ypd.h5" matches index "anc_ypd").
    metadata_path : str or Path
        Path to ConditionsMatrixCombined.csv.
    ancestry : {"anc", "evo", "combined"}
        Which subset of conditions to load.
        "anc"      → ancestry == 0
        "evo"      → ancestry == 1
        "combined" → all conditions
    output_path : str or Path or None
        If provided, saves the merged AnnData to this path as .h5ad.

    Returns
    -------
    AnnData
        Merged object with obs columns from the metadata CSV and a
        top-level "condition" label added by ad.concat.
    """
    h5_dir = Path(h5_dir)
    metadata_path = Path(metadata_path)

    if ancestry not in _ANCESTRY_VALUES:
        raise ValueError(
            f"ancestry must be one of {list(_ANCESTRY_VALUES.keys())}, got {ancestry!r}"
        )

    # ── Load and filter metadata ──────────────────────────────────────────────
    metadata = pd.read_csv(metadata_path, index_col=0)

    ancestry_val = _ANCESTRY_VALUES[ancestry]
    if ancestry_val is not None:
        metadata = metadata[metadata["ancestry"] == ancestry_val]

    print(f"[merge] ancestry='{ancestry}' → {len(metadata)} conditions selected")

    # ── Load each condition's h5 ──────────────────────────────────────────────
    adatas = {}
    missing = []

    for condition, row in metadata.iterrows():
        h5_path = h5_dir / f"{condition}.h5"
        if not h5_path.exists():
            missing.append(str(h5_path))
            continue

        adata_cond = sc.read_10x_h5(h5_path)

        # Prefix barcodes with condition name to guarantee uniqueness after concat
        adata_cond.obs_names = [f"{condition}_{bc}" for bc in adata_cond.obs_names]
        adata_cond.var_names_make_unique()

        # Annotate each cell with all metadata columns for this condition
        for col, val in row.items():
            adata_cond.obs[col] = val

        adatas[condition] = adata_cond

    if missing:
        print(f"[merge] WARNING: {len(missing)} .h5 file(s) not found:")
        for p in missing:
            print(f"  {p}")

    if not adatas:
        raise FileNotFoundError(
            f"No .h5 files matched the selected conditions in {h5_dir}"
        )

    # ── Concatenate ──────────────────────────────────────────────────────────
    merged = ad.concat(adatas, label="condition")

    n_dup = merged.obs_names.duplicated().sum()
    if n_dup > 0:
        print(f"[merge] WARNING: {n_dup} duplicate obs_names detected after concat")

    # Add shuffled null label: random binary per cell, used as MI null in MICDF.
    # Generated here so it travels with the AnnData through all downstream steps.
    # Not in ConditionsMatrixCombined.csv -- it is per-cell, not per-condition.
    rng = np.random.default_rng(RANDOM_SEED)
    merged.obs["shuffled"] = rng.integers(0, 2, size=merged.n_obs)

    print(f"[merge] Merged AnnData: {merged.n_obs} cells x {merged.n_vars} genes")
    print(f"[merge] Conditions loaded: {merged.obs['condition'].nunique()}")

    # ── Optionally save ───────────────────────────────────────────────────────
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        merged.write_h5ad(output_path)
        print(f"[merge] Saved → {output_path}")

    return merged
