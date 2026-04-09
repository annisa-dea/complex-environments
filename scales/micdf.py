"""
micdf.py
--------
MICDF (Mutual Information Cumulative Distribution Function) analysis.

Pipeline
--------
1. run_svd_folds()   — stratified k-fold split → SVD on log-norm HVG matrix
                       per fold; saves U, s, Vt and HVG gene names to disk
2. compute_micdf()   — loads SVD outputs; computes MI per spectral window per
                       condition; applies null subtraction; returns mean ± STE
3. plot_micdf()      — plots MICDF curves with shaded error bands
4. top_genes_per_pc() — returns the top-N genes by absolute loading for a
                        given list of PCs, based on the Vt matrix from SVD

References
----------
Zaydman et al. eLife 2022 — MICDF framework
Doran et al. 2025         — evolved yeast application

Mathematical note on SVD choice
--------------------------------
PCA is equivalent to SVD of the (centered/scaled) data matrix.
Running SVD on adata.obsm['X_pca'] would be a double decomposition.

Ground truth (SCALES_Kfold.ipynb): SVD is run on adata_fold.raw.X, which
contains log-normalized, HVG-subset, UNSCALED counts. adata.raw is frozen
before sc.pp.scale() is called in load_data_raw.ipynb, so raw.X deliberately
preserves the natural variance structure of the transcriptome. z-scoring
before SVD would equalize variance across genes and distort the spectral
decomposition used by MICDF. SCALE_BEFORE_SVD is therefore False by default.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scipy.sparse
from pathlib import Path
from sklearn.model_selection import StratifiedKFold
from scipy.spatial.distance import squareform

import spipy as sp


# ══════════════════════════════════════════════════════════════════════════════
# ── HARD-CODED SETTINGS ───────────────────────────────────────────────────────
# These define which environmental axes MI is computed for.
# Edit here if you need to add or remove conditions.
# 'shuffled' is always included as the null reference and is NOT plotted.

ENV_LIST = ["rich_media", "30c", "dglucose", "ros", "osmolytes"]
SHUFFLED_LABEL = "shuffled"

# Color map for plotting — add entries here if ENV_LIST is extended
ENV_COLORS = {
    "dglucose" : "#2ca02c",   # green
    "rich_media": "#ffdd00",  # yellow
    "ros"       : "#9467bd",  # purple
    "30c"       : "#d62728",  # red
    "osmolytes" : "#1f77b4",  # blue
}

# SVD / MI parameters
N_FOLDS       = 6
N_WINDOW      = 3     # number of PCs per spectral window
N_WINDOWS     = 40    # number of windows to compute — matches notebook range(40)
RANDOM_STATE  = 12345

# False = pull log-norm HVG counts from adata.raw.X
# True  = pull ALL genes from adata.raw.X (DEFAULT — matches notebook ground truth)
SCALE_BEFORE_SVD = False

# True  = SVD on all genes in raw.X (DEFAULT — matches notebook Cell 9)
# False = SVD on HVG subset only (diagnostic)
FULL_GENE_SVD = True
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. SVD PER FOLD ───────────────────────────────────────────────────────────

def run_svd_folds(
    adata,
    output_dir: str | Path,
) -> Path:
    """
    Split cells into N_FOLDS stratified folds (by condition), then run
    truncated SVD on the log-normalized HVG matrix for each fold.

    Saves per-fold .npz files (U, s, Vt) and a gene names .csv to output_dir.

    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData. Must have:
          - adata.raw.X  : log-normalized sparse matrix (all genes)
          - adata.var['highly_variable'] : HVG boolean mask
          - adata.obs['condition'] : condition labels
    output_dir : str or Path
        Directory to write SVD outputs.

    Returns
    -------
    Path
        Path to output_dir (for chaining into compute_micdf).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Extract the matrix for SVD ────────────────────────────────────────────
    # After preprocess.py, adata.X is the z-scored HVG matrix (all genes
    # already subset to HVGs, sc.pp.scale applied). gene names come from
    # adata.var_names directly.
    #
    # If SCALE_BEFORE_SVD=False, we fall back to pulling log-norm counts from
    # adata.raw.X and subsetting to HVGs manually — for unscaled SVD only.

    if SCALE_BEFORE_SVD:
        # adata.X is already z-scored HVGs — use directly
        gene_names = adata.var_names.tolist()
        X = adata.X
        if scipy.sparse.issparse(X):
            hvg_matrix_full = np.asarray(X.todense())
        else:
            hvg_matrix_full = np.array(X)
        print("[svd] Using z-scored adata.X (scaled HVGs from preprocess.py)")
    elif FULL_GENE_SVD:
        # Use all genes from raw.X — diagnostic mode, no HVG filtering
        gene_names = adata.raw.var_names.tolist()
        raw_X = adata.raw.X
        if scipy.sparse.issparse(raw_X):
            hvg_matrix_full = np.asarray(raw_X.todense())
        else:
            hvg_matrix_full = np.array(raw_X)
        print(f"[svd] Using full-gene log-norm matrix ({len(gene_names)} genes)")
    else:
        # Pull log-norm counts from raw and subset to HVGs.
        # adata.raw.var does not carry 'highly_variable' (frozen before HVG step).
        # Instead, reconstruct the mask by matching adata.var_names (the HVG list)
        # against adata.raw.var_names (all genes).
        hvg_set = set(adata.var_names)
        hvg_mask = np.array([g in hvg_set for g in adata.raw.var_names])
        gene_names = adata.var_names.tolist()
        raw_X = adata.raw.X
        if scipy.sparse.issparse(raw_X):
            hvg_matrix_full = np.asarray(raw_X[:, hvg_mask].todense())
        else:
            hvg_matrix_full = np.array(raw_X[:, hvg_mask])
        print("[svd] Using unscaled log-norm HVG matrix from adata.raw.X")

    # ── Stratified k-fold split ───────────────────────────────────────────────
    # Ground truth (notebook Cell 7): second argument to skf.split() is
    # np.zeros(n_samples) — i.e. no real stratification, pure random split.
    # This matches the notebook exactly; do NOT stratify by condition here.
    conditions = adata.obs["condition"].values
    skf = StratifiedKFold(n_splits=N_FOLDS, random_state=RANDOM_STATE, shuffle=True)
    fold_ids = np.zeros(adata.n_obs, dtype=int)
    for fold_i, (_, test_idx) in enumerate(skf.split(np.zeros(adata.n_obs), np.zeros(adata.n_obs))):
        fold_ids[test_idx] = fold_i
    adata.obs["fold_id"] = fold_ids

    # ── SVD per fold ─────────────────────────────────────────────────────────
    for fold_i in range(N_FOLDS):
        fold_mask = fold_ids == fold_i
        X_fold = hvg_matrix_full[fold_mask, :]

        print(f"[svd] Fold {fold_i}: {X_fold.shape[0]} cells × {X_fold.shape[1]} HVGs")
        U, s, Vt = np.linalg.svd(X_fold, full_matrices=False)

        np.savez_compressed(
            output_dir / f"svd_fold{fold_i}.npz",
            U=U, s=s, Vt=Vt,
        )

    # Save gene names and fold assignments
    pd.Series(gene_names).to_csv(output_dir / "hvg_gene_names.csv",
                                 index=False, header=["gene"])
    adata.obs[["condition", "fold_id"]].to_csv(output_dir / "fold_assignments.csv")

    print(f"[svd] SVD outputs saved → {output_dir}")
    return output_dir



# ── 2. COMPUTE MICDF ──────────────────────────────────────────────────────────

def compute_micdf(
    adata,
    svd_dir: str | Path,
    output_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load per-fold SVD outputs, compute MI per spectral window, normalize to
    CDF, and average across folds.

    MI formula (SCALES_Kfold.ipynb ground truth)
    ---------------------------------------------
        raw_MI[window, env]  = empiricalMI(U_window, label_env)
        cumulative_MI        = prepend_zero_row(cumsum(raw_MI, axis=0))
        MICDF[env]           = cumulative_MI[env] / cumulative_MI[env].iloc[-1]

    'shuffled' is treated as just another environment — it is computed and
    plotted alongside the real conditions. Its curve serves as a visual
    reference for the ambient MI floor; no subtraction is applied.

    Parameters
    ----------
    adata : AnnData
        Same preprocessed AnnData used in run_svd_folds().
        obs must contain 'fold_id', 'condition', all ENV_LIST columns,
        and SHUFFLED_LABEL.
    svd_dir : str or Path
        Directory containing svd_fold*.npz and fold_assignments.csv.
    output_dir : str or Path or None
        If provided, saves mean and STE CSVs.

    Returns
    -------
    mean_micdf : pd.DataFrame  (N_WINDOWS+1 × len(ENV_LIST) + 1)
    ste_micdf  : pd.DataFrame  (N_WINDOWS+1 × len(ENV_LIST) + 1)
        Columns include all ENV_LIST entries plus SHUFFLED_LABEL.
    """
    svd_dir = Path(svd_dir)
    all_envs = ENV_LIST + [SHUFFLED_LABEL]

    # ── Verify obs columns ────────────────────────────────────────────────────
    missing_cols = [e for e in all_envs if e not in adata.obs.columns]
    if missing_cols:
        raise ValueError(
            f"adata.obs is missing columns required for MI: {missing_cols}\n"
            f"These must match metadata column names from ConditionsMatrixCombined.csv"
        )

    micdf_per_fold = []

    for fold_i in range(N_FOLDS):
        svd_path = svd_dir / f"svd_fold{fold_i}.npz"
        svd = np.load(svd_path)
        U, s = svd["U"], svd["s"]

        fold_mask = adata.obs["fold_id"].values == fold_i
        adata_fold = adata[fold_mask]

        # ── Raw MI: shape (N_WINDOWS × len(all_envs)) ────────────────────────
        raw_mi = pd.DataFrame(
            0.0,
            index=range(N_WINDOWS),
            columns=all_envs,
        )

        for w in range(N_WINDOWS):
            u_corr = sp.calc_spcorr_mtx(U, s, range(w, w + N_WINDOW))
            u_corr_sq = squareform(u_corr, checks=False)

            for env in all_envs:
                env_labels = adata_fold.obs[env].values
                env_corr = squareform(
                    env_labels[:, None] == env_labels, checks=False
                )
                raw_mi.loc[w, env] = sp.empiricalMI_masked(u_corr_sq, env_corr)

        # ── Normalize to CDF ──────────────────────────────────────────────────
        # Ground truth (notebook Cell 27):
        #   cumulative_mi = cumulative_mi_unscaled.loc[0:40] / cumulative_mi_unscaled.loc[40]
        # Divide by the value at row index N_WINDOWS specifically (not iloc[-1],
        # which would be the zero-prepended row 0 if indexing goes wrong).
        zero_row = pd.DataFrame(
            {env: [0.0] for env in all_envs}
        )
        cumsum = pd.concat(
            [zero_row, raw_mi.cumsum(axis=0)], axis=0
        ).reset_index(drop=True)

        # Normalize by the value at row N_WINDOWS (index = N_WINDOWS after prepending zero)
        total = cumsum.loc[N_WINDOWS].copy()
        total[total == 0] = 1   # avoid divide-by-zero on flat columns
        micdf_fold = cumsum.loc[0:N_WINDOWS] / total

        micdf_per_fold.append(micdf_fold)
        print(f"[micdf] Fold {fold_i} complete")

    # ── Average across folds ──────────────────────────────────────────────────
    concat_all = pd.concat(micdf_per_fold)
    mean_micdf = concat_all.groupby(concat_all.index).mean()
    std_micdf  = concat_all.groupby(concat_all.index).std()
    ste_micdf  = std_micdf / N_FOLDS   # notebook Cell 30: ste = std / 6

    # ── Save ─────────────────────────────────────────────────────────────────
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        mean_micdf.to_csv(output_dir / "micdf_mean.csv")
        ste_micdf.to_csv(output_dir / "micdf_ste.csv")
        print(f"[micdf] CSVs saved → {output_dir}")

    return mean_micdf, ste_micdf


# ── 3. PLOT MICDF ─────────────────────────────────────────────────────────────

def plot_micdf(
    mean_micdf: pd.DataFrame,
    ste_micdf: pd.DataFrame,
    output_path: str | Path | None = None,
    title: str = "MICDF",
) -> plt.Figure:
    """
    Plot MICDF curves with shaded STE error bands.

    Parameters
    ----------
    mean_micdf : pd.DataFrame
    ste_micdf  : pd.DataFrame
    output_path : str or Path or None
        If provided, saves figure as PDF.
    title : str
        Plot title.

    Returns
    -------
    matplotlib Figure
    """
    fig, ax = plt.subplots(figsize=(6, 4))

    # Plot real environments + shuffled — shuffled is gray and dashed
    envs_to_plot = ENV_LIST + [SHUFFLED_LABEL]
    for env in envs_to_plot:
        if env not in mean_micdf.columns:
            print(f"[plot] WARNING: '{env}' not in mean_micdf columns, skipping")
            continue
        if env == SHUFFLED_LABEL:
            color, lw, ls = "#aaaaaa", 1.0, "--"
        else:
            color, lw, ls = ENV_COLORS.get(env, "#888888"), 1.5, "-"
        x = mean_micdf.index.values
        y = mean_micdf[env].values
        lower = y - ste_micdf[env].values
        upper = y + ste_micdf[env].values

        ax.plot(x, y, color=color, linewidth=lw, linestyle=ls, label=env)
        ax.fill_between(x, lower, upper, color=color, alpha=0.2)

    ax.set_xlabel("Spectral window (by 3 PCs)")
    ax.set_ylabel("Mutual Information CDF")
    ax.set_xlim(0, N_WINDOWS)
    ax.set_ylim(0, 1)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(title, fontsize=10)
    plt.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=150)
        print(f"[plot] Figure saved → {output_path}")

    return fig


# ── 4. TOP GENES PER PC ───────────────────────────────────────────────────────

def top_genes_per_pc(
    svd_dir: str | Path,
    pcs: list[int],
    n: int = 10,
    fold: int = 0,
) -> pd.DataFrame:
    """
    Return the top-N genes by absolute loading magnitude for each requested PC.

    Loadings are taken from the Vt matrix (right singular vectors) of a
    single SVD fold. Vt has shape (min(cells, HVGs) × HVGs); row i of Vt
    is the gene-space loading vector for PC i.

    Parameters
    ----------
    svd_dir : str or Path
        Directory containing svd_fold*.npz and hvg_gene_names.csv.
    pcs : list of int
        Zero-indexed PC numbers to report (e.g. [0, 1, 2, 3, 4] for top 5).
    n : int
        Number of top genes to return per PC. Default 10.
    fold : int
        Which fold's SVD to use. Default 0.

    Returns
    -------
    pd.DataFrame
        Columns: ['PC', 'rank', 'gene', 'loading']
        Sorted by PC then rank.
    """
    svd_dir = Path(svd_dir)
    svd = np.load(svd_dir / f"svd_fold{fold}.npz")
    Vt = svd["Vt"]   # shape: (n_components × n_HVGs)

    gene_names = pd.read_csv(svd_dir / "hvg_gene_names.csv")["gene"].tolist()

    if len(gene_names) != Vt.shape[1]:
        raise ValueError(
            f"Gene name count ({len(gene_names)}) != Vt columns ({Vt.shape[1]})"
        )

    rows = []
    for pc in pcs:
        if pc >= Vt.shape[0]:
            print(f"[top_genes] WARNING: PC {pc} exceeds Vt rows ({Vt.shape[0]}), skipping")
            continue
        loadings = Vt[pc, :]
        top_idx = np.argsort(np.abs(loadings))[::-1][:n]
        for rank, idx in enumerate(top_idx, start=1):
            rows.append({
                "PC"     : pc,
                "rank"   : rank,
                "gene"   : gene_names[idx],
                "loading": float(loadings[idx]),
            })

    result = pd.DataFrame(rows)
    return result
