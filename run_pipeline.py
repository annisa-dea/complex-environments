"""
run_pipeline.py
---------------
Top-level script: runs the full SCALES analysis pipeline from a directory
of CellRanger .h5 files to all figures and CSVs used in the paper.

Usage
-----
    python run_pipeline.py --ancestry anc --h5_dir /path/to/h5 --out_dir /path/to/output
    python run_pipeline.py --ancestry evo --h5_dir /path/to/h5 --out_dir /path/to/output
    python run_pipeline.py --ancestry combined --h5_dir /path/to/h5 --out_dir /path/to/output

Output layout
-------------
    out_dir/
    ├── merged.h5ad
    ├── processed.h5ad
    ├── qc/
    │   └── qc_prefilter.pdf
    ├── umap/
    │   ├── umap_condition.pdf
    │   └── umap_coords.csv
    ├── svd/
    │   ├── svd_fold0.npz .. svd_fold5.npz
    │   ├── hvg_gene_names.csv
    │   └── fold_assignments.csv
    └── micdf/
        ├── micdf_mean.csv
        ├── micdf_ste.csv
        ├── micdf_plot.pdf
        └── top_genes_pc0-4.csv
"""

import warnings #expected output when performing functions on sparse matrices
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", message=".*densified.*")

import argparse
from pathlib import Path

from scales.merge      import merge_h5
from scales.preprocess import preprocess
from scales.umap       import run_umap, plot_scree
from scales.micdf      import (
    run_svd_folds,
    compute_micdf,
    plot_micdf,
    top_genes_per_pc,
)

# Path to the metadata file — update if repo layout changes
METADATA_PATH = Path(__file__).parent / "data" / "metadata" / "ConditionsMatrixCombined.csv"


def main():
    parser = argparse.ArgumentParser(
        description="SCALES scRNA-seq pipeline: h5 file directory to figures and CSVs"
    )
    parser.add_argument(
        "--ancestry",
        choices=["anc", "evo", "combined"],
        required=True,
        help="Which conditions to run: 'anc' (ancestry==0), 'evo' (ancestry==1), or 'combined'",
    )
    parser.add_argument(
        "--h5_dir",
        required=True,
        help="Directory containing per-condition CellRanger .h5 files",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Root output directory for all results",
    )
    parser.add_argument(
        "--skip_umap",
        action="store_true",
        help="Skip UMAP computation (useful if only re-running MICDF)",
    )
    parser.add_argument(
        "--skip_svd",
        action="store_true",
        help="Skip SVD step and load existing svd/ outputs (for re-running MI only)",
    )
    parser.add_argument(
        "--scale_before_svd",
        action="store_true",
        help="Z-score genes before SVD (diagnostic; not used in published analysis)",
    )
    parser.add_argument(
        "--full_gene_svd",
        action="store_true",
        help="Run SVD on all genes instead of HVGs only (diagnostic)",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)

    # Apply diagnostic SVD flags by overriding micdf module constants at runtime
    import scales.micdf as _micdf
    if args.scale_before_svd:
        _micdf.SCALE_BEFORE_SVD = True
        print("[config] SCALE_BEFORE_SVD = True (z-score before SVD)")
    if args.full_gene_svd:
        _micdf.FULL_GENE_SVD = True
        print("[config] FULL_GENE_SVD = True (all genes, not just HVGs)")

    # ── Step 1: Merge ────────────────────────────────────────────────────────
    print("\n═══ Step 1: Merge h5 files ═══")
    adata = merge_h5(
        h5_dir        = args.h5_dir,
        metadata_path = METADATA_PATH,
        ancestry      = args.ancestry,
        output_path   = out_dir / "merged.h5ad",
    )

    # ── Step 2: Preprocess ───────────────────────────────────────────────────
    print("\n═══ Step 2: Preprocess ═══")
    adata = preprocess(
        adata,
        output_path = out_dir / "processed.h5ad",
        plot        = True,
        plot_dir    = out_dir / "qc",
    )

    # ── Step 3: UMAP ─────────────────────────────────────────────────────────
    if not args.skip_umap:
        print("\n═══ Step 3: UMAP ═══")
        adata = run_umap(
            adata,
            output_dir = out_dir / "umap",
            color_by   = ["condition"],
        )
        print("\n═══ Step 3b: Scree plot ═══")
        plot_scree(adata, output_dir=out_dir / "umap")
    else:
        print("\n═══ Step 3: UMAP (skipped) ═══")

    # ── Step 4: SVD per fold ─────────────────────────────────────────────────
    svd_dir = out_dir / "svd"
    if not args.skip_svd:
        print("\n═══ Step 4: SVD per fold ═══")
        run_svd_folds(adata, output_dir=svd_dir)
    else:
        print("\n═══ Step 4: SVD (skipped — loading existing outputs) ═══")

    # ── Step 5: MICDF ────────────────────────────────────────────────────────
    print("\n═══ Step 5: Compute MICDF ═══")
    mean_micdf, ste_micdf = compute_micdf(
        adata,
        svd_dir    = svd_dir,
        output_dir = out_dir / "micdf",
    )

    # ── Step 6: Plot ─────────────────────────────────────────────────────────
    print("\n═══ Step 6: Plot MICDF ═══")
    plot_micdf(
        mean_micdf,
        ste_micdf,
        output_path = out_dir / "micdf" / "micdf_plot.pdf",
        title       = (
            f"MICDF -- {args.ancestry}"
            + (" [z-scored]" if args.scale_before_svd else "")
            + (" [full-gene]" if args.full_gene_svd else "")
        ),
    )

    # ── Step 7: Top genes per PC ─────────────────────────────────────────────
    print("\n═══ Step 7: Top genes per PC ═══")
    top_genes = top_genes_per_pc(
        svd_dir = svd_dir,
        pcs     = [0, 1, 2, 3, 4],
        n       = 10,
        fold    = 0,
    )
    top_genes_path = out_dir / "micdf" / "top_genes_pc0-4.csv"
    top_genes.to_csv(top_genes_path, index=False)
    print(f"[top_genes] Saved → {top_genes_path}")
    print(top_genes.to_string(index=False))

    print(f"\n✓ Pipeline complete. All outputs in: {out_dir}")


if __name__ == "__main__":
    main()
