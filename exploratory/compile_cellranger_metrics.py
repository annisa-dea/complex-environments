"""
compile_cellranger_metrics.py
------------------------------
Walks metrics_summary/<ancestry>/<condition>/outs/metrics_summary.csv
and compiles all samples into a single CSV, with condition names prefixed
by ancestry (anc_ or evo_).

Usage
-----
    python compile_cellranger_metrics.py \
        --metrics_dir /path/to/metrics_summary \
        --output /path/to/cellranger_metrics_all.csv
"""

import argparse
import pandas as pd
from pathlib import Path


def compile_metrics(metrics_dir: Path, output_path: Path) -> pd.DataFrame:
    rows = []
    missing = []

    for ancestry_dir in sorted(metrics_dir.iterdir()):
        if not ancestry_dir.is_dir():
            continue
        ancestry = ancestry_dir.name  # "anc" or "evo"

        for condition_dir in sorted(ancestry_dir.iterdir()):
            if not condition_dir.is_dir():
                continue
            condition = condition_dir.name  # e.g. "ypd", "sdc-kcl"

            csv_path = condition_dir / "outs" / "metrics_summary.csv"
            if not csv_path.exists():
                missing.append(str(csv_path))
                continue

            df = pd.read_csv(csv_path)
            row = df.iloc[0].to_dict()
            row["condition"] = f"{ancestry}_{condition}"
            rows.append(row)

    if missing:
        print(f"WARNING: {len(missing)} metrics_summary.csv file(s) not found:")
        for p in missing:
            print(f"  {p}")

    compiled = pd.DataFrame(rows)

    # Put condition column first
    cols = ["condition"] + [c for c in compiled.columns if c != "condition"]
    compiled = compiled[cols]

    # Sort: all anc_ first, then evo_, alphabetically within each group
    compiled["_sort_key"] = compiled["condition"].apply(
        lambda x: (0 if x.startswith("anc_") else 1, x)
    )
    compiled = compiled.sort_values("_sort_key").drop(columns="_sort_key")
    compiled = compiled.reset_index(drop=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    compiled.to_csv(output_path, index=False)
    print(f"Saved {len(compiled)} conditions -> {output_path}")
    print(f"Columns: {list(compiled.columns)}")

    return compiled


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics_dir", required=True,
                        help="Path to metrics_summary/ folder")
    parser.add_argument("--output", default="cellranger_metrics_all.csv",
                        help="Output CSV path")
    args = parser.parse_args()

    compile_metrics(
        metrics_dir  = Path(args.metrics_dir),
        output_path  = Path(args.output),
    )
