"""Export RCTD CSV outputs into a cell2location-compatible sp.h5ad."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def normalize_rows(weights: pd.DataFrame) -> pd.DataFrame:
    """Row-normalize weight matrix; zero rows stay zero."""
    arr = weights.to_numpy(dtype=float)
    row_sums = arr.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    return pd.DataFrame(arr / row_sums, index=weights.index, columns=weights.columns)


def dominant_from_weights(weights: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Return dominant cell type and abundance per spot."""
    vals = weights.to_numpy(dtype=float)
    order = np.argsort(-vals, axis=1)
    cols = weights.columns.to_numpy()
    dominant = pd.Series(cols[order[:, 0]], index=weights.index, name="dominant_celltype")
    abundance = pd.Series(vals[np.arange(len(vals)), order[:, 0]], index=weights.index, name="dominant_abundance")
    return dominant, abundance


def dominant_from_doublet(weights: pd.DataFrame, results_df: pd.DataFrame) -> tuple[pd.Series, pd.Series, pd.Series]:
    """Derive dominant labels from doublet-mode RCTD outputs."""
    aligned = results_df.reindex(weights.index)
    spot_class = aligned["spot_class"].fillna("reject").astype(str)
    dominant = aligned["first_type"].fillna("reject").astype(str)
    abundance = pd.Series(
        [
            float(weights.loc[idx, dominant[idx]])
            if dominant[idx] in weights.columns
            else np.nan
            for idx in weights.index
        ],
        index=weights.index,
        name="dominant_abundance",
    )
    return dominant.rename("dominant_celltype"), abundance, spot_class.rename("rctd_spot_class")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rctd_dir", type=str, required=True, help="RCTD output directory")
    parser.add_argument("--meta_h5ad", type=str, required=True, help="Preprocessed spatial h5ad for coords/metadata")
    parser.add_argument("--out_h5ad", type=str, required=True, help="Output sp.h5ad path")
    parser.add_argument("--label", type=str, required=True)
    parser.add_argument("--ref_label", type=str, required=True)
    parser.add_argument("--ref_tag", type=str, required=True)
    parser.add_argument("--doublet_mode", type=str, default="doublet")
    parser.add_argument("--min_counts_per_bin", type=int, default=0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rctd_dir = Path(args.rctd_dir)
    weights = pd.read_csv(rctd_dir / "rctd_weights.csv", index_col=0)
    weights = normalize_rows(weights)

    results_df = None
    results_path = rctd_dir / "rctd_results_df.csv"
    if results_path.exists():
        results_df = pd.read_csv(results_path, index_col=0)

    if args.doublet_mode == "full":
        dominant, abundance = dominant_from_weights(weights)
        spot_class = pd.Series("full", index=weights.index, name="rctd_spot_class")
    else:
        if results_df is None:
            raise FileNotFoundError(f"Missing doublet results: {results_path}")
        dominant, abundance, spot_class = dominant_from_doublet(weights, results_df)

    meta = sc.read_h5ad(args.meta_h5ad)
    common = weights.index.intersection(meta.obs_names)
    if len(common) == 0:
        raise ValueError("No barcode overlap between RCTD weights and meta h5ad")
    weights = weights.loc[common]
    dominant = dominant.loc[common]
    abundance = abundance.loc[common]
    spot_class = spot_class.loc[common]

    adata = meta[common].copy()
    adata.obs["dominant_celltype"] = dominant.astype(str)
    adata.obs["dominant_abundance"] = abundance.astype(float)
    adata.obs["rctd_spot_class"] = spot_class.astype(str)
    adata.obsm["rctd_weights"] = weights.to_numpy(dtype=float)

    run_meta_path = rctd_dir / "rctd_run_meta.json"
    run_meta = {}
    if run_meta_path.exists():
        run_meta = json.loads(run_meta_path.read_text())
    adata.uns["rctd_run"] = {
        "label": args.label,
        "ref_label": args.ref_label,
        "ref_tag": args.ref_tag,
        "doublet_mode": args.doublet_mode,
        "min_counts_per_bin": args.min_counts_per_bin,
        "rctd_dir": str(rctd_dir),
        "cell_type_names": weights.columns.tolist(),
        **run_meta,
    }

    out_path = Path(args.out_h5ad)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(out_path)
    print(f"Wrote {out_path} with {adata.n_obs} spots and {len(weights.columns)} cell types")


if __name__ == "__main__":
    main()
