#!/usr/bin/env python
"""Compute Moran's I spatial autocorrelation on a preprocessed spatial h5ad.

This is intentionally separate from spatial_pipeline.py so the critical path
(preprocess → cell2location / Seurat) is not blocked by a slow SVG step.
"""

from __future__ import annotations

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, help="Input spatial h5ad")
    parser.add_argument("--label", required=True, help="Sample label (for output names)")
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory (e.g. data/processed/Tuyu_Apr2026)",
    )
    parser.add_argument("--n_perms", type=int, default=100, help="Moran's I permutations")
    parser.add_argument("--n_jobs", type=int, default=4, help="Parallel jobs for Moran's I")
    parser.add_argument("--n_top", type=int, default=6, help="Top SVGs to plot")
    parser.add_argument(
        "--n_genes",
        type=int,
        default=0,
        help="If >0, restrict to this many highly variable genes before Moran's I "
        "(much faster on Visium HD). 0 = all genes.",
    )
    parser.add_argument(
        "--ref_suffix",
        default="",
        help="Optional ref_tag filename suffix (e.g. _merged_nontumor_v1_subclass)",
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    adata = sc.read_h5ad(args.h5ad)
    print(f"Loaded {adata.n_obs} spots × {adata.n_vars} genes from {args.h5ad}")

    if "spatial_connectivities" not in adata.obsp:
        print("Building spatial neighbors (grid, n_neighs=6)")
        sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)

    if args.n_genes > 0 and adata.n_vars > args.n_genes:
        print(f"Selecting top {args.n_genes} highly variable genes for Moran's I")
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_genes, flavor="seurat")
        adata = adata[:, adata.var["highly_variable"]].copy()
        if "spatial_connectivities" not in adata.obsp:
            sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)

    print(f"Running Moran's I (n_perms={args.n_perms}, n_jobs={args.n_jobs})")
    sq.gr.spatial_autocorr(
        adata,
        mode="moran",
        n_perms=args.n_perms,
        n_jobs=args.n_jobs,
    )

    moran = adata.uns["moranI"].sort_values("I", ascending=False)
    csv_path = os.path.join(
        args.out_dir, f"spatial_{args.label}{args.ref_suffix}_moranI.csv"
    )
    moran.to_csv(csv_path)
    print(f"Saved: {csv_path}")
    print(moran.head(20))

    top_svgs = moran.head(args.n_top).index.tolist()
    if top_svgs:
        sq.pl.spatial_scatter(adata, color=top_svgs, size=1.5, ncols=3)
        png_path = os.path.join(
            args.out_dir, f"spatial_svgs_{args.label}{args.ref_suffix}.png"
        )
        plt.savefig(png_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Saved: {png_path}")

    print("Moran's I complete")


if __name__ == "__main__":
    main()
