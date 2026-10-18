"""Prepare RCTD reference and spatial inputs as Matrix Market + metadata for R."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix, issparse

from rctd_preprocess import (
    clean_labels,
    ensure_raw_counts,
    filter_reference_genes,
    filter_sparse_bins,
    harmonize_genes_with_ref,
    load_raw_spatial,
)

LABEL_KEY = "cell_type_post"
REF_BATCH_KEY = "batch"


def ensure_dir(path: str | Path) -> None:
    os.makedirs(path, exist_ok=True)


def write_mtx_bundle(
    out_dir: Path,
    counts_cells_by_genes,
    barcodes: pd.Index,
    genes: pd.Index,
    prefix: str,
) -> None:
    """Write genes x cells sparse matrix and sidecar files for R."""
    ensure_dir(out_dir)
    if issparse(counts_cells_by_genes):
        mat = counts_cells_by_genes.T.tocsr()
    else:
        mat = csr_matrix(np.asarray(counts_cells_by_genes).T)
    mmwrite(out_dir / f"{prefix}_counts.mtx", mat)
    pd.Series(genes, name="gene").to_csv(out_dir / f"{prefix}_genes.txt", index=False, header=False)
    pd.Series(barcodes, name="barcode").to_csv(out_dir / f"{prefix}_barcodes.txt", index=False, header=False)


def build_reference_rds(bundle_dir: Path, reference_rds: Path, rscript: str) -> None:
    """Call R helper to assemble Reference inputs into an RDS cache."""
    script = Path(__file__).resolve().parent / "build_rctd_reference.R"
    subprocess.run(
        [
            rscript,
            str(script),
            "--bundle_dir",
            str(bundle_dir),
            "--out_rds",
            str(reference_rds),
        ],
        check=True,
    )


def prepare_reference(
    ref_file: str,
    ref_label: str,
    ref_tag: str,
    ref_dir: Path,
    rscript: str,
) -> Path:
    """Load, filter, and cache reference inputs for RCTD."""
    reference_rds = ref_dir / "reference.rds"
    bundle_dir = ref_dir / "reference_bundle"
    if reference_rds.exists():
        print(f"Using cached reference RDS: {reference_rds}")
        return reference_rds

    adata_ref = sc.read_h5ad(ref_file)
    if adata_ref.var_names.duplicated().any():
        adata_ref.var_names_make_unique()
    adata_ref = ensure_raw_counts(adata_ref, "reference")
    adata_ref = clean_labels(adata_ref, LABEL_KEY, REF_BATCH_KEY)
    adata_ref = filter_reference_genes(adata_ref)

    write_mtx_bundle(
        bundle_dir,
        adata_ref.X,
        adata_ref.obs_names,
        adata_ref.var_names,
        "reference",
    )
    cell_types = pd.DataFrame(
        {"barcode": adata_ref.obs_names, "cell_type": adata_ref.obs[LABEL_KEY].astype(str)},
    )
    cell_types.to_csv(bundle_dir / "reference_cell_types.csv", index=False)
    meta = {
        "label_key": LABEL_KEY,
        "ref_label": ref_label,
        "ref_tag": ref_tag,
        "n_cells": int(adata_ref.n_obs),
        "n_genes": int(adata_ref.n_vars),
        "cell_types": sorted(adata_ref.obs[LABEL_KEY].astype(str).unique().tolist()),
    }
    (bundle_dir / "reference_meta.json").write_text(json.dumps(meta, indent=2))
    build_reference_rds(bundle_dir, reference_rds, rscript)
    return reference_rds


def prepare_spatial(
    sr_out: str,
    label: str,
    bin_um: str,
    meta_h5ad: str,
    min_counts_per_bin: int,
    spatial_dir: Path,
    reference_rds: Path,
    rscript: str,
) -> Path:
    """Load spatial counts, harmonize with reference genes, export spatial bundle."""
    spatial_rds = spatial_dir / "spatial_input.rds"
    bundle_dir = spatial_dir / "spatial_bundle"
    if spatial_rds.exists():
        print(f"Using cached spatial input RDS: {spatial_rds}")
        return spatial_rds

    adata = load_raw_spatial(sr_out, bin_size=bin_um, meta_h5ad=meta_h5ad)
    adata = filter_sparse_bins(adata, min_counts_per_bin)
    if adata.n_obs == 0:
        raise ValueError(
            f"Spatial AnnData is empty after loading/filtering (bin_um={bin_um}). "
            "Cannot run RCTD."
        )

    ref_bundle = reference_rds.parent / "reference_bundle"
    genes = pd.read_csv(ref_bundle / "reference_genes.txt", header=None).iloc[:, 0].astype(str)
    barcodes = pd.read_csv(ref_bundle / "reference_barcodes.txt", header=None).iloc[:, 0].astype(str)
    ref_mat = mmread(ref_bundle / "reference_counts.mtx").tocsr().T  # cells x genes
    adata_ref = sc.AnnData(ref_mat)
    adata_ref.obs_names = barcodes
    adata_ref.var_names = genes
    adata, adata_ref = harmonize_genes_with_ref(adata, adata_ref)

    write_mtx_bundle(
        bundle_dir,
        adata.X,
        adata.obs_names,
        adata.var_names,
        "spatial",
    )
    if "spatial" in adata.obsm:
        coords = pd.DataFrame(
            adata.obsm["spatial"],
            index=adata.obs_names,
            columns=["x", "y"],
        )
    else:
        raise ValueError("Spatial coordinates missing from adata.obsm['spatial']")
    coords.to_csv(bundle_dir / "spatial_coords.csv")
    meta = {
        "label": label,
        "bin_um": bin_um,
        "min_counts_per_bin": min_counts_per_bin,
        "n_pixels": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "meta_h5ad": meta_h5ad,
    }
    (bundle_dir / "spatial_meta.json").write_text(json.dumps(meta, indent=2))

    script = Path(__file__).resolve().parent / "build_rctd_spatial_input.R"
    subprocess.run(
        [
            rscript,
            str(script),
            "--bundle_dir",
            str(bundle_dir),
            "--out_rds",
            str(spatial_rds),
        ],
        check=True,
    )
    return spatial_rds


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["reference", "spatial"], required=True)
    parser.add_argument("--ref_out", type=str, help="Reference input directory")
    parser.add_argument("--ref_label", type=str, help="Reference label")
    parser.add_argument("--ref_file", type=str, help="Reference h5ad filename")
    parser.add_argument("--ref_tag", type=str, default="v2")
    parser.add_argument("--ref_dir", type=str, help="Output directory for reference RDS")
    parser.add_argument("--sr_out", type=str, help="Space Ranger outs directory")
    parser.add_argument("--label", type=str, help="Spatial sample label")
    parser.add_argument("--bin_um", type=str, default="008um")
    parser.add_argument("--meta_h5ad", type=str, default=None)
    parser.add_argument("--min_counts_per_bin", type=int, default=0)
    parser.add_argument("--spatial_dir", type=str, help="Output directory for spatial RCTD run")
    parser.add_argument("--reference_rds", type=str, default=None)
    parser.add_argument(
        "--rscript",
        type=str,
        default="Rscript",
        help="Path to Rscript binary",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode == "reference":
        if not all([args.ref_out, args.ref_label, args.ref_file, args.ref_dir]):
            raise ValueError("reference mode requires --ref_out, --ref_label, --ref_file, --ref_dir")
        ref_file = str(Path(args.ref_out) / args.ref_file)
        prepare_reference(
            ref_file=ref_file,
            ref_label=args.ref_label,
            ref_tag=args.ref_tag,
            ref_dir=Path(args.ref_dir),
            rscript=args.rscript,
        )
        return

    if not all([args.sr_out, args.label, args.spatial_dir, args.reference_rds]):
        raise ValueError("spatial mode requires --sr_out, --label, --spatial_dir, --reference_rds")
    #meta_h5ad = args.meta_h5ad
    #if not meta_h5ad:
    #    raise ValueError("--meta_h5ad is required for spatial mode")
    if args.mode == "spatial":
        meta_h5ad = args.meta_h5ad
        if not meta_h5ad:
            raise ValueError("--meta_h5ad is required for spatial mode")
        prepare_spatial(
            sr_out=args.sr_out,
            label=args.label,
            bin_um=args.bin_um or "008um",
            meta_h5ad=meta_h5ad,
            min_counts_per_bin=args.min_counts_per_bin,
            spatial_dir=Path(args.spatial_dir),
            reference_rds=Path(args.reference_rds),
            rscript=args.rscript,
        )


if __name__ == "__main__":
    main()
