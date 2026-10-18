"""Shared preprocessing for cell2location and RCTD deconvolution."""

from __future__ import annotations

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import scipy.sparse as sp


def ensure_raw_counts(adata, name: str = "adata"):
    """Ensure AnnData ``X`` holds integer raw counts."""
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    elif adata.raw is not None:
        adata = adata.raw.to_adata()
    if sp.issparse(adata.X):
        adata.X.data = np.round(adata.X.data).astype(np.int32)
        adata.X.eliminate_zeros()
    else:
        adata.X = np.round(adata.X).astype(np.int32)
    assert adata.X.max() > 20, f"{name}: X does not look like raw counts"
    if sp.issparse(adata.X):
        assert np.all(adata.X.data % 1 == 0), f"{name}: counts not integer"
    else:
        assert np.all(adata.X % 1 == 0), f"{name}: counts not integer"
    return adata


def clean_labels(adata, label_key: str, batch_key: str):
    """Drop invalid labels and cell types with fewer than five cells."""
    lab = adata.obs[label_key].astype(str).str.strip()
    batch = adata.obs[batch_key].astype(str).str.strip()
    keep = ~lab.isin(["nan", "NaN", "None", ""]) & lab.notna()
    adata = adata[keep].copy()
    lab = lab[keep]
    batch = batch[keep]
    vc = lab.value_counts()
    keep_types = vc[vc >= 5].index
    adata = adata[lab.isin(keep_types)].copy()
    adata.obs[label_key] = pd.Categorical(
        adata.obs.loc[adata.obs_names, label_key],
        categories=sorted(keep_types),
    )
    adata.obs[batch_key] = pd.Categorical(adata.obs[batch_key].astype(str))
    return adata


def harmonize_genes(adata, inf_aver):
    """Subset spatial AnnData and reference signatures to shared genes."""
    if adata.var_names.duplicated().any():
        print(f"Deduplicating {adata.var_names.duplicated().sum()} spatial gene names")
        adata.var_names_make_unique()
    inf_aver = inf_aver[~inf_aver.index.duplicated(keep="first")].copy()
    shared = pd.Index(adata.var_names).intersection(inf_aver.index)
    if len(shared) == 0:
        raise ValueError("No shared genes between spatial and reference")
    print(f"Shared genes: {len(shared)}")
    adata = adata[:, shared].copy()
    inf_aver = inf_aver.loc[shared, :].copy()
    return adata, inf_aver


def harmonize_genes_with_ref(adata, adata_ref):
    """Subset spatial and reference AnnData objects to shared genes."""
    if adata.var_names.duplicated().any():
        print(f"Deduplicating {adata.var_names.duplicated().sum()} spatial gene names")
        adata.var_names_make_unique()
    if adata_ref.var_names.duplicated().any():
        adata_ref.var_names_make_unique()
    shared = pd.Index(adata.var_names).intersection(adata_ref.var_names)
    if len(shared) == 0:
        raise ValueError("No shared genes between spatial and reference")
    print(f"Shared genes: {len(shared)}")
    return adata[:, shared].copy(), adata_ref[:, shared].copy()


def load_raw_spatial(sr_out, bin_size: str = "008um", meta_h5ad: str | None = None):
    """Load raw Visium HD binned counts and optionally intersect with meta h5ad."""
    bin_dir = f"{sr_out}/binned_outputs/square_{bin_size}"
    adata = sq.read.visium(bin_dir, counts_file="filtered_feature_bc_matrix.h5")
    adata = ensure_raw_counts(adata, "spatial")
    if meta_h5ad is not None:
        meta = sc.read_h5ad(meta_h5ad)
        common = adata.obs_names.intersection(meta.obs_names)
        print(f"Spatial overlap with metadata h5ad: {len(common)} / {adata.n_obs}")
        if len(common) == 0:
            raise ValueError(
                f"No barcode overlap between spatial bins ({bin_size}, n={adata.n_obs}) "
                f"and metadata h5ad ({meta_h5ad}, n={meta.n_obs}). "
                "Check that --bin_um matches the bin size used to build meta_h5ad."
            )
        adata = adata[common].copy()
        for col in meta.obs.columns:
            if col not in adata.obs.columns:
                adata.obs[col] = meta.obs.loc[common, col]
        if "spatial" in meta.obsm and "spatial" not in adata.obsm:
            adata.obsm["spatial"] = meta[common].obsm["spatial"].copy()
    adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
    adata = adata[:, ~adata.var["mt"]].copy()
    return adata


def filter_sparse_bins(adata, min_counts: int):
    """Filter spatial bins with fewer than ``min_counts`` UMIs."""
    if min_counts <= 0:
        return adata
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=min_counts)
    print("filtering spatial adata")
    print(f"Bins with >= {min_counts} UMIs: {adata.n_obs} / {n_before}")
    return adata


def filter_reference_genes(adata_ref):
    """Apply the same gene filter used by cell2location reference training."""
    from cell2location.utils.filtering import filter_genes

    selected = filter_genes(
        adata_ref,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12,
    )
    return adata_ref[:, selected].copy()
