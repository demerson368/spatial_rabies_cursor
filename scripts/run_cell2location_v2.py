"""Train cell2location on Visium HD 008um bins using a scRNA-seq reference.

Loads or trains reference signatures (RegressionModel), deconvolves spatial
counts (Cell2location), and writes ``sp.h5ad`` with posterior abundances and
``dominant_celltype`` / ``dominant_abundance`` in ``obs``. Run params are stored
in ``adata.uns['cell2location_run']``.
"""
from __future__ import annotations

import os
# Let LSF set CUDA_VISIBLE_DEVICES for the allocated GPU; do not override.
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")

import torch
#Enable Tensor Core optimization for A100
torch.set_float32_matmul_precision('high')

import matplotlib
matplotlib.use("Agg")   # non-interactive backend, MUST be set before importing pyplot
import matplotlib.pyplot as plt

import gc
import argparse
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.cm as cm
import scipy.sparse as sp
from scipy.sparse import issparse
import squidpy as sq
import math
import cell2location
from cell2location.models import RegressionModel
from cell2location.models import Cell2location

from rctd_preprocess import (
    clean_labels,
    ensure_raw_counts,
    filter_reference_genes,
    filter_sparse_bins,
    harmonize_genes,
    load_raw_spatial,
)




def free_gpu_memory():
    gc.collect()
    try:
        import torch

        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            torch.cuda.synchronize()
    except ImportError:
        pass


def ensure_dir(path):

    os.makedirs(path, exist_ok=True)

def extract_inf_aver(adata_ref):

    factor_names = adata_ref.uns["mod"]["factor_names"]
    if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
            [f"means_per_cluster_mu_fg_{i}" for i in factor_names]
        ].copy()
    else:
        inf_aver = adata_ref.var[
            [f"means_per_cluster_mu_fg_{i}" for i in factor_names]
        ].copy()
    inf_aver.columns = factor_names

    return inf_aver


def c2l_n_cells_dir_suffix(n_cells_per_location: int) -> str:
    """Directory suffix when ``N_cells_per_location`` differs from legacy default 2.

    Args:
        n_cells_per_location: Cell2location ``N_cells_per_location`` hyperparameter.

    Returns:
        str: Empty string for 2; otherwise ``_n{n}`` (e.g. ``_n1``).
    """
    if n_cells_per_location == 2:
        return ""
    return f"_n{n_cells_per_location}"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for cell2location spatial deconvolution.

    Returns:
        argparse.Namespace: Parsed paths, sample label, and model hyperparameters.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref_out", type=str, help="ref input directory")
    parser.add_argument("--ref_label", type=str, help="ref name label")
    parser.add_argument(
        "--ref_file",
        type=str,
        help="ref input .h5ad. e.g. integrated_ref_rawcounts_with_labels.h5ad",
    )
    parser.add_argument(
        "--ref_tag",
        type=str,
        default="v2",
        help="tag for ref/spatial output dirs (matches config ref_tag)",
    )
    parser.add_argument("--sr_out", type=str, help="spatial input spot dir")
    parser.add_argument(
        "--label",
        type=str,
        help="spatial name label. Assumed in Tuyu_Apr2026 subdirectory",
    )
    parser.add_argument(
        "--cohort",
        type=str,
        default=None,
        help="Cohort subdirectory under results_cell2location/ and data/processed/",
    )
    parser.add_argument(
        "--min_counts_per_bin",
        type=int,
        default=0,
        help="minimum high confidence UMI per spot (aim for 80k - 120k spots)",
    )
    parser.add_argument(
        "--bin_um",
        type=str,
        default=None,
        help="Visium HD bin size folder (e.g. 008um, 016um). Defaults to 008um.",
    )
    parser.add_argument(
        "--meta_h5ad",
        type=str,
        default=None,
        help="Preprocess or filtered h5ad for spot barcode intersection",
    )
    parser.add_argument(
        "--spatial_dir",
        type=str,
        default=None,
        help="Output directory for spatial cell2location model (parent of sp.h5ad)",
    )
    parser.add_argument(
        "--n_cells_per_location",
        type=int,
        default=2,
        help="Cell2location N_cells_per_location (use 1 for single-type 8um bins)",
    )
    parser.add_argument(
        "--detection_alpha",
        type=float,
        default=200.0,
        help="Cell2location detection_alpha hyperparameter",
    )
    return parser.parse_args()


def main() -> None:
    """Load reference signatures, train Cell2location on spatial counts, write sp.h5ad."""
    args = parse_args()

    label = args.label
    ref_label = args.ref_label
    n_cells_per_location = args.n_cells_per_location
    detection_alpha = args.detection_alpha

    print("ref_label: ", ref_label)
    print("label: ", label)
    print("n_cells_per_location: ", n_cells_per_location)
    print("detection_alpha: ", detection_alpha)

    if args.meta_h5ad:
        meta_h5ad = args.meta_h5ad
    else:
        meta_h5ad = (
            "/research/groups/mackgrp/home/demerson/spatial_rabies_cursor/"
            f"data/processed/Tuyu_Apr2026/spatial_{label}.h5ad"
        )

    cohort = args.cohort
    if not cohort and "/processed/" in meta_h5ad:
        cohort = meta_h5ad.split("/processed/")[1].split("/")[0]
    if not cohort:
        raise ValueError(
            "Could not infer cohort; pass --cohort or use --meta_h5ad under data/processed/{cohort}/"
        )

    ref_out = args.ref_out
    sr_out = args.sr_out

    print("cohort: ", cohort)
    print("ref_out: ", ref_out)
    print("sr_out: ", sr_out)

    ref_file = args.ref_out + args.ref_file
    ref_tag = args.ref_tag

    print("ref_file: ", ref_file)
    print("ref_tag: ", ref_tag)

    ref_dir = f"../results_cell2location/{cohort}/ref_{ref_label}_{ref_tag}"
    n_cells_suffix = c2l_n_cells_dir_suffix(n_cells_per_location)
    if args.spatial_dir:
        spatial_dir = args.spatial_dir
    else:
        spatial_dir = (
            f"../results_cell2location/{cohort}/c2l_spatial_{label}_{ref_tag}{n_cells_suffix}"
        )

    print("ref_dir: ", ref_dir)
    print("spatial_dir: ", spatial_dir)

    min_counts_per_bin = args.min_counts_per_bin
    bin_um = args.bin_um or "008um"

    ensure_dir(ref_dir)
    ensure_dir(spatial_dir)

    #Keep collective only for the single spatial slide.  Batch = dJ33,J51,D87,etc.
    label_key = "cell_type_post"
    ref_batch_key = "batch"          
    spatial_batch_key = "collective"
    train_kwargs_gpu = {"accelerator": "gpu"}
    sample_kwargs_gpu = {"accelerator": "gpu"}

    print(
        f"CUDA available={torch.cuda.is_available()}, "
        f"device_count={torch.cuda.device_count()}, "
        f"CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES', 'unset')}"
    )

    # ----- REFERENCE -----
    
    #Notes: original model performed per nontumor batch (J58, J59, etc.).  Model at ref_dir.  

    adata_ref = sc.read_h5ad(ref_file)
    if adata_ref.var_names.duplicated().any():
        adata_ref.var_names_make_unique()

    adata_ref = ensure_raw_counts(adata_ref, "reference")
    adata_ref = clean_labels(adata_ref, label_key, ref_batch_key)
    adata_ref = filter_reference_genes(adata_ref)

    ref_h5ad = f"{ref_dir}/sc.h5ad"
    if os.path.exists(ref_h5ad):
        print("Loading existing reference signatures:", ref_h5ad)
        adata_ref = sc.read_h5ad(ref_h5ad)
        inf_aver = extract_inf_aver(adata_ref)
    else:
        RegressionModel.setup_anndata(
            adata=adata_ref,
            batch_key=ref_batch_key,
            labels_key=label_key,
            categorical_covariate_keys=None,
        )
        mod_ref = RegressionModel(adata_ref)
        mod_ref.train(max_epochs=250, batch_size=2500,datasplitter_kwargs={"num_workers": 4}, **train_kwargs_gpu)

        fig, ax = plt.subplots(figsize=(8, 5))
        mod_ref.plot_history(20,ax=ax)
        fig.savefig(f"{ref_dir}/ref_training_history.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        adata_ref = mod_ref.export_posterior(
            adata_ref,
            sample_kwargs={"num_samples": 1000, "batch_size": None, **sample_kwargs_gpu},  ##reduce from 2500 to 1000 for cuda memory?
        )
        mod_ref.save(ref_dir, overwrite=True)
        adata_ref.write(ref_h5ad)
        inf_aver = extract_inf_aver(adata_ref)
        inf_aver.to_csv(f"{ref_dir}/reference_signatures.csv")

    # ----- SPATIAL -----
    adata = load_raw_spatial(sr_out, bin_size=bin_um, meta_h5ad=meta_h5ad)
    adata.obs[spatial_batch_key] = pd.Categorical(["nontumor"] * adata.n_obs)
    adata, inf_aver = harmonize_genes(adata, inf_aver)
    adata = filter_sparse_bins(adata, min_counts_per_bin)
    if adata.n_obs == 0:
        raise ValueError(
            f"Spatial AnnData is empty after loading/filtering (bin_um={bin_um}). "
            "Cannot run Cell2location."
        )

    Cell2location.setup_anndata(adata=adata, batch_key=spatial_batch_key)
    mod = Cell2location(
        adata,
        cell_state_df=inf_aver,
        N_cells_per_location=n_cells_per_location,
        detection_alpha=detection_alpha,
    )
    mod.train(max_epochs=30000, batch_size=None, train_size=1,datasplitter_kwargs={"num_workers": 4}, **train_kwargs_gpu)

    fig, ax = plt.subplots(figsize=(8, 5))
    mod.plot_history(1000,ax=ax)
    #plt.gcf().savefig(f"{spatial_dir}/spatial_training_history.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{spatial_dir}/spatial_training_history.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    #torch.cuda.empty_cache()
    free_gpu_memory()
    adata = mod.export_posterior(
        adata,
        sample_kwargs={
            "num_samples": 1000,
            "batch_size": min(1000, adata.n_obs), #reduce from 2000 to 1000 for cuda memory?
            **sample_kwargs_gpu,
        },
    )
    mod.save(spatial_dir, overwrite=True)
    adata.write(f"{spatial_dir}/sp.h5ad")

    # ----- EXPORT ABUNDANCES -----
    factor_names = adata.uns["mod"]["factor_names"]
    abund_df = pd.DataFrame(
        np.asarray(adata.obsm["means_cell_abundance_w_sf"]),
        index=adata.obs_names,
        columns=factor_names,
    )
    print("NaN fraction:", abund_df.isna().mean().mean())
    adata.obs["dominant_celltype"] = abund_df.idxmax(axis=1)
    adata.obs["dominant_abundance"] = abund_df.max(axis=1)
    adata.uns["cell2location_run"] = {
        "cohort": cohort,
        "label": label,
        "ref_label": ref_label,
        "ref_tag": ref_tag,
        "n_cells_per_location": n_cells_per_location,
        "detection_alpha": detection_alpha,
        "min_counts_per_bin": min_counts_per_bin,
        "spatial_dir": spatial_dir,
    }
    adata.write(f"{spatial_dir}/sp.h5ad")

if __name__ == '__main__':
    main()



