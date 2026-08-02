#!/usr/bin/env python3
"""Build a slim harmonized AnnData with embeddings + Loupe annotations.

Reads the large combined normalized h5ad in backed mode (no full X matrix load),
merges per-sample for_loupe CSV annotations on (sample_id, barcode), and writes
a lightweight object suitable for labeling / clustering on X_pca_harmony.

Combined obs_names use anndata.concat(..., index_unique='-'), e.g.:
  s_008um_00469_00254-1-846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1
Loupe CSV Barcode column uses the plain Visium barcode:
  s_008um_00469_00254-1

Join key: obs['sample_id'] + barcode_base  <->  loupe sample_id + Barcode
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

DEFAULT_LABELS = [
    "Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um",
    "Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um",
    "Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J74_3_A1_GFP_R1_008um",
    "Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J75_3_D1_GFP_R1_008um",
]

DEFAULT_ANNOT_COLS = [
    "cell2location",
    "predicted.id",
    "cellestimate",
    "tlcelltype_assigned",
    "celltype_assigned_norm",
    "nnls_celltype_assigned",
    "norm_tlscore_gfp_probe",
    "norm_tltumscore_C11orf95_relafus1",
]

DEFAULT_OBSM_KEYS = ["X_pca_harmony", "X_umap", "spatial", "X_pca"]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def barcode_base_from_obs(names: pd.Index, sample_ids: pd.Series) -> pd.Series:
    """Strip '-{sample_id}' suffix added by anndata.concat(index_unique='-')."""
    bases = []
    for name, sid in zip(names.astype(str), sample_ids.astype(str)):
        suffix = f"-{sid}"
        if not name.endswith(suffix):
            raise ValueError(f"obs_name {name!r} does not end with sample suffix {suffix!r}")
        bases.append(name[: -len(suffix)])
    return pd.Series(bases, index=names, name="barcode_base")


def load_samples_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def label_to_sample_id(label: str, samples: pd.DataFrame) -> str:
    row = samples.loc[samples["label"] == label]
    if row.empty:
        raise KeyError(f"Label {label!r} not found in {samples}")
    return str(row.iloc[0]["sr_sample_dir"])


def load_loupe_annotations(
    loupe_dir: Path,
    labels: list[str],
    samples_tsv: Path,
    annot_cols: list[str],
) -> pd.DataFrame:
    samples = load_samples_tsv(samples_tsv)
    frames: list[pd.DataFrame] = []

    for label in labels:
        csv_path = loupe_dir / f"spatial_{label}_for_loupe.csv"
        if not csv_path.is_file():
            raise FileNotFoundError(csv_path)

        df = pd.read_csv(csv_path)
        if "Barcode" not in df.columns:
            raise ValueError(f"{csv_path} missing 'Barcode' column")

        sample_id = label_to_sample_id(label, samples)
        keep_cols = ["Barcode"] + [c for c in annot_cols if c in df.columns]
        missing = [c for c in annot_cols if c not in df.columns]
        if missing:
            print(f"Warning: {csv_path.name} missing columns: {missing}")

        chunk = df[keep_cols].copy()
        chunk["sample_id"] = sample_id
        chunk = chunk.rename(columns={"Barcode": "barcode_base"})
        frames.append(chunk)

    loupe = pd.concat(frames, ignore_index=True)
    loupe = loupe.drop_duplicates(subset=["sample_id", "barcode_base"], keep="first")
    return loupe


def report_merge_stats(obs: pd.DataFrame, annot_cols: list[str]) -> None:
    print("\nMerge summary by sample_id:")
    for sid, grp in obs.groupby("sample_id", sort=False):
        n = len(grp)
        if "cell2location" in grp.columns:
            matched = grp["cell2location"].notna().sum()
        else:
            matched = grp[annot_cols[0]].notna().sum() if annot_cols else 0
        print(f"  {sid}: {matched:,}/{n:,} spots annotated ({100 * matched / n:.1f}%)")


def build_slim_harmonized(
    input_h5ad: Path,
    output_h5ad: Path,
    loupe_dir: Path,
    samples_tsv: Path,
    labels: list[str],
    obsm_keys: list[str],
    annot_cols: list[str],
) -> ad.AnnData:
    print(f"Reading (backed): {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad, backed="r")

    if "sample_id" not in adata.obs.columns:
        raise KeyError("Combined h5ad obs must contain 'sample_id'")

    obs = adata.obs.copy()
    obs["barcode_base"] = barcode_base_from_obs(adata.obs_names, obs["sample_id"])

    loupe = load_loupe_annotations(loupe_dir, labels, samples_tsv, annot_cols)
    merge_cols = ["sample_id", "barcode_base"] + [
        c for c in annot_cols if c in loupe.columns
    ]
    obs = obs.merge(loupe[merge_cols], on=["sample_id", "barcode_base"], how="left")

    present_obsm = [k for k in obsm_keys if k in adata.obsm]
    missing_obsm = [k for k in obsm_keys if k not in adata.obsm]
    if missing_obsm:
        print(f"Warning: obsm keys not found and skipped: {missing_obsm}")
    if not present_obsm:
        raise KeyError(f"No requested obsm keys found in {input_h5ad}")

    print(f"Loading obsm arrays: {present_obsm}")
    obsm = {k: np.asarray(adata.obsm[k]) for k in present_obsm}

    slim = ad.AnnData(obs=obs, obsm=obsm)
    slim.uns["source_h5ad"] = str(input_h5ad)
    slim.uns["loupe_labels"] = labels
    slim.uns["annot_cols"] = [c for c in annot_cols if c in obs.columns]

    report_merge_stats(obs, annot_cols)

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nWriting slim h5ad ({slim.n_obs:,} cells, obsm={list(slim.obsm.keys())}): {output_h5ad}")
    slim.write_h5ad(output_h5ad)
    return slim


def parse_args() -> argparse.Namespace:
    root = repo_root()
    default_input = (
        root
        / "data/processed/Tuyu_Apr2026/spatial_combined_normalized_D90_7_D1_D90_5_A1_J74_3_A1_J75_3_D1.h5ad"
    )
    default_output = (
        root
        / "data/processed/Tuyu_Apr2026/spatial_combined_slim_harmony_annotated_D90_7_D1_D90_5_A1_J74_3_A1_J75_3_D1.h5ad"
    )

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-h5ad", type=Path, default=default_input)
    p.add_argument("--output-h5ad", type=Path, default=default_output)
    p.add_argument("--loupe-dir", type=Path, default=root / "for_loupe")
    p.add_argument("--samples-tsv", type=Path, default=root / "config/samples.tsv")
    p.add_argument(
        "--labels",
        nargs="+",
        default=DEFAULT_LABELS,
        help="Snakemake-style labels used in for_loupe CSV filenames",
    )
    p.add_argument(
        "--obsm-keys",
        nargs="+",
        default=DEFAULT_OBSM_KEYS,
        help="obsm arrays to copy into slim object",
    )
    p.add_argument(
        "--annot-cols",
        nargs="+",
        default=DEFAULT_ANNOT_COLS,
        help="Loupe CSV columns to merge into obs",
    )
    p.add_argument(
        "--cell2location-only",
        action="store_true",
        help="Only merge cell2location column (shorthand for --annot-cols cell2location)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    annot_cols = ["cell2location"] if args.cell2location_only else args.annot_cols
    build_slim_harmonized(
        input_h5ad=args.input_h5ad,
        output_h5ad=args.output_h5ad,
        loupe_dir=args.loupe_dir,
        samples_tsv=args.samples_tsv,
        labels=args.labels,
        obsm_keys=args.obsm_keys,
        annot_cols=annot_cols,
    )


if __name__ == "__main__":
    main()
