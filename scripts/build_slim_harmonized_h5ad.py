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
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
DEFAULT_OUTPUT_LABEL = "D90_5_A1_D90_7_D1_J74_3_A1_J75_3_D1"
DEFAULT_UMAP_COLOR_COLS = ("cell2location", "cellestimate", "sample_id")


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
    """Load sample metadata table.

    Args:
        path: TSV path with at least ``label`` and ``sr_sample_dir`` columns.

    Returns:
        pd.DataFrame: Sample metadata rows.
    """
    return pd.read_csv(path, sep="\t")


def resolve_samples_row(loupe_label: str, samples: pd.DataFrame) -> pd.Series:
    """Map a for_loupe filename stem to a ``samples.tsv`` row.

    Resolution order:
    1. Exact match on ``label``
    2. Exact match on optional ``loupe_label`` column
    3. Unique prefix match where ``loupe_label`` equals or extends ``label``
       (e.g. ``..._008um_merged_nontumor_v1_broad`` -> ``..._008um``)

    Args:
        loupe_label: Stem between ``spatial_`` and ``_for_loupe.csv``.
        samples: Sample metadata loaded from ``config/samples.tsv``.

    Returns:
        pd.Series: Matching metadata row.

    Raises:
        KeyError: If no row matches or multiple rows tie for the longest prefix.
    """
    loupe_label = loupe_label.strip()
    exact = samples.loc[samples["label"] == loupe_label]
    if len(exact) == 1:
        return exact.iloc[0]
    if len(exact) > 1:
        raise KeyError(f"Multiple samples.tsv rows for label {loupe_label!r}")

    if "loupe_label" in samples.columns:
        by_loupe = samples.loc[samples["loupe_label"].astype(str) == loupe_label]
        if len(by_loupe) == 1:
            return by_loupe.iloc[0]
        if len(by_loupe) > 1:
            raise KeyError(f"Multiple samples.tsv rows for loupe_label {loupe_label!r}")

    candidates: list[tuple[int, pd.Series]] = []
    for _, row in samples.iterrows():
        base = str(row["label"])
        if loupe_label == base or loupe_label.startswith(f"{base}_"):
            candidates.append((len(base), row))

    if not candidates:
        raise KeyError(
            f"Label {loupe_label!r} not found in samples.tsv "
            f"(no exact or prefix match on label/loupe_label)"
        )

    candidates.sort(key=lambda item: item[0], reverse=True)
    best_len = candidates[0][0]
    best_rows = [row for length, row in candidates if length == best_len]
    if len(best_rows) > 1:
        bases = [str(row["label"]) for row in best_rows]
        raise KeyError(
            f"Ambiguous loupe label {loupe_label!r}: equally matches samples.tsv labels {bases}"
        )
    return best_rows[0]


def label_to_sample_id(loupe_label: str, samples: pd.DataFrame) -> str:
    """Return ``sr_sample_dir`` for a Loupe or spatial sample label.

    Args:
        loupe_label: Loupe CSV stem or canonical ``samples.tsv`` label.
        samples: Sample metadata loaded from ``config/samples.tsv``.

    Returns:
        str: Space Ranger sample directory id used in combined h5ad ``sample_id``.
    """
    row = resolve_samples_row(loupe_label, samples)
    return str(row["sr_sample_dir"])


def resolve_loupe_csv_path(loupe_dir: Path, label: str, ref_tag: str | None) -> Path:
    """Resolve a Loupe export CSV for a base or ref-tagged sample label.

    Args:
        loupe_dir: Directory containing ``spatial_*_for_loupe.csv`` files.
        label: Base ``samples.tsv`` label or full Loupe filename stem.
        ref_tag: Optional ref suffix (e.g. ``merged_nontumor_v1_broad``).

    Returns:
        Path: Existing Loupe CSV path.

    Raises:
        FileNotFoundError: If no matching Loupe CSV exists.
    """
    label = label.strip()
    candidates: list[Path] = []

    if ref_tag:
        candidates.append(loupe_dir / f"spatial_{label}_{ref_tag}_for_loupe.csv")
    candidates.append(loupe_dir / f"spatial_{label}_for_loupe.csv")

    seen: set[Path] = set()
    for path in candidates:
        if path in seen:
            continue
        seen.add(path)
        if path.is_file():
            return path

    raise FileNotFoundError(
        f"No Loupe CSV found for label {label!r} in {loupe_dir} "
        f"(tried ref_tag={ref_tag!r})"
    )


def load_loupe_annotations(
    loupe_dir: Path,
    labels: list[str],
    samples_tsv: Path,
    annot_cols: list[str],
    ref_tag: str | None = None,
) -> pd.DataFrame:
    """Load and concatenate Loupe CSV annotations for requested samples.

    Args:
        loupe_dir: Directory with ``spatial_*_for_loupe.csv`` exports.
        labels: Base ``samples.tsv`` labels or full Loupe filename stems.
        samples_tsv: Sample metadata table path.
        annot_cols: Loupe columns to merge into obs.
        ref_tag: Optional ref suffix when ``labels`` are base sample labels.

    Returns:
        pd.DataFrame: Concatenated annotations keyed by ``sample_id`` and
        ``barcode_base``.
    """
    samples = load_samples_tsv(samples_tsv)
    frames: list[pd.DataFrame] = []

    for label in labels:
        csv_path = resolve_loupe_csv_path(loupe_dir, label, ref_tag)
        loupe_label = csv_path.name[len("spatial_") : -len("_for_loupe.csv")]

        df = pd.read_csv(csv_path)
        if "Barcode" not in df.columns:
            raise ValueError(f"{csv_path} missing 'Barcode' column")

        sample_id = label_to_sample_id(loupe_label, samples)
        keep_cols = ["Barcode"] + [c for c in annot_cols if c in df.columns]
        missing = [c for c in annot_cols if c not in df.columns]
        if missing:
            print(f"Warning: {csv_path.name} missing columns: {missing}")

        chunk = df[keep_cols].copy()
        chunk["sample_id"] = sample_id
        chunk = chunk.rename(columns={"Barcode": "barcode_base"})
        frames.append(chunk)
        print(f"Loupe merge: {csv_path.name} -> sample_id={sample_id!r}")

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
    ref_tag: str | None = None,
) -> ad.AnnData:
    print(f"Reading (backed): {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad, backed="r")

    if "sample_id" not in adata.obs.columns:
        raise KeyError("Combined h5ad obs must contain 'sample_id'")

    obs = adata.obs.copy()
    obs["barcode_base"] = barcode_base_from_obs(adata.obs_names, obs["sample_id"])

    loupe = load_loupe_annotations(
        loupe_dir, labels, samples_tsv, annot_cols, ref_tag=ref_tag
    )
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


def relabel_umap_legend_with_counts(
    ax: plt.Axes,
    adata: ad.AnnData,
    col: str,
) -> None:
    """Append spot counts to UMAP legend labels without changing category order.

    Args:
        ax: Matplotlib axes returned by :func:`scanpy.pl.umap`.
        adata: AnnData used for the plot (counts reflect current obs rows).
        col: ``obs`` column whose legend entries receive ``(n=...)`` suffixes.

    Returns:
        None
    """
    legend = ax.get_legend()
    if legend is None:
        print(f"Warning: no legend found for UMAP colored by {col!r}")
        return

    counts = adata.obs[col].value_counts()
    for text in legend.get_texts():
        label = text.get_text()
        if label in counts.index:
            text.set_text(f"{label} (n={counts[label]:,})")


def plot_umap_with_count_legend(
    adata: ad.AnnData,
    col: str,
    out_path: Path,
    *,
    title: str | None = None,
    figsize: tuple[float, float] = (8.0, 6.0),
    dpi: int = 200,
) -> None:
    """Save a UMAP panel with per-category spot counts in the legend.

    Args:
        adata: AnnData with ``X_umap`` and the requested ``obs`` column.
        col: Metadata column to color UMAP points.
        out_path: Output PNG path.
        title: Figure title; defaults to ``col``.
        figsize: Matplotlib figure size (width, height) in inches.
        dpi: PNG resolution.

    Returns:
        None

    Raises:
        KeyError: If ``col`` is absent from ``adata.obs``.
    """
    if col not in adata.obs.columns:
        raise KeyError(f"Column {col!r} not in adata.obs: {list(adata.obs.columns)}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=figsize)
    sc.pl.umap(
        adata,
        color=col,
        ax=ax,
        show=False,
        title=title or col,
        frameon=False,
        legend_loc="right margin",
    )
    relabel_umap_legend_with_counts(ax, adata, col)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def parse_args() -> argparse.Namespace:
    root = repo_root()
    default_input = (
        root
        / "data/processed/Tuyu_Apr2026/spatial_combined_normalized_D90_7_D1_D90_5_A1_J74_3_A1_J75_3_D1.h5ad"
    )
    default_output_dir = (
        root
        / "data/processed/Tuyu_Apr2026/"
    )

    default_output_h5ad = (
        "spatial_combined_slim_harmony_annotated_D90_7_D1_D90_5_A1_J74_3_A1_J75_3_D1.h5ad"
    )

    default_output_label = (
        "combined_slim_harmony_annotated_D90_7_D1_D90_5_A1_J74_3_A1_J75_3_D1"
    )

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-h5ad", type=Path, default=default_input)
    p.add_argument("--output-dir", type=Path, default=default_output_dir)
    p.add_argument("--output-h5ad", type=Path, default=default_output_h5ad)
    p.add_argument("--loupe-dir", type=Path, default=root / "for_loupe")
    p.add_argument("--samples-tsv", type=Path, default=root / "config/samples.tsv")
    p.add_argument(
        "--labels",
        nargs="+",
        default=DEFAULT_LABELS,
        help="samples.tsv labels or full for_loupe filename stems "
        "(e.g. with _merged_nontumor_v1_broad suffix)",
    )
    p.add_argument(
        "--ref-tag",
        default=None,
        help="Loupe ref suffix when --labels are base samples.tsv labels "
        "(e.g. merged_nontumor_v1_broad; matches config ref_tag)",
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
    p.add_argument(
        "--output_label",
        default=default_output_label,
        help="label for output umap .png",
    )

    return p.parse_args()


def main() -> None:
    args = parse_args()
    #annot_cols = ["cell2location"] if args.cell2location_only else args.annot_cols
    annot_cols = ["cell2location","cellestimate"]

    args.output_h5ad = args.output_dir / args.output_h5ad
    print("output h5ad: ",args.output_h5ad)
    adata = build_slim_harmonized(
        input_h5ad=args.input_h5ad,
        output_h5ad=args.output_h5ad,
        loupe_dir=args.loupe_dir,
        samples_tsv=args.samples_tsv,
        labels=args.labels,
        obsm_keys=args.obsm_keys,
        annot_cols=annot_cols,
        ref_tag=args.ref_tag,
    )

    adata_filtered = adata[adata.obs["cell2location"] != "Unknown"].copy()
    sc.pp.neighbors(adata_filtered, use_rep="X_pca_harmony")
    sc.tl.umap(adata_filtered)

    output_dir = Path(args.output_dir / "figures_umap")
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"output_dir: {output_dir}")

    for col in DEFAULT_UMAP_COLOR_COLS:
        out_name = f"spatial_{args.output_label}_umap_{col}.png"
        plot_umap_with_count_legend(
            adata_filtered,
            col,
            output_dir / out_name,
            title=f"Harmony UMAP — {col}",
        )


if __name__ == "__main__":
    main()
