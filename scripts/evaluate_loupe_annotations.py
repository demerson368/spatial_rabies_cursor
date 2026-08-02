#!/usr/bin/env python
"""Evaluate Loupe annotation exports for spatial coherence and cross-method agreement.

Reads ``for_loupe/spatial_*_for_loupe.csv``, parses Visium HD grid coordinates from
barcodes, and computes per-sample:

  - Moran's I and Geary's C on one-hot label indicators (spatial clustering)
  - Neighborhood enrichment z-scores (squidpy)
  - Mean neighbor label consistency (fraction of spatial neighbors sharing label)

Optionally reports pairwise agreement (fraction, Cohen's kappa, ARI, NMI) between
annotation columns.
"""

from __future__ import annotations

import argparse
import re
import sys
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq
from sklearn.metrics import adjusted_rand_score, cohen_kappa_score, normalized_mutual_info_score

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from gfp_neighbor_cell2location import discover_input_files, infer_bin_um, infer_label

COHORT = "Tuyu_Apr2026"
DEFAULT_BARCODE_COL = "Barcode"
DEFAULT_LABEL_COLS = [
    "tlcelltype_assigned",
    "celltype_assigned_norm",
    "nnls_celltype_assigned",
    "cellestimate",
    "predicted.id",
    "cell2location",
]
PREFIXES_TO_STRIP = ("tlscore_", "nnls_")


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input_dir",
        type=Path,
        default=repo_root / "for_loupe",
        help="Directory containing spatial_*_for_loupe.csv files",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=repo_root / "results" / COHORT / "loupe_eval",
        help="Directory for evaluation tables and figures",
    )
    parser.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="Optional sample labels (stems between spatial_ and _for_loupe). "
        "Default: all matching CSVs in input_dir.",
    )
    parser.add_argument(
        "--label_cols",
        nargs="+",
        default=list(DEFAULT_LABEL_COLS),
        help="Annotation columns to evaluate spatially",
    )
    parser.add_argument(
        "--bin_um",
        type=float,
        default=None,
        help="Bin size in microns. Default: inferred from filename (008um -> 8).",
    )
    parser.add_argument(
        "--n_neighs",
        type=int,
        default=6,
        help="Grid neighbors for spatial graph (Visium HD default: 6)",
    )
    parser.add_argument(
        "--n_perms",
        type=int,
        default=99,
        help="Permutations for Moran's I, Geary's C, and nhood enrichment",
    )
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=4,
        help="Parallel jobs for squidpy spatial statistics",
    )
    parser.add_argument(
        "--min_label_spots",
        type=int,
        default=50,
        help="Drop label categories with fewer spots before autocorrelation / enrichment",
    )
    parser.add_argument(
        "--max_labels",
        type=int,
        default=40,
        help="Keep at most this many label categories (by spot count) for autocorrelation",
    )
    parser.add_argument(
        "--subsample",
        type=int,
        default=0,
        help="If >0, randomly subsample this many spots (for faster exploratory runs)",
    )
    parser.add_argument(
        "--agreement",
        action="store_true",
        help="Compute pairwise agreement metrics between label columns",
    )
    parser.add_argument(
        "--no_normalize_for_agreement",
        action="store_false",
        dest="normalize_for_agreement",
        help="Compare raw label strings for agreement metrics (default strips tlscore_/nnls_ prefixes)",
    )
    parser.set_defaults(normalize_for_agreement=True)
    parser.add_argument(
        "--skip_autocorr",
        action="store_true",
        help="Skip Moran's I and Geary's C (faster)",
    )
    parser.add_argument(
        "--skip_enrichment",
        action="store_true",
        help="Skip neighborhood enrichment (faster)",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Write heatmaps for neighborhood enrichment and coherence bar chart",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for subsampling",
    )
    return parser.parse_args()


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def barcode_pattern(bin_um: int) -> re.Pattern[str]:
    width = 3 if bin_um < 100 else len(str(int(bin_um)))
    tag = f"{bin_um:0{width}d}um"
    return re.compile(rf"s_{tag}_(\d+)_(\d+)-")


def load_loupe_annotations(
    path: Path,
    bin_um: float,
    label_cols: list[str],
    barcode_col: str = DEFAULT_BARCODE_COL,
    subsample: int = 0,
    seed: int = 0,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Load Loupe CSV, parse grid coords, and retain requested label columns."""
    df = pd.read_csv(path)
    if barcode_col not in df.columns:
        raise ValueError(f"{path}: missing column {barcode_col!r}")

    present_cols = [c for c in label_cols if c in df.columns]
    missing_cols = sorted(set(label_cols) - set(present_cols))
    if missing_cols:
        print(f"  Warning: {path.name} missing label columns: {missing_cols}")
    if not present_cols:
        raise ValueError(f"{path}: none of the requested label columns are present")

    pattern = barcode_pattern(int(bin_um))
    coords = df[barcode_col].astype(str).str.extract(pattern)
    valid = coords[0].notna()
    if not valid.all():
        print(f"  Warning: dropping {(~valid).sum()} spots with unparseable barcodes")

    keep = [barcode_col, *present_cols]
    out = df.loc[valid, keep].copy()
    for col in present_cols:
        out[col] = out[col].fillna("Unknown").astype(str)

    xs = coords.loc[valid, 0].astype(np.int32).to_numpy()
    ys = coords.loc[valid, 1].astype(np.int32).to_numpy()

    if subsample > 0 and len(out) > subsample:
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(out), size=subsample, replace=False)
        out = out.iloc[idx].reset_index(drop=True)
        xs = xs[idx]
        ys = ys[idx]
        print(f"  Subsampled to {subsample:,} spots")

    return out.reset_index(drop=True), xs, ys


def filter_label_categories(
    labels: pd.Series,
    min_spots: int,
    max_labels: int,
) -> tuple[pd.Series, list[str]]:
    """Keep frequent categories; rare types become ``Other_rare``."""
    counts = labels.value_counts()
    keep = counts[counts >= min_spots].index.tolist()
    if max_labels > 0 and len(keep) > max_labels:
        keep = counts.head(max_labels).index.tolist()
    keep_set = set(keep)
    filtered = labels.where(labels.isin(keep_set), other="Other_rare")
    return filtered.astype(str), sorted(keep_set)


def build_spatial_adata(
    obs: pd.DataFrame,
    xs: np.ndarray,
    ys: np.ndarray,
    label_col: str,
    n_neighs: int,
) -> ad.AnnData:
    adata = ad.AnnData(
        X=np.zeros((obs.shape[0], 1), dtype=np.float32),
        obs=obs[[label_col]].copy(),
    )
    adata.obs[label_col] = adata.obs[label_col].astype("category")
    adata.obsm["spatial"] = np.column_stack([xs.astype(float), ys.astype(float)])
    sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=n_neighs)
    return adata


def neighbor_coherence(adata: ad.AnnData, label_col: str) -> float:
    """Mean fraction of spatial neighbors sharing the same label."""
    conn = adata.obsp["spatial_connectivities"]
    labels = adata.obs[label_col].astype(str).to_numpy()
    scores: list[float] = []
    for i in range(adata.n_obs):
        nbrs = conn[i].nonzero()[1]
        if len(nbrs) == 0:
            continue
        scores.append(float((labels[nbrs] == labels[i]).mean()))
    return float(np.mean(scores)) if scores else float("nan")


def label_onehot_autocorr(
    adata: ad.AnnData,
    label_col: str,
    min_spots: int,
    max_labels: int,
    n_perms: int,
    n_jobs: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Moran's I and Geary's C on one-hot indicators for each label category."""
    filtered, _ = filter_label_categories(
        adata.obs[label_col].astype(str),
        min_spots=min_spots,
        max_labels=max_labels,
    )
    dummies = pd.get_dummies(filtered, dtype=np.float32)
    if dummies.shape[1] == 0:
        raise ValueError(f"No label categories passed filters for {label_col}")

    dummy = ad.AnnData(
        X=dummies.to_numpy(),
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=dummies.columns),
    )
    dummy.obsp = adata.obsp.copy()
    dummy.obsm = adata.obsm.copy()
    dummy.uns = adata.uns.copy()

    sq.gr.spatial_autocorr(dummy, mode="moran", n_perms=n_perms, n_jobs=n_jobs)
    moran = dummy.uns["moranI"].copy()
    moran.index.name = "label_category"
    moran = moran.reset_index()

    sq.gr.spatial_autocorr(dummy, mode="geary", n_perms=n_perms, n_jobs=n_jobs)
    geary = dummy.uns["gearyC"].copy()
    geary.index.name = "label_category"
    geary = geary.reset_index()

    return moran, geary


def nhood_enrichment_matrix(
    adata: ad.AnnData,
    label_col: str,
    min_spots: int,
    max_labels: int,
    n_perms: int,
    n_jobs: int,
) -> pd.DataFrame:
    filtered, _ = filter_label_categories(
        adata.obs[label_col].astype(str),
        min_spots=min_spots,
        max_labels=max_labels,
    )
    eval_obs = adata.obs.copy()
    eval_obs[label_col] = filtered.astype("category")
    eval_adata = ad.AnnData(
        X=np.zeros((adata.n_obs, 1), dtype=np.float32),
        obs=eval_obs[[label_col]].copy(),
    )
    eval_adata.obsp = adata.obsp.copy()
    eval_adata.obsm = adata.obsm.copy()
    eval_adata.uns = adata.uns.copy()
    eval_adata.obs[label_col] = eval_adata.obs[label_col].astype("category")

    sq.gr.nhood_enrichment(
        eval_adata,
        cluster_key=label_col,
        n_perms=n_perms,
        n_jobs=n_jobs,
        seed=0,
    )
    key = f"{label_col}_nhood_enrichment"
    zscore = eval_adata.uns[key]["zscore"]
    categories = list(eval_adata.obs[label_col].cat.categories)
    return pd.DataFrame(zscore, index=categories, columns=categories)


def strip_prefix(label: str) -> str:
    for prefix in PREFIXES_TO_STRIP:
        if label.startswith(prefix):
            return label[len(prefix) :]
    return label


def pairwise_agreement(
    df: pd.DataFrame,
    cols: list[str],
    normalize: bool,
) -> pd.DataFrame:
    rows: list[dict] = []
    for col_a, col_b in combinations(cols, 2):
        a = df[col_a].astype(str)
        b = df[col_b].astype(str)
        if normalize:
            a = a.map(strip_prefix)
            b = b.map(strip_prefix)
        rows.append(
            {
                "col_a": col_a,
                "col_b": col_b,
                "n_spots": len(df),
                "agreement": float((a == b).mean()),
                "cohen_kappa": float(cohen_kappa_score(a, b)),
                "ari": float(adjusted_rand_score(a, b)),
                "nmi": float(normalized_mutual_info_score(a, b)),
            }
        )
    return pd.DataFrame(rows)


def plot_enrichment_heatmap(zscore: pd.DataFrame, title: str, out_path: Path, max_types: int = 25) -> None:
    if zscore.empty:
        return
    if len(zscore) > max_types:
        top = zscore.index[:max_types]
        zscore = zscore.loc[top, top]

    fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(zscore)), max(5, 0.35 * len(zscore))))
    im = ax.imshow(zscore.to_numpy(), cmap="coolwarm", vmin=-3, vmax=3, aspect="auto")
    ax.set_xticks(range(len(zscore.columns)))
    ax.set_yticks(range(len(zscore.index)))
    ax.set_xticklabels(zscore.columns, rotation=90, fontsize=7)
    ax.set_yticklabels(zscore.index, fontsize=7)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="Z-score")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_coherence_bar(coherence_df: pd.DataFrame, label: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4))
    order = coherence_df.sort_values("neighbor_coherence", ascending=False)
    ax.bar(order["annotation_col"], order["neighbor_coherence"], color="steelblue")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Neighbor label coherence")
    ax.set_title(f"Spatial label coherence — {label}")
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def analyze_sample(
    path: Path,
    output_dir: Path,
    label_cols: list[str],
    bin_um_override: float | None,
    n_neighs: int,
    n_perms: int,
    n_jobs: int,
    min_label_spots: int,
    max_labels: int,
    subsample: int,
    seed: int,
    agreement: bool,
    normalize_for_agreement: bool,
    skip_autocorr: bool,
    skip_enrichment: bool,
    plot: bool,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    label = infer_label(path)
    bin_um = infer_bin_um(path, bin_um_override)
    print(f"\nSample: {label} ({path.name}, bin={bin_um:g} µm)")

    df, xs, ys = load_loupe_annotations(
        path,
        bin_um=bin_um,
        label_cols=label_cols,
        subsample=subsample,
        seed=seed,
    )
    present_cols = [c for c in label_cols if c in df.columns]
    print(f"  Spots: {len(df):,}; evaluating columns: {present_cols}")

    coherence_rows: list[dict] = []
    autocorr_frames: list[pd.DataFrame] = []
    geary_frames: list[pd.DataFrame] = []

    plot_dir = output_dir / "figures"

    for col in present_cols:
        adata = build_spatial_adata(df, xs, ys, col, n_neighs=n_neighs)
        coherence = neighbor_coherence(adata, col)
        n_types = df[col].nunique()
        coherence_rows.append(
            {
                "sample": label,
                "annotation_col": col,
                "n_spots": len(df),
                "n_label_categories": int(n_types),
                "neighbor_coherence": coherence,
            }
        )
        print(f"  {col}: coherence={coherence:.3f} ({n_types} categories)")

        if not skip_autocorr:
            moran, geary = label_onehot_autocorr(
                adata,
                col,
                min_spots=min_label_spots,
                max_labels=max_labels,
                n_perms=n_perms,
                n_jobs=n_jobs,
            )
            moran.insert(0, "annotation_col", col)
            moran.insert(0, "sample", label)
            geary.insert(0, "annotation_col", col)
            geary.insert(0, "sample", label)
            autocorr_frames.append(moran)
            geary_frames.append(geary)

        if not skip_enrichment:
            zscore = nhood_enrichment_matrix(
                adata,
                col,
                min_spots=min_label_spots,
                max_labels=max_labels,
                n_perms=n_perms,
                n_jobs=n_jobs,
            )
            enrich_path = output_dir / f"spatial_{label}_{col}_nhood_enrichment.csv"
            zscore.to_csv(enrich_path)
            print(f"  Saved {enrich_path.name}")

            if plot:
                plot_enrichment_heatmap(
                    zscore,
                    title=f"{label} — {col} nhood enrichment",
                    out_path=plot_dir / f"spatial_{label}_{col}_nhood_enrichment.png",
                )

    coherence_df = pd.DataFrame(coherence_rows)
    write_csv(coherence_df, output_dir / f"spatial_{label}_neighbor_coherence.csv")

    if autocorr_frames:
        moran_all = pd.concat(autocorr_frames, ignore_index=True)
        geary_all = pd.concat(geary_frames, ignore_index=True)
        write_csv(moran_all, output_dir / f"spatial_{label}_label_moranI.csv")
        write_csv(geary_all, output_dir / f"spatial_{label}_label_gearyC.csv")
        print(f"  Saved Moran's I / Geary's C tables")

    if plot:
        plot_coherence_bar(
            coherence_df,
            label=label,
            out_path=plot_dir / f"spatial_{label}_neighbor_coherence.png",
        )

    agreement_df = None
    if agreement and len(present_cols) >= 2:
        agreement_df = pairwise_agreement(df, present_cols, normalize=normalize_for_agreement)
        agreement_df.insert(0, "sample", label)
        write_csv(agreement_df, output_dir / f"spatial_{label}_annotation_agreement.csv")
        print(f"  Saved annotation agreement table")

    return coherence_df, agreement_df


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    input_files = discover_input_files(args.input_dir, args.samples)
    if not input_files:
        raise FileNotFoundError(f"No spatial_*_for_loupe.csv files found in {args.input_dir}")

    print(f"Found {len(input_files)} Loupe CSV(s) in {args.input_dir}")

    all_coherence: list[pd.DataFrame] = []
    all_agreement: list[pd.DataFrame] = []

    for path in input_files:
        coherence_df, agreement_df = analyze_sample(
            path=path,
            output_dir=args.output_dir,
            label_cols=args.label_cols,
            bin_um_override=args.bin_um,
            n_neighs=args.n_neighs,
            n_perms=args.n_perms,
            n_jobs=args.n_jobs,
            min_label_spots=args.min_label_spots,
            max_labels=args.max_labels,
            subsample=args.subsample,
            seed=args.seed,
            agreement=args.agreement,
            normalize_for_agreement=args.normalize_for_agreement,
            skip_autocorr=args.skip_autocorr,
            skip_enrichment=args.skip_enrichment,
            plot=args.plot,
        )
        all_coherence.append(coherence_df)
        if agreement_df is not None:
            all_agreement.append(agreement_df)

    write_csv(pd.concat(all_coherence, ignore_index=True), args.output_dir / "all_samples_neighbor_coherence.csv")
    print(f"\nSaved combined coherence: {args.output_dir / 'all_samples_neighbor_coherence.csv'}")

    if all_agreement:
        agreement_all = pd.concat(all_agreement, ignore_index=True)
        write_csv(agreement_all, args.output_dir / "all_samples_annotation_agreement.csv")
        print(f"Saved combined agreement: {args.output_dir / 'all_samples_annotation_agreement.csv'}")
        print(agreement_all.to_string(index=False))


if __name__ == "__main__":
    main()
