#!/usr/bin/env python
"""Score ZR fusion on preprocess h5ad and subset to non-tumor spots.

Loads the base preprocess AnnData, applies the 93-gene C11orf95-RELA fusion
signature score, keeps spots at or below a configurable threshold, and writes a
filtered h5ad plus an audit report and spatial QC plot.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "scripts"))

from score_gaba_glut_signatures import (  # noqa: E402
    SIGNATURE_SETS,
    score_signature_set,
)

COHORT = "Tuyu_Apr2026"
DEFAULT_MAX_SCORE = 0.0
ZR_FUSION_KEY = "zr_fusion"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for non-tumor spot filtering.

    Returns:
        argparse.Namespace: Parsed input paths, threshold, and output options.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=None,
        help="Base preprocess h5ad (spatial_{label}.h5ad)",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Sample label for default path lookup and output filenames",
    )
    parser.add_argument(
        "--cohort",
        default=COHORT,
        help="Cohort subdirectory under data/processed/",
    )
    parser.add_argument(
        "--max_score",
        type=float,
        default=DEFAULT_MAX_SCORE,
        help="Keep spots with score_93_gene_zr_fusion <= this value (default 0)",
    )
    parser.add_argument(
        "--zr_fusion_genes",
        type=Path,
        default=_REPO_ROOT / "reference" / SIGNATURE_SETS[ZR_FUSION_KEY][0],
        help="93-gene ZR fusion reference list",
    )
    parser.add_argument(
        "--out_h5ad",
        type=Path,
        default=None,
        help="Filtered h5ad output (default spatial_{label}_nontumor.h5ad)",
    )
    parser.add_argument(
        "--report_csv",
        type=Path,
        default=None,
        help="Per-spot filter report CSV",
    )
    parser.add_argument(
        "--qc_plot",
        type=Path,
        default=None,
        help="Spatial QC PNG (kept vs dropped)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="QC plot resolution",
    )
    return parser.parse_args()


def resolve_input_h5ad(args: argparse.Namespace) -> Path:
    """Resolve base preprocess h5ad from CLI arguments.

    Args:
        args: Parsed namespace with ``h5ad``, ``label``, and ``cohort``.

    Returns:
        Path: Existing preprocess h5ad file.

    Raises:
        FileNotFoundError: If the h5ad file does not exist.
        ValueError: If neither ``--h5ad`` nor ``--label`` is provided.
    """
    if args.h5ad is not None:
        h5ad = args.h5ad.expanduser().resolve()
        if not h5ad.is_file():
            raise FileNotFoundError(f"h5ad not found: {h5ad}")
        return h5ad

    if not args.label:
        raise ValueError("Provide --h5ad or --label")

    processed_dir = _REPO_ROOT / "data" / "processed" / args.cohort
    candidate = processed_dir / f"spatial_{args.label}.h5ad"
    if not candidate.is_file():
        raise FileNotFoundError(f"Preprocess h5ad not found: {candidate}")
    return candidate


def default_outputs(
    args: argparse.Namespace,
    label: str,
) -> tuple[Path, Path, Path]:
    """Build default output paths when not set on the CLI.

    Args:
        args: Parsed namespace with ``cohort``.
        label: Sample label.

    Returns:
        tuple[Path, Path, Path]: ``out_h5ad``, ``report_csv``, ``qc_plot`` paths.
    """
    processed = _REPO_ROOT / "data" / "processed" / args.cohort
    results = _REPO_ROOT / "results" / args.cohort
    out_h5ad = args.out_h5ad or (processed / f"spatial_{label}_nontumor.h5ad")
    report_csv = args.report_csv or (
        processed / f"spatial_{label}_nontumor_filter_report.csv"
    )
    qc_plot = args.qc_plot or (
        results / f"spatial_{label}_zr_fusion_nontumor_qc.png"
    )
    return out_h5ad, report_csv, qc_plot


def plot_filter_qc(
    coords: np.ndarray,
    kept: np.ndarray,
    out_path: Path,
    title: str,
    dpi: int = 200,
) -> None:
    """Save a spatial scatter of kept vs dropped spots.

    Args:
        coords: Spot coordinates shaped ``(n_spots, 2)``.
        kept: Boolean mask; True = kept (non-tumor).
        out_path: Output PNG path.
        title: Figure title.
        dpi: PNG resolution.

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    dropped = ~kept
    if dropped.any():
        ax.scatter(
            coords[dropped, 0],
            coords[dropped, 1],
            c="#cccccc",
            s=1.5,
            alpha=0.5,
            linewidths=0,
            label="dropped (tumor)",
        )
    if kept.any():
        ax.scatter(
            coords[kept, 0],
            coords[kept, 1],
            c="#2166ac",
            s=1.5,
            alpha=0.8,
            linewidths=0,
            label="kept (non-tumor)",
        )
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title)
    ax.legend(loc="upper right", markerscale=4, frameon=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def filter_nontumor_spots(
    adata: sc.AnnData,
    gene_path: Path,
    max_score: float,
) -> tuple[sc.AnnData, pd.DataFrame, str]:
    """Score ZR fusion and subset to non-tumor spots.

    Args:
        adata: Preprocess AnnData with log-normalized ``.X``.
        gene_path: Path to the 93-gene fusion reference list.
        max_score: Keep spots with fusion score <= this value.

    Returns:
        tuple[sc.AnnData, pd.DataFrame, str]: Filtered adata, report frame, score column.

    Raises:
        RuntimeError: If scoring fails (too few matched genes).
    """
    _, score_col = SIGNATURE_SETS[ZR_FUSION_KEY]
    result = score_signature_set(
        adata,
        set_key=ZR_FUSION_KEY,
        score_col=score_col,
        gene_path=gene_path,
    )
    if not result.scored:
        raise RuntimeError(
            f"ZR fusion scoring failed: only {len(result.matched_genes)} genes matched "
            f"(need >= 3). Missing: {result.missing_genes[:10]}..."
        )

    scores = adata.obs[score_col].astype(float)
    kept = scores <= max_score
    report = pd.DataFrame(
        {
            "barcode": adata.obs_names,
            "score_93_gene_zr_fusion": scores.values,
            "max_score_threshold": max_score,
            "is_nontumor": kept.values,
        }
    )
    filtered = adata[kept].copy()
    return filtered, report, score_col


def main() -> None:
    """Score ZR fusion, subset spots, and write filtered h5ad plus QC outputs.

    Returns:
        None
    """
    args = parse_args()
    h5ad_path = resolve_input_h5ad(args)
    label = args.label or h5ad_path.stem.replace("spatial_", "").replace("_nontumor", "")
    out_h5ad, report_csv, qc_plot = default_outputs(args, label)

    print(f"Loading {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    n_before = adata.n_obs

    filtered, report, score_col = filter_nontumor_spots(
        adata,
        gene_path=args.zr_fusion_genes,
        max_score=args.max_score,
    )
    n_kept = filtered.n_obs
    n_dropped = n_before - n_kept
    print(
        f"Threshold: {score_col} <= {args.max_score}; "
        f"kept {n_kept}/{n_before} spots ({n_dropped} dropped)"
    )

    if "spatial" not in filtered.obsm:
        raise ValueError(f"{h5ad_path} missing obsm['spatial']")

    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    filtered.write_h5ad(out_h5ad)
    print(f"Saved: {out_h5ad}")

    report_csv.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(report_csv, index=False)
    print(f"Saved: {report_csv}")

    coords = adata.obsm["spatial"]
    kept_mask = report["is_nontumor"].to_numpy()
    plot_filter_qc(
        coords,
        kept_mask,
        qc_plot,
        title=(
            f"{label}\n{score_col} <= {args.max_score}: "
            f"{n_kept} kept / {n_dropped} dropped"
        ),
        dpi=args.dpi,
    )
    print(f"Saved: {qc_plot}")


if __name__ == "__main__":
    main()
