#!/usr/bin/env python
"""Add ``cest_w_tumor`` and ``cest_w_tumor_rabies`` to spatial h5ad ``obs``.

``cest_w_tumor``: fusion score above threshold → tumor label, else cellestimate.
``cest_w_tumor_rabies``: rabies-positive bins → ``rabies``; else same as ``cest_w_tumor``.

Input should include ``rabies_assignment`` (e.g. ``*_rabies_assignment.h5ad``).
Writes a new h5ad by default (suffix ``_cest_w_tumor``).

Example::

    python add_cest_w_tumor_column.py \\
      --h5ad ../data/processed/Tuyu_Apr2026/spatial_*_rabies_assignment.h5ad \\
      --tumor_threshold 0.1
"""

from __future__ import annotations

import argparse
import warnings
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import scanpy as sc

DEFAULT_SCORE_COL = "tltumscore_C11orf95_relafus1"
DEFAULT_CELLESTIMATE_COL = "cellestimate"
DEFAULT_TUMOR_LABEL = "tumor"
DEFAULT_RABIES_COL = "rabies_assignment"
DEFAULT_RABIES_POSITIVE = "positive"
DEFAULT_RABIES_LABEL = "rabies"
OUTPUT_COLUMN = "cest_w_tumor"
RABIES_COLUMN = "cest_w_tumor_rabies"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed paths, thresholds, and column names.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--h5ad", type=Path, required=True, help="Input spatial h5ad")
    parser.add_argument(
        "--output_h5ad",
        type=Path,
        default=None,
        help=f"Output h5ad (default: {{stem}}_{OUTPUT_COLUMN}.h5ad)",
    )
    parser.add_argument(
        "--overwrite-inplace",
        action="store_true",
        help="Write annotated h5ad back to --h5ad instead of a new file",
    )
    parser.add_argument(
        "--tumor_threshold",
        type=float,
        default=0.1,
        help="Fusion score strictly above this value labels a bin as tumor",
    )
    parser.add_argument(
        "--score_col",
        default=DEFAULT_SCORE_COL,
        help=f"Fusion score column in adata.obs (default {DEFAULT_SCORE_COL})",
    )
    parser.add_argument(
        "--cellestimate_col",
        default=DEFAULT_CELLESTIMATE_COL,
        help=f"Cellestimate column when score is at or below threshold (default {DEFAULT_CELLESTIMATE_COL})",
    )
    parser.add_argument(
        "--tumor_label",
        default=DEFAULT_TUMOR_LABEL,
        help=f"Label for bins above threshold (default {DEFAULT_TUMOR_LABEL})",
    )
    parser.add_argument(
        "--rabies_col",
        default=DEFAULT_RABIES_COL,
        help=f"Rabies assignment column (default {DEFAULT_RABIES_COL})",
    )
    parser.add_argument(
        "--rabies_positive",
        default=DEFAULT_RABIES_POSITIVE,
        help=f"Value in --rabies_col that maps to --rabies_label (default {DEFAULT_RABIES_POSITIVE})",
    )
    parser.add_argument(
        "--rabies_label",
        default=DEFAULT_RABIES_LABEL,
        help=f"Label when --rabies_col equals --rabies_positive (default {DEFAULT_RABIES_LABEL})",
    )
    return parser.parse_args()


def build_cest_w_tumor(
    score: pd.Series,
    cellestimate: pd.Series,
    threshold: float,
    tumor_label: str = DEFAULT_TUMOR_LABEL,
) -> pd.Series:
    """Build per-bin ``cest_w_tumor`` labels from fusion score and cellestimate.

    Args:
        score: Per-bin fusion score (e.g. ``tltumscore_C11orf95_relafus1``).
        cellestimate: Per-bin cellestimate label when not tumor.
        threshold: Bins with ``score > threshold`` become ``tumor_label``.
        tumor_label: String assigned to tumor bins.

    Returns:
        pd.Series: Combined labels aligned to ``score.index``.
    """
    score_f = pd.to_numeric(score, errors="coerce")
    is_tumor = score_f > threshold
    base = cellestimate.astype(str)
    out = base.where(~is_tumor, tumor_label)
    return out.astype(object)


def build_cest_w_tumor_rabies(
    rabies_assignment: pd.Series,
    cest_w_tumor: pd.Series,
    positive_value: str = DEFAULT_RABIES_POSITIVE,
    rabies_label: str = DEFAULT_RABIES_LABEL,
) -> pd.Series:
    """Build ``cest_w_tumor_rabies``: rabies-positive bins override with ``rabies_label``.

    Args:
        rabies_assignment: Per-bin rabies call (e.g. ``positive`` / ``negative``).
        cest_w_tumor: Per-bin tumor/cellestimate label when not rabies-positive.
        positive_value: ``rabies_assignment`` value treated as rabies-positive.
        rabies_label: Label assigned when rabies-positive.

    Returns:
        pd.Series: Labels aligned to ``rabies_assignment.index``.
    """
    is_rabies = rabies_assignment.astype(str) == positive_value
    return cest_w_tumor.astype(str).where(~is_rabies, rabies_label).astype(object)


def resolve_output_path(
    h5ad_path: Path,
    output_h5ad: Path | None,
    overwrite_inplace: bool,
) -> Path:
    """Resolve the output h5ad path from CLI flags."""
    if overwrite_inplace:
        return h5ad_path
    if output_h5ad is not None:
        return output_h5ad.expanduser()
    return h5ad_path.parent / f"{h5ad_path.stem}_{OUTPUT_COLUMN}.h5ad"


def main() -> None:
    """Load h5ad, add ``cest_w_tumor`` and ``cest_w_tumor_rabies``, write object."""
    args = parse_args()
    h5ad_path = args.h5ad.expanduser().resolve()
    if not h5ad_path.is_file():
        raise FileNotFoundError(f"h5ad not found: {h5ad_path}")

    adata = sc.read_h5ad(h5ad_path)
    for col in (args.score_col, args.cellestimate_col, args.rabies_col):
        if col not in adata.obs.columns:
            preview = ", ".join(list(adata.obs.columns.astype(str))[:25])
            raise KeyError(
                f"Missing obs column {col!r}. Available (first 25): {preview}"
            )

    score = adata.obs[args.score_col]
    cellestimate = adata.obs[args.cellestimate_col]
    rabies = adata.obs[args.rabies_col]
    n_nan_score = int(pd.to_numeric(score, errors="coerce").isna().sum())
    if n_nan_score > 0:
        warnings.warn(
            f"{n_nan_score} bins have non-numeric/NaN {args.score_col!r}; "
            f"they will not be labeled tumor.",
            stacklevel=2,
        )

    for col_name in (OUTPUT_COLUMN, RABIES_COLUMN):
        if col_name in adata.obs.columns:
            warnings.warn(
                f"Overwriting existing obs column {col_name!r}.",
                stacklevel=2,
            )

    cest_w_tumor = build_cest_w_tumor(
        score,
        cellestimate,
        args.tumor_threshold,
        tumor_label=args.tumor_label,
    )
    adata.obs[OUTPUT_COLUMN] = pd.Series(
        cest_w_tumor.to_numpy(dtype=object),
        index=adata.obs_names,
        dtype=object,
    )

    cest_w_tumor_rabies = build_cest_w_tumor_rabies(
        rabies,
        cest_w_tumor,
        positive_value=args.rabies_positive,
        rabies_label=args.rabies_label,
    )
    adata.obs[RABIES_COLUMN] = pd.Series(
        cest_w_tumor_rabies.to_numpy(dtype=object),
        index=adata.obs_names,
        dtype=object,
    )

    score_f = pd.to_numeric(score, errors="coerce")
    n_tumor = int((score_f > args.tumor_threshold).sum())
    n_rabies = int((rabies.astype(str) == args.rabies_positive).sum())
    n_total = int(adata.n_obs)
    if n_rabies == 0:
        warnings.warn(
            f"No bins with {args.rabies_col}={args.rabies_positive!r}.",
            stacklevel=2,
        )

    adata.uns[OUTPUT_COLUMN] = {
        "score_col": args.score_col,
        "cellestimate_col": args.cellestimate_col,
        "tumor_threshold": float(args.tumor_threshold),
        "tumor_label": args.tumor_label,
        "rule": "score > threshold -> tumor_label; else cellestimate_col",
        "n_tumor": n_tumor,
        "n_total": n_total,
        "n_nan_score": n_nan_score,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    adata.uns[RABIES_COLUMN] = {
        "rabies_col": args.rabies_col,
        "rabies_positive": args.rabies_positive,
        "rabies_label": args.rabies_label,
        "depends_on": OUTPUT_COLUMN,
        "rule": "rabies_positive -> rabies_label; else cest_w_tumor",
        "n_rabies": n_rabies,
        "n_total": n_total,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }

    write_path = resolve_output_path(
        h5ad_path, args.output_h5ad, args.overwrite_inplace
    )
    write_path = write_path.resolve()
    write_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(write_path)

    print(f"{OUTPUT_COLUMN}: {n_tumor} tumor / {n_total} bins ({100.0 * n_tumor / n_total:.2f}%)")
    print(f"{RABIES_COLUMN}: {n_rabies} rabies / {n_total} bins ({100.0 * n_rabies / n_total:.2f}%)")
    print(f"threshold={args.tumor_threshold} on {args.score_col!r}")
    print(f"{RABIES_COLUMN} value counts:")
    print(adata.obs[RABIES_COLUMN].astype(str).value_counts().to_string())
    print(f"Wrote: {write_path}")


if __name__ == "__main__":
    main()
