#!/usr/bin/env python
"""Summarize cell-type annotations in neighborhoods of GFP+ spots.

Reads Loupe export CSVs from for_loupe/, parses grid coordinates from Visium HD
barcodes (e.g. s_008um_00788_00395-1), finds neighbors within one or more radii,
and reports pooled and per-seed cell-type distributions around GFP+ seeds.

The annotation source is configurable (default: ``cellestimate``); other Loupe
columns such as ``cell2location``, ``predicted.id``, or ``nnls_celltype_assigned``
may be used via ``--annotation_col``.

GFP+ seeds are spots where norm_tlscore_gfp_probe == "norm_tlscore_gfp_probe"
(score > 0.07 in assemble_loupe_file.py).
"""

from __future__ import annotations

import argparse
import math
import re
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_RADII_UM = (8, 16, 24, 48, 80, 120)
DEFAULT_GFP_VALUE = "norm_tlscore_gfp_probe"
DEFAULT_GFP_COL = "norm_tlscore_gfp_probe"
DEFAULT_ANNOTATION_COL = "cellestimate"
DEFAULT_C2L_COL = "cell2location"  # deprecated alias only
DEFAULT_BARCODE_COL = "Barcode"
CELLTYPE_COL = "celltype"
COHORT = "Tuyu_Apr2026"


def annotation_slug(annotation_col: str) -> str:
    """Convert a Loupe annotation column name to a filesystem-safe slug.

    Args:
        annotation_col: Column name (e.g. ``cellestimate``, ``predicted.id``).

    Returns:
        str: Lowercase slug (e.g. ``cellestimate``, ``predicted_id``).
    """
    return re.sub(r"[^A-Za-z0-9]+", "_", annotation_col).strip("_").lower()


def resolve_annotation_col(annotation_col: str, c2l_col: str | None) -> str:
    """Resolve the annotation column, honoring deprecated ``--c2l_col``.

    Args:
        annotation_col: Value from ``--annotation_col``.
        c2l_col: Optional deprecated override from ``--c2l_col``.

    Returns:
        str: Annotation column name to use.
    """
    if c2l_col is not None:
        warnings.warn(
            "--c2l_col is deprecated; use --annotation_col instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return c2l_col
    return annotation_col


def global_result_path(output_dir: Path, label: str, annotation_col: str) -> Path:
    """Path for tissue-wide annotation distribution CSV."""
    slug = annotation_slug(annotation_col)
    return output_dir / f"spatial_{label}_gfp_neighbor_global_{slug}.csv"


def pooled_all_radii_path(output_dir: Path, label: str, annotation_col: str) -> Path:
    """Path for pooled neighbor results across all radii."""
    slug = annotation_slug(annotation_col)
    return output_dir / f"spatial_{label}_gfp_neighbor_{slug}_all_radii.csv"


def per_seed_all_radii_path(output_dir: Path, label: str, annotation_col: str) -> Path:
    """Path for per-seed neighbor summaries across all radii."""
    slug = annotation_slug(annotation_col)
    return output_dir / f"spatial_{label}_gfp_neighbor_{slug}_per_seed_all_radii.csv"


def summary_path(output_dir: Path, annotation_col: str) -> Path:
    """Path for combined summary CSV across samples."""
    slug = annotation_slug(annotation_col)
    return output_dir / f"all_samples_gfp_neighbor_{slug}_summary.csv"


def bar_plot_path(plot_dir: Path, sample: str, radius_um: float, annotation_col: str) -> Path:
    """Path for a single-radius bar plot PNG."""
    slug = annotation_slug(annotation_col)
    return plot_dir / f"spatial_{sample}_gfp_neighbor_{slug}_bar_{int(radius_um)}um.png"


def perm_radius_path(output_dir: Path, label: str, radius_um: float, annotation_col: str) -> Path:
    """Path for permutation results at one radius."""
    slug = annotation_slug(annotation_col)
    return output_dir / f"spatial_{label}_gfp_neighbor_{slug}_perm_{int(radius_um)}um.csv"


def perm_all_radii_path(output_dir: Path, label: str, annotation_col: str) -> Path:
    """Path for combined permutation results across radii for one sample."""
    slug = annotation_slug(annotation_col)
    return output_dir / f"spatial_{label}_gfp_neighbor_{slug}_perm_all_radii.csv"


def perm_summary_path(output_dir: Path, annotation_col: str) -> Path:
    """Path for combined permutation summary across samples."""
    slug = annotation_slug(annotation_col)
    return output_dir / f"all_samples_gfp_neighbor_{slug}_perm_summary.csv"


def build_celltype_color_map(cell_types: list[str] | set[str]) -> dict[str, str | tuple[float, float, float, float]]:
    """Assign stable colors to cell-type labels for plotting.

    Args:
        cell_types: Cell type names to include in the map. Sorted alphabetically
            before colors are assigned so the same label always gets the same
            color across samples and radii.

    Returns:
        dict[str, color]: Mapping from cell-type label to an RGBA color.
    """
    unique = sorted(set(cell_types))
    palette = (
        list(plt.get_cmap("tab20").colors)
        + list(plt.get_cmap("tab20b").colors)
        + list(plt.get_cmap("tab20c").colors)
    )
    return {cell_type: palette[i % len(palette)] for i, cell_type in enumerate(unique)}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the GFP neighbor analysis script.

    Returns:
        argparse.Namespace: Parsed arguments including input/output paths,
            sample filters, radii, GFP/annotation column names, and plot options.
    """
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
        default=repo_root / "results" / COHORT,
        help="Directory for summary tables",
    )
    parser.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="Optional sample labels (filename stems without spatial_/ _for_loupe). "
        "Default: all spatial_*_for_loupe.csv files in input_dir.",
    )
    parser.add_argument(
        "--radii_um",
        type=float,
        nargs="+",
        default=list(DEFAULT_RADII_UM),
        help="Neighborhood radii in microns",
    )
    parser.add_argument(
        "--bin_um",
        type=float,
        default=None,
        help="Bin size in microns. Default: inferred from filename (008um -> 8).",
    )
    parser.add_argument(
        "--gfp_col",
        default=DEFAULT_GFP_COL,
        help="Column used to identify GFP+ seed spots",
    )
    parser.add_argument(
        "--gfp_value",
        default=DEFAULT_GFP_VALUE,
        help="Value in gfp_col that marks GFP+ seed spots",
    )
    parser.add_argument(
        "--annotation_col",
        default=DEFAULT_ANNOTATION_COL,
        help="Loupe column for cell-type annotation (default: cellestimate)",
    )
    parser.add_argument(
        "--c2l_col",
        default=None,
        help="Deprecated alias for --annotation_col (default: cell2location when set)",
    )
    parser.add_argument(
        "--write_per_seed",
        action="store_true",
        help="Also write per-seed neighbor summary tables (one file per radius).",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Write bar plots of top neighbor cell types for each radius.",
    )
    parser.add_argument(
        "--plot_only",
        action="store_true",
        help="Skip analysis and plot from existing *_gfp_neighbor_{slug}_all_radii.csv files.",
    )
    parser.add_argument(
        "--plot_dir",
        type=Path,
        default=None,
        help="Directory for figure output (default: <output_dir>/figures).",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=10,
        help="Number of top cell types to show in bar plots.",
    )
    parser.add_argument(
        "--plot_metric",
        choices=("fraction", "count"),
        default="fraction",
        help="Bar height metric for neighbor cell-type plots.",
    )
    return parser.parse_args()


def write_csv(df: pd.DataFrame, path: Path) -> None:
    """Write a DataFrame to CSV, creating parent directories as needed.

    Args:
        df: Table to write.
        path: Destination CSV path. Parent directories are created if missing.

    Returns:
        None
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def infer_bin_um(path: Path, override: float | None) -> float:
    """Infer Visium HD bin size in microns from a filename.

    Args:
        path: Input file path containing a token like ``_008um_``.
        override: If provided, this value is returned instead of parsing.

    Returns:
        float: Bin size in microns (e.g. 8 for ``008um``).

    Raises:
        ValueError: If bin size cannot be inferred and no override is given.
    """
    if override is not None:
        return override
    match = re.search(r"_(\d{3})um_", path.name)
    if not match:
        raise ValueError(f"Could not infer bin size from filename: {path.name}")
    return int(match.group(1))


def infer_label(path: Path) -> str:
    """Extract the sample label from a Loupe CSV filename.

    Args:
        path: Path to ``spatial_{label}_for_loupe.csv`` or similar.

    Returns:
        str: Sample label (stem between ``spatial_`` and ``_for_loupe.csv``),
            or ``path.stem`` if the pattern does not match.
    """
    name = path.name
    if name.startswith("spatial_") and name.endswith("_for_loupe.csv"):
        return name[len("spatial_") : -len("_for_loupe.csv")]
    return path.stem


def barcode_pattern(bin_um: int) -> re.Pattern[str]:
    """Build a regex that extracts grid coordinates from Visium HD barcodes.

    Args:
        bin_um: Bin size in microns used in the barcode tag (e.g. 8 -> ``008um``).

    Returns:
        re.Pattern[str]: Compiled pattern with two capture groups for x and y
            grid indices (e.g. ``s_008um_00788_00395-1``).
    """
    width = 3 if bin_um < 100 else len(str(int(bin_um)))
    tag = f"{bin_um:0{width}d}um"
    return re.compile(rf"s_{tag}_(\d+)_(\d+)-")


def load_loupe_csv(
    path: Path,
    bin_um: float,
    gfp_col: str,
    gfp_value: str,
    annotation_col: str,
    barcode_col: str,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load a Loupe export CSV and parse spatial coordinates from barcodes.

    Args:
        path: Path to a ``spatial_*_for_loupe.csv`` file.
        bin_um: Bin size in microns for barcode parsing.
        gfp_col: Column name indicating GFP probe score category.
        gfp_value: Value in ``gfp_col`` that marks a GFP+ seed spot.
        annotation_col: Column name for cell-type assignments.
        barcode_col: Column name for spot barcodes.

    Returns:
        tuple containing:
            - pd.DataFrame: Spot metadata (barcode, GFP column, annotation,
              ``is_gfp_seed`` flag).
            - np.ndarray: x grid coordinates (int32), one per spot.
            - np.ndarray: y grid coordinates (int32), one per spot.
            - np.ndarray: cell-type labels (str), one per spot.
            - np.ndarray: Boolean mask, True for GFP+ seed spots.

    Raises:
        ValueError: If required columns are missing from the CSV.
    """
    df = pd.read_csv(path)
    required = {barcode_col, gfp_col, annotation_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")

    pattern = barcode_pattern(int(bin_um))
    coords = df[barcode_col].astype(str).str.extract(pattern)
    valid = coords[0].notna()
    if not valid.all():
        print(f"  Warning: dropping {(~valid).sum()} spots with unparseable barcodes")

    out = df.loc[valid, [barcode_col, gfp_col, annotation_col]].copy()
    out[annotation_col] = out[annotation_col].fillna("Unknown").astype(str)
    out["is_gfp_seed"] = out[gfp_col] == gfp_value

    xs = coords.loc[valid, 0].astype(np.int32).to_numpy()
    ys = coords.loc[valid, 1].astype(np.int32).to_numpy()
    labels = out[annotation_col].to_numpy()
    is_seed = out["is_gfp_seed"].to_numpy()
    return out.reset_index(drop=True), xs, ys, labels, is_seed


def build_spatial_index(xs: np.ndarray, ys: np.ndarray, cell_size: int = 20) -> dict[tuple[int, int], list[int]]:
    """Bucket spot indices into a coarse grid for fast neighbor lookup.

    Args:
        xs: x grid coordinates for all spots (barcode grid units).
        ys: y grid coordinates for all spots (barcode grid units).
        cell_size: Width of each spatial hash bucket in grid units. Larger
            buckets mean fewer buckets but more spots checked per bucket
            during neighbor search. Does not affect which neighbors are found.

    Returns:
        dict[tuple[int, int], list[int]]: Mapping from bucket ``(bx, by)`` to
            list of spot indices whose coordinates fall in that bucket.
    """
    grid: dict[tuple[int, int], list[int]] = defaultdict(list)
    for idx in range(len(xs)):
        grid[(int(xs[idx]) // cell_size, int(ys[idx]) // cell_size)].append(idx)
    return grid


def neighbor_indices(
    xs: np.ndarray,
    ys: np.ndarray,
    grid: dict[tuple[int, int], list[int]],
    seed_idx: int,
    radius_um: float,
    bin_um: float,
    cell_size: int = 20,
) -> list[int]:
    """Return indices of all spots within a Euclidean radius of a seed spot.

    Args:
        xs: x grid coordinates for all spots.
        ys: y grid coordinates for all spots.
        grid: Spatial index from :func:`build_spatial_index`.
        seed_idx: Index of the seed spot whose neighbors are sought.
        radius_um: Search radius in microns.
        bin_um: Bin size in microns (converts radius to grid units).
        cell_size: Bucket width used when ``grid`` was built; must match the
            value passed to :func:`build_spatial_index`.

    Returns:
        list[int]: Spot indices (excluding the seed itself) within ``radius_um``
            microns of the seed location.
    """
    sx, sy = int(xs[seed_idx]), int(ys[seed_idx])
    radius_bins_sq = (radius_um / bin_um) ** 2
    cx0, cy0 = sx // cell_size, sy // cell_size
    reach = int(math.ceil(radius_um / bin_um / cell_size)) + 1
    found: list[int] = []

    for cx in range(cx0 - reach, cx0 + reach + 1):
        for cy in range(cy0 - reach, cy0 + reach + 1):
            for idx in grid.get((cx, cy), []):
                if idx == seed_idx:
                    continue
                dx = int(xs[idx]) - sx
                dy = int(ys[idx]) - sy
                if dx * dx + dy * dy <= radius_bins_sq:
                    found.append(idx)
    return found


def global_distribution(labels: np.ndarray) -> pd.DataFrame:
    """Compute tissue-wide cell-type label counts and fractions.

    Args:
        labels: Array of cell-type assignment strings, one per spot.

    Returns:
        pd.DataFrame: Columns ``celltype``, ``count``, and ``fraction``
            sorted by descending count.
    """
    counts = pd.Series(labels).value_counts()
    total = len(labels)
    out = counts.rename("count").reset_index()
    out.columns = [CELLTYPE_COL, "count"]
    out["fraction"] = out["count"] / total
    return out


def pooled_neighbor_distribution(
    df: pd.DataFrame,
    xs: np.ndarray,
    ys: np.ndarray,
    labels: np.ndarray,
    grid: dict[tuple[int, int], list[int]],
    seed_indices: np.ndarray,
    radius_um: float,
    bin_um: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pool cell-type counts across neighborhoods of GFP+ seed spots.

    Args:
        df: Spot metadata DataFrame from :func:`load_loupe_csv`.
        xs: x grid coordinates for all spots.
        ys: y grid coordinates for all spots.
        labels: Cell-type labels for all spots.
        grid: Spatial index from :func:`build_spatial_index`.
        seed_indices: Indices of GFP+ seed spots.
        radius_um: Neighborhood radius in microns.
        bin_um: Bin size in microns.

    Returns:
        tuple containing:
            - pd.DataFrame: Pooled neighbor counts/fractions with columns
              ``celltype``, ``count``, ``fraction``.
            - pd.DataFrame: Per-seed summary with barcode, coordinates, seed
              cell type, neighbor count, and dominant neighbor label.
    """
    counter: Counter[str] = Counter()
    per_seed_rows: list[dict] = []

    for seed_idx in seed_indices:
        nbr_idxs = neighbor_indices(xs, ys, grid, int(seed_idx), radius_um, bin_um)
        local = Counter(labels[nbr_idxs])
        counter.update(local)
        dominant = local.most_common(1)[0][0] if local else "None"
        per_seed_rows.append(
            {
                DEFAULT_BARCODE_COL: df.at[seed_idx, DEFAULT_BARCODE_COL],
                "seed_x": int(xs[seed_idx]),
                "seed_y": int(ys[seed_idx]),
                "seed_celltype": labels[seed_idx],
                "n_neighbors": len(nbr_idxs),
                "dominant_neighbor_celltype": dominant,
            }
        )

    total = sum(counter.values())
    pooled = pd.DataFrame(
        [
            {CELLTYPE_COL: ct, "count": n, "fraction": n / total if total else 0.0}
            for ct, n in counter.most_common()
        ]
    )
    return pooled, pd.DataFrame(per_seed_rows)


def add_enrichment(neighbor_df: pd.DataFrame, global_df: pd.DataFrame) -> pd.DataFrame:
    """Add global tissue fractions and fold-change vs global to neighbor table.

    Args:
        neighbor_df: Pooled neighbor distribution with ``celltype``,
            ``count``, and ``fraction`` columns.
        global_df: Tissue-wide distribution from :func:`global_distribution`.

    Returns:
        pd.DataFrame: Input columns plus ``global_fraction`` and
            ``fold_change_vs_global`` (neighbor fraction / global fraction).
    """
    merged = neighbor_df.merge(
        global_df[[CELLTYPE_COL, "fraction"]].rename(columns={"fraction": "global_fraction"}),
        on=CELLTYPE_COL,
        how="left",
    )
    merged["global_fraction"] = merged["global_fraction"].fillna(0.0)
    merged["fold_change_vs_global"] = merged.apply(
        lambda row: row["fraction"] / row["global_fraction"]
        if row["global_fraction"] > 0
        else float("inf"),
        axis=1,
    )
    return merged


def infer_label_from_pooled_path(path: Path, annotation_col: str) -> str:
    """Extract sample label from a pooled neighbor results filename.

    Args:
        path: Path to ``spatial_{label}_gfp_neighbor_{slug}_all_radii.csv``.
        annotation_col: Annotation column used to build the filename slug.

    Returns:
        str: Sample label, or ``path.stem`` if the pattern does not match.
    """
    slug = annotation_slug(annotation_col)
    name = path.name
    prefix = "spatial_"
    suffix = f"_gfp_neighbor_{slug}_all_radii.csv"
    if name.startswith(prefix) and name.endswith(suffix):
        return name[len(prefix) : -len(suffix)]
    return path.stem


def discover_pooled_result_files(
    output_dir: Path,
    annotation_col: str,
    samples: list[str] | None,
) -> list[Path]:
    """Find pooled neighbor result CSVs in an output directory.

    Args:
        output_dir: Directory containing pooled result CSVs for the annotation.
        annotation_col: Annotation column used to build the filename slug.
        samples: Optional list of sample labels to restrict discovery. If
            ``None``, all matching files are returned.

    Returns:
        list[Path]: Sorted paths to pooled result CSV files.

    Raises:
        FileNotFoundError: If ``samples`` is set but a requested label has no file.
    """
    slug = annotation_slug(annotation_col)
    all_files = sorted(output_dir.glob(f"spatial_*_gfp_neighbor_{slug}_all_radii.csv"))
    if samples:
        wanted = {s.strip() for s in samples}
        files = [p for p in all_files if infer_label_from_pooled_path(p, annotation_col) in wanted]
        missing = wanted - {infer_label_from_pooled_path(p, annotation_col) for p in files}
        if missing:
            raise FileNotFoundError(
                f"No pooled neighbor CSV found for sample labels: {sorted(missing)}"
            )
        return files
    return all_files


def plot_neighbor_celltypes_bar(
    pooled_df: pd.DataFrame,
    sample: str,
    radius_um: float,
    top_n: int,
    plot_metric: str,
    plot_path: Path,
    annotation_col: str,
    global_df: pd.DataFrame | None = None,
    color_map: dict[str, str | tuple[float, float, float, float]] | None = None,
) -> None:
    """Save a horizontal bar chart of top neighbor cell types.

    Args:
        pooled_df: Pooled neighbor results with ``sample``, ``radius_um``,
            ``celltype``, ``count``, and ``fraction`` columns.
        sample: Sample label to plot.
        radius_um: Neighborhood radius in microns.
        top_n: Number of top cell types to include.
        plot_metric: Bar length column, either ``fraction`` or ``count``.
        plot_path: Output PNG path.
        annotation_col: Loupe annotation column name (for plot title).
        global_df: Optional global tissue distribution; when provided and
            ``plot_metric`` is ``fraction``, global values are overlaid as
            open circles using the same cell-type colors.
        color_map: Optional mapping from cell-type label to color. Built
            from the plotted labels if ``None``.

    Returns:
        None
    """
    sub = pooled_df[(pooled_df["sample"] == sample) & (pooled_df["radius_um"] == radius_um)].copy()
    if sub.empty:
        print(f"  Warning: no rows for {sample} at radius={radius_um:g} µm; skipping plot")
        return

    sub = sub.sort_values("count", ascending=False).head(top_n)
    metric_col = plot_metric
    y_labels = sub[CELLTYPE_COL].tolist()
    neighbor_vals = sub[metric_col].to_numpy()
    if color_map is None:
        color_map = build_celltype_color_map(y_labels)
    bar_colors = [color_map.get(ct, "#AAAAAA") for ct in y_labels]

    fig_h = max(3.5, 0.35 * len(y_labels) + 1.5)
    fig, ax = plt.subplots(figsize=(8, fig_h))
    y_pos = np.arange(len(y_labels))
    ax.barh(y_pos, neighbor_vals, color=bar_colors)

    if global_df is not None and plot_metric == "fraction":
        global_sub = global_df[global_df[CELLTYPE_COL].isin(y_labels)].set_index(CELLTYPE_COL)
        global_vals = [global_sub.loc[ct, "fraction"] if ct in global_sub.index else 0.0 for ct in y_labels]
        for val, y, ct in zip(global_vals, y_pos, y_labels):
            ax.scatter(
                val,
                y,
                facecolors="white",
                edgecolors=color_map.get(ct, "#AAAAAA"),
                linewidths=1.5,
                s=55,
                zorder=3,
            )
        ax.scatter([], [], facecolors="white", edgecolors="black", linewidths=1.5, s=55, label="Global tissue")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels)
    ax.invert_yaxis()
    if plot_metric == "fraction":
        ax.set_xlabel("Fraction of neighbor spots")
        ax.set_xlim(0, min(1.0, max(neighbor_vals.max() * 1.15, 0.05)))
    else:
        ax.set_xlabel("Neighbor spot count")
    ax.set_title(
        f"{sample}\nTop {annotation_col} neighbors within {radius_um:g} µm of GFP+ spots"
    )
    if global_df is not None and plot_metric == "fraction":
        ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_neighbor_celltypes_all_radii(
    pooled_df: pd.DataFrame,
    sample: str,
    radii_um: list[float],
    top_n: int,
    plot_metric: str,
    plot_dir: Path,
    annotation_col: str,
    global_df: pd.DataFrame | None = None,
    color_map: dict[str, str | tuple[float, float, float, float]] | None = None,
) -> list[Path]:
    """Write bar plots for one sample across multiple radii.

    Args:
        pooled_df: Pooled neighbor results table.
        sample: Sample label to plot.
        radii_um: List of neighborhood radii in microns.
        top_n: Number of top cell types per plot.
        plot_metric: Bar length column (``fraction`` or ``count``).
        plot_dir: Directory for output PNG files.
        annotation_col: Loupe annotation column name (for filenames and titles).
        global_df: Optional global tissue distribution for overlay.
        color_map: Optional stable cell-type color mapping shared across radii.

    Returns:
        list[Path]: Paths to PNG files that were successfully written.
    """
    saved: list[Path] = []
    for radius_um in radii_um:
        plot_path = bar_plot_path(plot_dir, sample, radius_um, annotation_col)
        plot_neighbor_celltypes_bar(
            pooled_df,
            sample,
            radius_um,
            top_n,
            plot_metric,
            plot_path,
            annotation_col,
            global_df=global_df,
            color_map=color_map,
        )
        if plot_path.exists():
            saved.append(plot_path)
            print(f"  Saved {plot_path.name}")
    return saved


def plot_from_pooled_results(
    output_dir: Path,
    plot_dir: Path,
    annotation_col: str,
    samples: list[str] | None,
    radii_um: list[float] | None,
    top_n: int,
    plot_metric: str,
) -> None:
    """Generate bar plots from existing pooled neighbor CSV files.

    Args:
        output_dir: Directory containing pooled and global result CSVs.
        plot_dir: Directory where PNG figures are written.
        annotation_col: Annotation column used to locate result files.
        samples: Optional sample labels to plot; ``None`` plots all found.
        radii_um: Radii to plot; ``None`` uses all radii present in each CSV.
        top_n: Number of top cell types per bar chart.
        plot_metric: Bar length column (``fraction`` or ``count``).

    Returns:
        None

    Raises:
        FileNotFoundError: If no pooled result CSVs are found.
    """
    pooled_files = discover_pooled_result_files(output_dir, annotation_col, samples)
    slug = annotation_slug(annotation_col)
    if not pooled_files:
        raise FileNotFoundError(
            f"No spatial_*_gfp_neighbor_{slug}_all_radii.csv files found in {output_dir}"
        )

    print(f"Plotting from {len(pooled_files)} pooled result file(s) ({annotation_col})")
    plot_dir.mkdir(parents=True, exist_ok=True)

    all_cell_types: set[str] = set()
    pooled_by_label: dict[str, pd.DataFrame] = {}
    for pooled_path in pooled_files:
        pooled_df = pd.read_csv(pooled_path)
        label = infer_label_from_pooled_path(pooled_path, annotation_col)
        pooled_by_label[label] = pooled_df
        all_cell_types.update(pooled_df[CELLTYPE_COL].astype(str))
    color_map = build_celltype_color_map(all_cell_types)

    for label, pooled_df in pooled_by_label.items():
        global_path = global_result_path(output_dir, label, annotation_col)
        global_df = pd.read_csv(global_path) if global_path.exists() else None

        sample_radii = radii_um
        if sample_radii is None:
            sample_radii = sorted(pooled_df["radius_um"].unique())

        print(f"\nSample: {label}")
        plot_neighbor_celltypes_all_radii(
            pooled_df,
            label,
            sample_radii,
            top_n,
            plot_metric,
            plot_dir,
            annotation_col,
            global_df=global_df,
            color_map=color_map,
        )


def resolve_loupe_csv(input_dir: Path, sample: str) -> Path:
    """Resolve a sample label to a single Loupe export CSV path.

    Prefers the canonical ``spatial_{sample}_for_loupe.csv`` file. Ref-tagged
    variants (e.g. ``_compv1_mmv2``) are only returned when the sample name
    explicitly includes that suffix.

    Args:
        input_dir: Directory containing ``spatial_*_for_loupe.csv`` files.
        sample: Sample label (stem between ``spatial_`` and ``_for_loupe.csv``).

    Returns:
        Path: Resolved Loupe CSV path.

    Raises:
        FileNotFoundError: If no matching file exists, or only ref-tagged
            variants exist for a base sample name.
    """
    sample = sample.strip()
    canonical = input_dir / f"spatial_{sample}_for_loupe.csv"
    if canonical.is_file():
        return canonical

    exact = [p for p in input_dir.glob("spatial_*_for_loupe.csv") if infer_label(p) == sample]
    if len(exact) == 1:
        return exact[0]

    variants = sorted(
        p for p in input_dir.glob("spatial_*_for_loupe.csv")
        if infer_label(p).startswith(f"{sample}_")
    )
    if variants:
        raise FileNotFoundError(
            f"No base Loupe CSV for {sample!r}. Expected {canonical.name}. "
            f"Ref-tagged variants exist but were not selected: {[p.name for p in variants]}"
        )
    raise FileNotFoundError(f"No Loupe CSV found for sample {sample!r} in {input_dir}")


def discover_input_files(input_dir: Path, samples: list[str] | None) -> list[Path]:
    """Find Loupe export CSVs in an input directory.

    Args:
        input_dir: Directory containing ``spatial_*_for_loupe.csv`` files.
        samples: Optional sample labels to restrict discovery. If ``None``, all
            matching files are returned.

    Returns:
        list[Path]: Sorted paths to Loupe CSV files.

    Raises:
        FileNotFoundError: If ``samples`` is set but a requested label has no file.
    """
    if samples:
        files = [resolve_loupe_csv(input_dir, s) for s in samples]
        for path in files:
            print(f"Using Loupe input: {path}")
        return files
    return sorted(input_dir.glob("spatial_*_for_loupe.csv"))


def analyze_sample(
    path: Path,
    output_dir: Path,
    radii_um: list[float],
    bin_um_override: float | None,
    gfp_col: str,
    gfp_value: str,
    annotation_col: str,
    write_per_seed: bool,
) -> list[dict]:
    """Run GFP neighborhood analysis for one sample and write result CSVs.

    Args:
        path: Path to a ``spatial_*_for_loupe.csv`` file.
        output_dir: Directory for global, pooled, and optional per-seed outputs.
        radii_um: Neighborhood radii in microns to evaluate.
        bin_um_override: Optional bin size override; inferred from filename if
            ``None``.
        gfp_col: Column marking GFP probe status.
        gfp_value: Value in ``gfp_col`` that defines GFP+ seeds.
        annotation_col: Loupe column for cell-type annotation.
        write_per_seed: If ``True``, also write per-seed neighbor summaries.

    Returns:
        list[dict]: One summary dict per radius with sample metadata and top
            neighbor cell type statistics.
    """
    label = infer_label(path)
    bin_um = infer_bin_um(path, bin_um_override)
    print(f"\nSample: {label} ({path.name}, bin={bin_um:g} µm, annotation={annotation_col})")

    df, xs, ys, labels, is_seed = load_loupe_csv(
        path, bin_um, gfp_col, gfp_value, annotation_col, DEFAULT_BARCODE_COL
    )
    grid = build_spatial_index(xs, ys)
    seed_indices = np.flatnonzero(is_seed)
    print(f"  Spots: {len(df):,}; GFP+ seeds: {len(seed_indices):,}")

    global_df = global_distribution(labels)
    global_df.insert(0, "sample", label)
    global_df.insert(1, "annotation_col", annotation_col)

    summary_rows: list[dict] = []
    pooled_all: list[pd.DataFrame] = []
    per_seed_all: list[pd.DataFrame] = []
    for radius_um in radii_um:
        pooled, per_seed = pooled_neighbor_distribution(
            df, xs, ys, labels, grid, seed_indices, radius_um, bin_um
        )
        pooled = add_enrichment(pooled, global_df.drop(columns=["sample", "annotation_col"]))
        pooled.insert(0, "radius_um", radius_um)
        pooled.insert(0, "sample", label)
        pooled.insert(1, "annotation_col", annotation_col)
        pooled_all.append(pooled)

        per_seed.insert(0, "radius_um", radius_um)
        per_seed.insert(0, "sample", label)
        per_seed.insert(1, "annotation_col", annotation_col)
        per_seed_all.append(per_seed)

        print(
            f"  radius={radius_um:g} µm: mean neighbors/seed="
            f"{per_seed['n_neighbors'].mean():.1f}"
        )

        top = pooled.sort_values("count", ascending=False).head(1)
        summary_rows.append(
            {
                "sample": label,
                "annotation_col": annotation_col,
                "radius_um": radius_um,
                "n_spots": len(df),
                "n_gfp_seeds": len(seed_indices),
                "mean_neighbors_per_seed": per_seed["n_neighbors"].mean(),
                "top_neighbor_celltype": top.iloc[0][CELLTYPE_COL] if len(top) else None,
                "top_neighbor_fraction": top.iloc[0]["fraction"] if len(top) else None,
                "top_neighbor_fold_change": top.iloc[0]["fold_change_vs_global"] if len(top) else None,
            }
        )

    global_path = global_result_path(output_dir, label, annotation_col)
    pooled_path = pooled_all_radii_path(output_dir, label, annotation_col)
    write_csv(global_df, global_path)
    write_csv(pd.concat(pooled_all, ignore_index=True), pooled_path)
    print(f"  Saved {global_path.name}, {pooled_path.name}")

    if write_per_seed:
        per_seed_path = per_seed_all_radii_path(output_dir, label, annotation_col)
        write_csv(pd.concat(per_seed_all, ignore_index=True), per_seed_path)
        print(f"  Saved {per_seed_path.name}")

    return summary_rows


def main() -> None:
    """Entry point: analyze GFP neighborhoods and optionally plot bar charts.

    Parses CLI arguments, discovers input files, runs :func:`analyze_sample`
    for each, writes a combined summary CSV, and optionally generates plots
    from pooled results.

    Returns:
        None
    """
    args = parse_args()
    annotation_col = resolve_annotation_col(args.annotation_col, args.c2l_col)
    plot_dir = args.plot_dir or (args.output_dir / "figures")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.plot_only:
        if not args.plot:
            args.plot = True
        plot_from_pooled_results(
            args.output_dir,
            plot_dir,
            annotation_col,
            args.samples,
            args.radii_um,
            args.top_n,
            args.plot_metric,
        )
        return

    input_files = discover_input_files(args.input_dir, args.samples)
    if not input_files:
        raise FileNotFoundError(f"No spatial_*_for_loupe.csv files found in {args.input_dir}")

    print(f"Found {len(input_files)} Loupe CSV(s) in {args.input_dir}")
    print(f"Annotation column: {annotation_col}")
    all_summary: list[dict] = []
    for path in input_files:
        all_summary.extend(
            analyze_sample(
                path,
                args.output_dir,
                args.radii_um,
                args.bin_um,
                args.gfp_col,
                args.gfp_value,
                annotation_col,
                args.write_per_seed,
            )
        )

    summary_df = pd.DataFrame(all_summary)
    combined_summary_path = summary_path(args.output_dir, annotation_col)
    write_csv(summary_df, combined_summary_path)
    print(f"\nSaved combined summary: {combined_summary_path}")
    print(summary_df.to_string(index=False))

    if args.plot:
        plot_from_pooled_results(
            args.output_dir,
            plot_dir,
            annotation_col,
            args.samples,
            args.radii_um,
            args.top_n,
            args.plot_metric,
        )


if __name__ == "__main__":
    main()
