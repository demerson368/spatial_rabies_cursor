#!/usr/bin/env python
"""Interactive HTML viewer: gene/fusion score with cellestimate highlight.

Loads an extended spatial h5ad, merges ``cellestimate`` labels from a
``for_loupe/*.csv`` file, scores a gene or signature set (reusing
``score_gaba_glut_signatures.py``), and writes a self-contained Plotly HTML
with a dropdown to highlight bins by cell type. Layers (bottom to top): IF
image, faded cellestimate highlight borders, gene/fusion score scatter.

Requires: pip install plotly

Example::

    python plot_score_celltype_html.py \\
      --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um \\
      --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/ \\
      --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif \\
      --gene Ptprc \\
      --loupe_csv for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_for_loupe.csv \\
      --highlight_celltype Glut
"""

from __future__ import annotations

import argparse
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from PIL import Image
from skimage import transform

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "scripts"))

from score_gaba_glut_signatures import (  # noqa: E402
    gene_obs_slug,
    infer_label_from_h5ad,
    parse_gene_args,
    prepare_if_background,
    resolve_h5ad_path,
    resolve_if_image_path,
    score_display_mask,
    score_signature_set,
    score_single_gene,
)

SIGNATURE_SETS: dict[str, tuple[str, str]] = {
    "gaba_glut_intersection": (
        "gaba_glut_intersection_genes_allen_yao_et_al.txt",
        "score_gaba_glut_intersection_allen",
    ),
    "gaba_unique": (
        "gaba_unique_genes_allen_yao_et_al.txt",
        "score_gaba_unique_allen",
    ),
    "glut_unique": (
        "glut_unique_genes_allen_yao_et_al.txt",
        "score_glut_unique_allen",
    ),
    "zr_fusion": (
        "93_core_gene_set_commonly_upregulated_in_human_and_mouse_C11orf95_RELAFUS1_driven_ependymoma.txt",
        "score_93_gene_zr_fusion",
    ),

}
COHORT = "Tuyu_Apr2026"
DEFAULT_IF_DIR = (
    "/research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/"
    "common/Tuyu_Apr2026/IF_images"
)
DEFAULT_BARCODE_COL = "Barcode"
DEFAULT_CELLTYPE_COL = "cellestimate"
DEFAULT_SCORE_COLORMAP = "Purples" #"magma" #"Yellows" #"Greens"
CELLTYPE_COL_ALIASES = {"cell_estimate": "cellestimate"}
UBCLASS_COLS: dict[str, str] = {
    "IT-ET Glut": "#F21D1D",
    "NP-CT-L6b Glut": "#F23F1D",
    "OB-CR Glut": "#F2601D",
    "DG-IMN Glut": "#F2821D",
    "OB-IMN GABA": "#F2A41D",
    "CTX-CGE GABA": "#F2C51D",
    "CTX-MGE GABA": "#F2E71D",
    "CNU-MGE GABA": "#DCF21D",
    "CNU-LGE GABA": "#BAF21D",
    "CNU-HYa GABA": "#98F21D",
    "LSX GABA": "#77F21D",
    "HY GABA": "#55F21D",
    "CNU-HYa Glut": "#34F21D",
    "HY Glut": "#1DF228",
    "HY Gnrh1 Glut": "#1DF24A",
    "HY MM Glut": "#1DF26C",
    "MH-LH Glut": "#1DF28D",
    "TH Glut": "#1DF2AF",
    "MB Glut": "#1DF2D1",
    "MB GABA": "#1DF2F2",
    "MB Dopa": "#1DD1F2",
    "MB-HB Sero": "#AFAFAF",
    "P Glut": "#1D8DF2",
    "MY Glut": "#1D6CF2",
    "Pineal Glut": "#1D4AF2",
    "P GABA": "#1D28F2",
    "MY GABA": "#341DF2",
    "CB GABA": "#551DF2",
    "CB Glut": "#771DF2",
    "Astro-Epen": "#981DF2",
    "OPC-Oligo": "#BA1DF2",
    "OEC": "#ff7034",
    "Vascular": "#964B00",
    "Microglia NN": "#F21DC5",
    "BAM NN": "#000000",
    "Monocytes NN": "#F21D82",
    "DC NN": "#F21D60",
    "Lymphoid NN": "#F21D3F",
    "GABA": "#1D28F2",
    "Glut": "#34F21D",
}
UBCLASS_FALLBACK_COLOR = "#808080"


def build_ubclass_color_map(
    cell_types: list[str] | set[str],
) -> dict[str, tuple[float, float, float]]:
    """Assign Allen subclass colors to cell-type labels for plotting.

    Args:
        cell_types: Cell type names to include in the map.

    Returns:
        dict[str, tuple[float, float, float]]: Label to RGB tuple in [0, 1].
    """
    color_map: dict[str, tuple[float, float, float]] = {}
    for label in sorted(set(cell_types)):
        hex_color = UBCLASS_COLS.get(label, UBCLASS_FALLBACK_COLOR)
        if label not in UBCLASS_COLS:
            warnings.warn(
                f"No UBCLASS_COLS entry for {label!r}; using gray",
                stacklevel=2,
            )
        color_map[label] = to_rgb(hex_color)
    return color_map


@dataclass(frozen=True)
class ScoreOverlay:
    """Resolved score column and metadata for plotting."""

    score_col: str
    score_label: str
    score_slug: str
    vmin: float
    vmax: float


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the interactive HTML viewer.

    Returns:
        argparse.Namespace: Parsed paths, scoring mode, and plot options.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--h5ad", type=Path, default=None, help="Extended h5ad path")
    parser.add_argument("--label", default=None, help="Sample label for defaults")
    parser.add_argument("--cohort", default=COHORT, help="Cohort under data/processed/")
    parser.add_argument(
        "--loupe_csv",
        type=Path,
        default=None,
        help="Loupe export CSV with cellestimate labels (auto-resolve if omitted)",
    )
    parser.add_argument(
        "--barcode_col",
        default=DEFAULT_BARCODE_COL,
        help="Barcode column in Loupe CSV",
    )
    parser.add_argument(
        "--celltype_col",
        default=DEFAULT_CELLTYPE_COL,
        help="Cell-type column in Loupe CSV (default: cellestimate)",
    )
    parser.add_argument(
        "--highlight_celltype",
        default=None,
        help="Initial cell type to highlight (e.g. Glut)",
    )
    parser.add_argument(
        "--highlight_opacity",
        type=float,
        default=0.45,
        help="Opacity of cellestimate highlight borders (0-1)",
    )
    parser.add_argument("--if_image", type=Path, default=None, help="IF TIFF/PNG path")
    parser.add_argument(
        "--if_dir",
        type=Path,
        default=Path(DEFAULT_IF_DIR),
        help="IF image directory",
    )
    parser.add_argument(
        "--if_channel",
        choices=("red", "green", "gray", "rgb"),
        default="rgb",
        help="IF channel to display (default rgb for full-color panel)",
    )
    parser.add_argument(
        "--if_flip_y",
        action="store_true",
        help="Flip IF vertically before overlay",
    )
    parser.add_argument(
        "--if_flip_x",
        action="store_true",
        help="Flip IF horizontally before overlay",
    )
    parser.add_argument(
        "--if_image_html",
        type=Path,
        default=None,
        help="Optional lower-res IF PNG/TIF for HTML background only",
    )
    parser.add_argument(
        "--if_downsample",
        type=float,
        default=0.25,
        help="Downsample factor for HTML IF layer only (1.0=full, 0.25=quarter side)",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=None,
        help="Output directory (default results/{cohort})",
    )
    parser.add_argument(
        "--signature",
        choices=tuple(SIGNATURE_SETS.keys()),
        default=None,
        help="Signature set to score (default zr_fusion when --gene omitted)",
    )
    parser.add_argument(
        "--gene",
        action="append",
        default=None,
        help="Gene symbol to plot (exclusive with --signature)",
    )
    parser.add_argument(
        "--expression_layer",
        choices=("normalized", "counts"),
        default="normalized",
        help="Expression layer for --gene mode",
    )
    parser.add_argument("--img_key", default="hires", help="Spatial image key in h5ad")
    parser.add_argument("--spot_size", type=float, default=4.0, help="Base marker size")
    parser.add_argument(
        "--score_opacity",
        type=float,
        default=0.8,
        help="Opacity of score markers so IF shows through (0-1)",
    )
    parser.add_argument(
        "--score_vmin",
        type=float,
        default=None,
        help="Minimum score for colormap; spots below this value are hidden (default: data min)",
    )
    parser.add_argument(
        "--score_vmax",
        type=float,
        default=None,
        help="Maximum score for colormap; spots above this value are hidden (default: data max)",
    )
    parser.add_argument(
        "--use_preprocess_h5ad",
        action="store_true",
        help="Resolve spatial_{label}.h5ad instead of _ext.h5ad",
    )
    return parser.parse_args()


def resolve_celltype_col(loupe: pd.DataFrame, requested: str) -> str:
    """Resolve the Loupe cell-type column name, honoring aliases.

    Args:
        loupe: Loupe export DataFrame.
        requested: CLI ``--celltype_col`` value.

    Returns:
        str: Resolved column name present in ``loupe``.

    Raises:
        ValueError: If no matching column exists.
    """
    canonical = CELLTYPE_COL_ALIASES.get(requested.lower(), requested)
    if canonical in loupe.columns:
        return canonical
    if requested in loupe.columns:
        return requested
    lower_map = {str(c).lower(): str(c) for c in loupe.columns}
    if canonical.lower() in lower_map:
        return lower_map[canonical.lower()]
    available = ", ".join(str(c) for c in loupe.columns)
    raise ValueError(
        f"Cell-type column {requested!r} not found in Loupe CSV. "
        f"Available: {available}"
    )


def resolve_loupe_csv_path(args: argparse.Namespace, label: str) -> Path:
    """Resolve Loupe CSV from CLI or auto-glob under for_loupe/.

    Args:
        args: Parsed CLI namespace.
        label: Sample label.

    Returns:
        Path: Existing Loupe CSV file.

    Raises:
        FileNotFoundError: If no CSV can be resolved.
    """
    if args.loupe_csv is not None:
        loupe = args.loupe_csv.expanduser().resolve()
        if not loupe.is_file():
            raise FileNotFoundError(f"Loupe CSV not found: {loupe}")
        return loupe

    loupe_dir = _REPO_ROOT / "for_loupe"
    if not loupe_dir.is_dir():
        raise FileNotFoundError(
            f"--loupe_csv not provided and {loupe_dir} missing"
        )

    direct = loupe_dir / f"spatial_{label}_for_loupe.csv"
    if direct.is_file():
        return direct

    matches = sorted(loupe_dir.glob(f"spatial_{label}*_for_loupe.csv"))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        names = ", ".join(p.name for p in matches)
        raise FileNotFoundError(
            f"Multiple Loupe CSV files match label {label!r}: {names}. "
            "Pass --loupe_csv explicitly."
        )

    raise FileNotFoundError(
        f"No Loupe CSV found for label {label!r} under {loupe_dir}"
    )


def merge_loupe_celltypes(
    adata: sc.AnnData,
    loupe_csv: Path,
    barcode_col: str,
    celltype_col: str,
) -> tuple[sc.AnnData, str]:
    """Attach Loupe cell-type labels to adata and subset to overlapping barcodes.

    Args:
        adata: Spatial AnnData object.
        loupe_csv: Loupe export CSV path.
        barcode_col: Barcode column name.
        celltype_col: Requested cell-type column (may be aliased).

    Returns:
        tuple[sc.AnnData, str]: Subset AnnData with ``obs['_cellestimate_plot']`` and
            the resolved column name.

    Raises:
        ValueError: If barcodes do not overlap or the column is missing.
    """
    loupe = pd.read_csv(loupe_csv)
    resolved_col = resolve_celltype_col(loupe, celltype_col)

    if barcode_col not in loupe.columns:
        raise ValueError(f"Column {barcode_col!r} not found in {loupe_csv}")

    loupe = loupe.set_index(barcode_col)
    overlap = adata.obs_names.intersection(loupe.index)
    if len(overlap) == 0:
        raise ValueError(f"No overlapping barcodes between {loupe_csv} and h5ad")

    adata = adata[overlap].copy()
    labels = loupe.loc[adata.obs_names, resolved_col].astype(str)
    labels = labels.replace({"nan": "", "None": "", "NaN": ""})
    adata.obs["_cellestimate_plot"] = labels
    return adata, resolved_col


def signature_gene_path(set_key: str) -> Path:
    """Resolve the reference gene list path for a signature set key.

    Args:
        set_key: Key in ``SIGNATURE_SETS`` (e.g. ``zr_fusion``).

    Returns:
        Path: Path to the reference gene list under ``reference/``.
    """
    filename = SIGNATURE_SETS[set_key][0]
    return _REPO_ROOT / "reference" / filename


def compute_score_overlay(
    adata: sc.AnnData,
    args: argparse.Namespace,
) -> ScoreOverlay:
    """Score a gene or signature set and return overlay metadata.

    Args:
        adata: Spatial AnnData object (modified in place).
        args: Parsed CLI namespace.

    Returns:
        ScoreOverlay: Score column name, label, slug, and value range.

    Raises:
        ValueError: If both or neither of ``--gene`` and ``--signature`` are set,
            or scoring fails.
    """
    genes = parse_gene_args(args.gene)
    if genes and args.signature:
        raise ValueError("Pass either --gene or --signature, not both")
    if genes:
        if len(genes) > 1:
            warnings.warn(
                f"Multiple --gene values provided; using first only: {genes[0]}",
                stacklevel=2,
            )
        result = score_single_gene(adata, genes[0], args.expression_layer)
        if not result.plotted or result.matched_var is None:
            raise ValueError(f"Gene {genes[0]!r} not found in adata.var_names")
        values = adata.obs[result.obs_col].to_numpy(dtype=float)
        vmin = args.score_vmin if args.score_vmin is not None else float(np.nanmin(values))
        vmax = args.score_vmax if args.score_vmax is not None else float(np.nanmax(values))
        return ScoreOverlay(
            score_col=result.obs_col,
            score_label=f"{result.matched_var} ({args.expression_layer})",
            score_slug=f"gene_{gene_obs_slug(result.matched_var)}",
            vmin=vmin,
            vmax=vmax,
        )

    set_key = args.signature or "zr_fusion"
    _, score_col = SIGNATURE_SETS[set_key]
    result = score_signature_set(
        adata,
        set_key=set_key,
        score_col=score_col,
        gene_path=signature_gene_path(set_key),
    )
    if not result.scored:
        raise ValueError(
            f"Signature {set_key!r} could not be scored "
            f"({len(result.matched_genes)} genes matched)"
        )
    values = adata.obs[score_col].to_numpy(dtype=float)
    vmin = args.score_vmin if args.score_vmin is not None else float(np.nanmin(values))
    vmax = args.score_vmax if args.score_vmax is not None else float(np.nanmax(values))
    return ScoreOverlay(
        score_col=score_col,
        score_label=set_key.replace("_", " "),
        score_slug=set_key,
        vmin=vmin,
        vmax=vmax,
    )


def apply_if_orientation(if_bg: np.ndarray, flip_x: bool, flip_y: bool) -> np.ndarray:
    """Apply optional horizontal/vertical flips to an IF background array.

    Args:
        if_bg: IF background array.
        flip_x: Flip left-right when True.
        flip_y: Flip top-bottom when True.

    Returns:
        np.ndarray: Oriented IF array.
    """
    oriented = if_bg
    if flip_y:
        oriented = np.flipud(oriented)
    if flip_x:
        oriented = np.fliplr(oriented)
    return oriented


def downsample_if_array(if_bg: np.ndarray, factor: float) -> np.ndarray:
    """Downsample an IF background array for lighter HTML embedding.

    Args:
        if_bg: IF background array in [0, 1], 2D or HxWxC.
        factor: Scale factor in (0, 1]; 1.0 returns the input unchanged.

    Returns:
        np.ndarray: Downsampled array with the same channel layout.
    """
    if factor >= 1.0:
        return if_bg

    height, width = if_bg.shape[:2]
    target_hw = (
        max(1, int(round(height * factor))),
        max(1, int(round(width * factor))),
    )
    if if_bg.ndim == 2:
        resized = transform.resize(
            if_bg,
            target_hw,
            preserve_range=True,
            anti_aliasing=True,
        )
    else:
        resized = transform.resize(
            if_bg,
            target_hw + (if_bg.shape[-1],),
            preserve_range=True,
            anti_aliasing=True,
        )
    return np.clip(resized, 0.0, 1.0)


def prepare_html_if_background(
    if_path: Path,
    hires_shape: tuple[int, ...],
    channel: str,
    downsample: float,
    flip_x: bool,
    flip_y: bool,
) -> np.ndarray:
    """Load, align, orient, and downsample an IF array for HTML display.

    Args:
        if_path: IF image path.
        hires_shape: Target hires canvas shape from spatial metadata.
        channel: IF channel selection passed to ``prepare_if_background``.
        downsample: Downsample factor for HTML embedding.
        flip_x: Flip horizontally after alignment.
        flip_y: Flip vertically after alignment.

    Returns:
        np.ndarray: Display-ready IF background array.
    """
    if_bg = prepare_if_background(if_path, hires_shape, channel)
    if_bg = apply_if_orientation(if_bg, flip_x=flip_x, flip_y=flip_y)
    return downsample_if_array(if_bg, downsample)


def resolve_if_html_path(args: argparse.Namespace, label: str, default_if_path: Path) -> Path:
    """Resolve the IF image path used for the HTML background layer.

    Args:
        args: Parsed CLI namespace.
        label: Sample label for auto-discovery fallback.
        default_if_path: Primary ``--if_image`` path when ``--if_image_html`` omitted.

    Returns:
        Path: Existing IF image file for HTML background preparation.
    """
    if args.if_image_html is None:
        return default_if_path

    raw = args.if_image_html.expanduser()
    candidates: list[Path] = [raw]
    if not raw.is_absolute():
        candidates.append(args.if_dir / raw)
        if raw.name == str(raw):
            candidates.insert(0, args.if_dir / raw.name)
    for candidate in candidates:
        resolved = candidate.expanduser().resolve()
        if resolved.is_file():
            return resolved
    tried = ", ".join(str(c) for c in candidates)
    raise FileNotFoundError(
        f"HTML IF image not found: {args.if_image_html}. Tried: {tried}"
    )


def plotly_if_image_kwargs(
    full_height: int,
    full_width: int,
    disp_height: int,
    disp_width: int,
) -> dict[str, float]:
    """Return Plotly ``go.Image`` kwargs scaled to a full coordinate extent.

    Args:
        full_height: Full hires canvas height used for axis limits.
        full_width: Full hires canvas width used for axis limits.
        disp_height: Display raster height embedded in the HTML.
        disp_width: Display raster width embedded in the HTML.

    Returns:
        dict[str, float]: Keyword args for ``x0``, ``y0``, ``dx``, and ``dy``.
    """
    return dict(
        x0=0,
        y0=0,
        dx=full_width / disp_width,
        dy=full_height / disp_height,
    )


def spatial_pixel_coords(adata: sc.AnnData, img_key: str) -> np.ndarray:
    """Return spot coordinates scaled to hires image pixel space.

    Args:
        adata: AnnData with spatial metadata.
        img_key: Image resolution key.

    Returns:
        np.ndarray: Coordinates shaped (n_spots, 2).
    """
    library_id = list(adata.uns["spatial"].keys())[0]
    scale = adata.uns["spatial"][library_id]["scalefactors"][f"tissue_{img_key}_scalef"]
    return np.asarray(adata.obsm["spatial"]) * scale


def rgba_to_plotly(color: tuple[float, ...] | str, alpha: float = 1.0) -> str:
    """Convert a matplotlib RGB or RGBA tuple to a Plotly rgba() string.

    Args:
        color: RGB or RGBA tuple in [0, 1], or an existing color string.
        alpha: Alpha override in [0, 1]; ignored for string colors.

    Returns:
        str: Plotly-compatible color string.
    """
    if isinstance(color, str):
        return color
    if len(color) == 3:
        r, g, b = color
    elif len(color) == 4:
        r, g, b, _ = color
    else:
        raise ValueError(f"Expected RGB or RGBA tuple, got length {len(color)}")
    return f"rgba({int(r * 255)},{int(g * 255)},{int(b * 255)},{alpha:.3f})"


def color_to_matplotlib_rgba(
    color: tuple[float, ...],
    alpha: float,
) -> tuple[float, float, float, float]:
    """Convert a matplotlib RGB or RGBA tuple to an RGBA tuple.

    Args:
        color: RGB or RGBA tuple in [0, 1].
        alpha: Alpha value in [0, 1] applied to the output.

    Returns:
        tuple[float, float, float, float]: RGBA values in [0, 1].
    """
    if len(color) == 3:
        r, g, b = color
    elif len(color) == 4:
        r, g, b, _ = color
    else:
        raise ValueError(f"Expected RGB or RGBA tuple, got length {len(color)}")
    return (float(r), float(g), float(b), float(alpha))


def if_array_to_uint8_rgb(if_bg: np.ndarray) -> np.ndarray:
    """Convert an IF background array to uint8 RGB for image display.

    Args:
        if_bg: Background array in [0, 1], 2D or HxWxC.

    Returns:
        np.ndarray: uint8 array shaped (H, W, 3).
    """
    if if_bg.ndim == 2:
        gray = (np.clip(if_bg, 0.0, 1.0) * 255).astype(np.uint8)
        return np.stack([gray, gray, gray], axis=-1)
    rgb = (np.clip(if_bg[..., :3], 0.0, 1.0) * 255).astype(np.uint8)
    if rgb.shape[-1] == 1:
        return np.concatenate([rgb, rgb, rgb], axis=-1)
    return rgb


def save_if_png(if_bg: np.ndarray, out_path: Path) -> None:
    """Save the aligned IF background as a PNG file.

    Args:
        if_bg: Background array in [0, 1], 2D or HxWxC.
        out_path: Output PNG path.

    Returns:
        None
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(if_array_to_uint8_rgb(if_bg)).save(out_path)


def output_if_png_path(output_dir: Path, label: str, score_slug: str) -> Path:
    """Build the IF PNG output path.

    Args:
        output_dir: Results directory.
        label: Sample label.
        score_slug: Filesystem slug for the score overlay.

    Returns:
        Path: Output IF PNG path.
    """
    return output_dir / f"spatial_{label}_{score_slug}_if_background.png"


def output_if_html_png_path(output_dir: Path, label: str, score_slug: str) -> Path:
    """Build the downsampled IF PNG path embedded in the HTML viewer.

    Args:
        output_dir: Results directory.
        label: Sample label.
        score_slug: Filesystem slug for the score overlay.

    Returns:
        Path: Output HTML IF PNG path.
    """
    return output_dir / f"spatial_{label}_{score_slug}_if_background_html.png"


def output_align_check_path(output_dir: Path, label: str, score_slug: str) -> Path:
    """Build the matplotlib alignment sanity-check PNG path.

    Args:
        output_dir: Results directory.
        label: Sample label.
        score_slug: Filesystem slug for the score overlay.

    Returns:
        Path: Output alignment-check PNG path.
    """
    return output_dir / f"spatial_{label}_{score_slug}_cellestimate_align_check.png"


def output_html_path(output_dir: Path, label: str, score_slug: str) -> Path:
    """Build the HTML output path.

    Args:
        output_dir: Results directory.
        label: Sample label.
        score_slug: Filesystem slug for the score overlay.

    Returns:
        Path: Output HTML path.
    """
    return output_dir / f"spatial_{label}_{score_slug}_cellestimate.html"


def output_csv_path(output_dir: Path, label: str, score_slug: str) -> Path:
    """Build the companion spot-level CSV path.

    Args:
        output_dir: Results directory.
        label: Sample label.
        score_slug: Filesystem slug for the score overlay.

    Returns:
        Path: Output CSV path.
    """
    return output_dir / f"spatial_{label}_{score_slug}_cellestimate_spots.csv"


def write_spot_table(
    adata: sc.AnnData,
    coords: np.ndarray,
    score_col: str,
    celltype_col: str,
    out_path: Path,
) -> None:
    """Write per-spot score and cell-type table for debugging.

    Args:
        adata: AnnData with score and cell-type columns.
        coords: Pixel coordinates (n_spots, 2).
        score_col: Score column in ``adata.obs``.
        celltype_col: Resolved Loupe cell-type column name.
        out_path: Output CSV path.

    Returns:
        None
    """
    df = pd.DataFrame(
        {
            "Barcode": adata.obs_names,
            "score": adata.obs[score_col].to_numpy(),
            celltype_col: adata.obs["_cellestimate_plot"].astype(str),
            "x": coords[:, 0],
            "y": coords[:, 1],
        }
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)


def plot_score_celltype_align_check(
    adata: sc.AnnData,
    coords: np.ndarray,
    if_bg: np.ndarray,
    overlay: ScoreOverlay,
    out_path: Path,
    title: str,
    highlight_celltype: str | None,
    color_map: dict[str, tuple[float, ...]],
    spot_size: float,
    highlight_opacity: float,
    if_alpha: float = 0.85,
    dpi: int = 300,
) -> None:
    """Save an alignment PNG with masked score overlay and optional highlight borders.

    Args:
        adata: AnnData with score and ``_cellestimate_plot`` columns.
        coords: Pixel coordinates shaped ``(n_spots, 2)``.
        if_bg: Aligned IF background array in [0, 1].
        overlay: Score overlay metadata including display limits.
        out_path: Output PNG path.
        title: Figure title.
        highlight_celltype: Cell type to outline, if any.
        color_map: Cell-type name to matplotlib RGB tuple.
        spot_size: Matplotlib scatter marker size for score bins.
        highlight_opacity: Alpha for highlight border rings in [0, 1].
        if_alpha: Alpha for the IF background image.
        dpi: PNG resolution.

    Returns:
        None

    Raises:
        ValueError: If ``highlight_opacity`` is outside [0, 1].
    """
    if not 0.0 <= highlight_opacity <= 1.0:
        raise ValueError(
            f"highlight_opacity must be in [0, 1], got {highlight_opacity}"
        )

    scores = adata.obs[overlay.score_col].to_numpy(dtype=float)
    celltypes = adata.obs["_cellestimate_plot"].astype(str)
    display_mask = score_display_mask(scores, overlay.vmin, overlay.vmax)

    fig, ax = plt.subplots(figsize=(10, 10))
    if if_bg.ndim == 2:
        ax.imshow(if_bg, cmap="gray", alpha=if_alpha)
    else:
        ax.imshow(if_bg, alpha=if_alpha)

    if highlight_celltype and highlight_celltype in color_map:
        highlight_mask = celltypes.to_numpy() == highlight_celltype
        edge_color = color_to_matplotlib_rgba(
            color_map[highlight_celltype],
            highlight_opacity,
        )
        ax.scatter(
            coords[highlight_mask, 0],
            coords[highlight_mask, 1],
            s=spot_size,
            facecolors="none",
            edgecolors=[edge_color],
            linewidths=0.8,
            zorder=3,
        )

    scatter = ax.scatter(
        coords[display_mask, 0],
        coords[display_mask, 1],
        c=scores[display_mask],
        cmap=DEFAULT_SCORE_COLORMAP,
        s=spot_size,
        alpha=0.9,
        linewidths=0,
        zorder=4,
        vmin=overlay.vmin,
        vmax=overlay.vmax,
    )
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(f"{overlay.score_label} [{overlay.vmin:g}, {overlay.vmax:g}]")

    ax.set_title(title, fontsize=12)
    ax.set_aspect("equal")
    ax.axis("off")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_interactive_figure(
    adata: sc.AnnData,
    coords: np.ndarray,
    if_bg_display: np.ndarray,
    extent_shape: tuple[int, int],
    overlay: ScoreOverlay,
    celltype_col: str,
    highlight_celltype: str | None,
    spot_size: float,
    score_opacity: float = 0.8,
    highlight_opacity: float = 0.45,
) -> "go.Figure":
    """Build a Plotly figure with IF, highlight borders, and score on top.

    Args:
        adata: AnnData with score and ``_cellestimate_plot`` columns.
        coords: Pixel coordinates (n_spots, 2) in full hires space.
        if_bg_display: Possibly downsampled IF raster for HTML embedding.
        extent_shape: Full hires canvas shape ``(height, width)`` for axis limits.
        overlay: Score overlay metadata.
        celltype_col: Resolved Loupe column name (for titles).
        highlight_celltype: Initial highlighted cell type, if any.
        spot_size: Base marker size for score layer.
        score_opacity: Opacity for score markers above highlight borders.
        highlight_opacity: Opacity for cellestimate highlight border rings.

    Returns:
        go.Figure: Interactive Plotly figure.

    Raises:
        ImportError: If plotly is not installed.
        ValueError: If opacity values are outside [0, 1].
    """
    if not 0.0 <= score_opacity <= 1.0:
        raise ValueError(f"score_opacity must be in [0, 1], got {score_opacity}")
    if not 0.0 <= highlight_opacity <= 1.0:
        raise ValueError(
            f"highlight_opacity must be in [0, 1], got {highlight_opacity}"
        )

    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise ImportError(
            "plotly is required for HTML output. Install with: pip install plotly"
        ) from exc

    extent_height, extent_width = extent_shape
    disp_height, disp_width = if_bg_display.shape[:2]
    scores = adata.obs[overlay.score_col].to_numpy(dtype=float)
    celltypes = adata.obs["_cellestimate_plot"].astype(str)
    unique_types = sorted({ct for ct in celltypes.unique() if ct and ct.lower() != "nan"})
    color_map = build_ubclass_color_map(unique_types)

    if highlight_celltype and highlight_celltype not in unique_types:
        warnings.warn(
            f"--highlight_celltype {highlight_celltype!r} not in Loupe labels; "
            f"available: {', '.join(unique_types)}",
            stacklevel=2,
        )
        highlight_celltype = unique_types[0] if unique_types else None

    hover = [
        f"Barcode: {bc}<br>"
        f"{overlay.score_label}: {val:.4g}<br>"
        f"{celltype_col}: {ct}<br>"
        f"x: {x:.1f}, y: {y:.1f}"
        for bc, val, ct, x, y in zip(
            adata.obs_names,
            scores,
            celltypes,
            coords[:, 0],
            coords[:, 1],
            strict=True,
        )
    ]

    initial_title = f"{overlay.score_label} — highlight {highlight_celltype} ({celltype_col})"
    if highlight_celltype is None:
        initial_title = f"{overlay.score_label} — no highlight"

    if_rgb = if_array_to_uint8_rgb(if_bg_display)
    score_mask = score_display_mask(scores, overlay.vmin, overlay.vmax)
    score_hover = [hover[i] for i in np.where(score_mask)[0]]
    fig = go.Figure()
    fig.add_trace(
        go.Image(
            z=if_rgb,
            **plotly_if_image_kwargs(
                extent_height,
                extent_width,
                disp_height,
                disp_width,
            ),
            hoverinfo="skip",
            name="IF",
        )
    )

    highlight_size = spot_size * 1.0
    for cell_type in unique_types:
        mask = celltypes.to_numpy() == cell_type
        border_color = rgba_to_plotly(color_map[cell_type], alpha=highlight_opacity)
        visible = cell_type == highlight_celltype
        fig.add_trace(
            go.Scatter(
                x=coords[mask, 0],
                y=coords[mask, 1],
                mode="markers",
                name=f"highlight {cell_type}",
                marker=dict(
                    size=highlight_size,
                    color="rgba(0,0,0,0)",
                    line=dict(color=border_color, width=2.0),
                ),
                text=[hover[i] for i in np.where(mask)[0]],
                hoverinfo="text",
                visible=visible,
                showlegend=True,
            )
        )

    fig.add_trace(
        go.Scatter(
            x=coords[score_mask, 0],
            y=coords[score_mask, 1],
            mode="markers",
            name=overlay.score_label,
            marker=dict(
                size=spot_size,
                color=scores[score_mask],
                colorscale=DEFAULT_SCORE_COLORMAP,
                cmin=overlay.vmin,
                cmax=overlay.vmax,
                opacity=score_opacity,
                line=dict(width=0),
                colorbar=dict(title=overlay.score_label),
            ),
            text=score_hover,
            hoverinfo="text",
        )
    )

    n_highlight = len(unique_types)
    visibility_none = [True] + [False] * n_highlight + [True]
    buttons: list[dict] = [
        dict(
            label="None",
            method="update",
            args=[
                {"visible": visibility_none},
                {"title": f"{overlay.score_label} — no highlight"},
            ],
        )
    ]
    for idx, cell_type in enumerate(unique_types):
        visibility = [True] + [False] * n_highlight + [True]
        visibility[1 + idx] = True
        buttons.append(
            dict(
                label=cell_type,
                method="update",
                args=[
                    {"visible": visibility},
                    {
                        "title": (
                            f"{overlay.score_label} — highlight {cell_type} "
                            f"({celltype_col})"
                        )
                    },
                ],
            )
        )

    fig.update_layout(
        title=initial_title,
        width=1000,
        height=1000,
        xaxis=dict(range=[0, extent_width], constrain="domain"),
        yaxis=dict(range=[extent_height, 0], scaleanchor="x", scaleratio=1),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=20, r=20, t=80, b=20),
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                showactive=True,
                x=0.01,
                xanchor="left",
                y=1.02,
                yanchor="bottom",
            )
        ],
    )

    return fig


def main() -> None:
    """Load data, score spots, and write interactive HTML plus spot table.

    Returns:
        None
    """
    args = parse_args()
    if not 0 < args.if_downsample <= 1.0:
        raise ValueError(
            f"--if_downsample must be in (0, 1], got {args.if_downsample}"
        )
    if not 0.0 <= args.score_opacity <= 1.0:
        raise ValueError(
            f"--score_opacity must be in [0, 1], got {args.score_opacity}"
        )
    if not 0.0 <= args.highlight_opacity <= 1.0:
        raise ValueError(
            f"--highlight_opacity must be in [0, 1], got {args.highlight_opacity}"
        )
    if args.score_vmin is not None and args.score_vmax is not None:
        if args.score_vmin >= args.score_vmax:
            raise ValueError(
                f"--score_vmin ({args.score_vmin}) must be less than "
                f"--score_vmax ({args.score_vmax})"
            )

    h5ad_path = resolve_h5ad_path(args)
    label = args.label or infer_label_from_h5ad(h5ad_path)
    loupe_csv = resolve_loupe_csv_path(args, label)
    output_dir = args.output_dir or (_REPO_ROOT / "results" / args.cohort)

    print(f"Loading {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    if "spatial" not in adata.obsm or "spatial" not in adata.uns:
        raise ValueError(f"{h5ad_path} missing spatial coordinates or images")

    print(f"Using Loupe CSV: {loupe_csv}")
    adata, celltype_col = merge_loupe_celltypes(
        adata,
        loupe_csv=loupe_csv,
        barcode_col=args.barcode_col,
        celltype_col=args.celltype_col,
    )
    n_types = adata.obs["_cellestimate_plot"].astype(str).nunique()
    print(f"Merged {len(adata)} spots with {n_types} {celltype_col} labels")

    if_path = resolve_if_image_path(args, label)
    print(f"Using IF image: {if_path}")

    library_id = list(adata.uns["spatial"].keys())[0]
    hires_shape = adata.uns["spatial"][library_id]["images"][args.img_key].shape
    if_bg_full = prepare_if_background(if_path, hires_shape, args.if_channel)
    if_bg_full = apply_if_orientation(
        if_bg_full,
        flip_x=args.if_flip_x,
        flip_y=args.if_flip_y,
    )

    html_if_path = resolve_if_html_path(args, label, if_path)
    if html_if_path != if_path:
        print(f"Using HTML IF image: {html_if_path}")
    if_bg_html = prepare_html_if_background(
        html_if_path,
        hires_shape,
        args.if_channel,
        args.if_downsample,
        flip_x=args.if_flip_x,
        flip_y=args.if_flip_y,
    )
    html_pixels = if_bg_html.shape[0] * if_bg_html.shape[1]
    full_pixels = if_bg_full.shape[0] * if_bg_full.shape[1]
    print(
        f"HTML IF raster: {if_bg_html.shape[1]}x{if_bg_html.shape[0]} "
        f"({html_pixels:,} px; {html_pixels / full_pixels:.1%} of full hires)"
    )

    overlay = compute_score_overlay(adata, args)
    print(
        f"Score overlay: {overlay.score_label} "
        f"[{overlay.vmin:g}, {overlay.vmax:g}]"
    )

    coords = spatial_pixel_coords(adata, args.img_key)
    celltypes = adata.obs["_cellestimate_plot"].astype(str)
    unique_types = sorted(
        {ct for ct in celltypes.unique() if ct and ct.lower() != "nan"}
    )
    color_map = build_ubclass_color_map(unique_types)
    highlight_celltype = args.highlight_celltype
    if highlight_celltype and highlight_celltype not in unique_types:
        warnings.warn(
            f"--highlight_celltype {highlight_celltype!r} not in Loupe labels; "
            f"available: {', '.join(unique_types)}",
            stacklevel=2,
        )
        highlight_celltype = unique_types[0] if unique_types else None

    fig = build_interactive_figure(
        adata,
        coords=coords,
        if_bg_display=if_bg_html,
        extent_shape=if_bg_full.shape[:2],
        overlay=overlay,
        celltype_col=celltype_col,
        highlight_celltype=highlight_celltype,
        spot_size=args.spot_size,
        score_opacity=args.score_opacity,
        highlight_opacity=args.highlight_opacity,
    )

    html_path = output_html_path(output_dir, label, overlay.score_slug)
    csv_path = output_csv_path(output_dir, label, overlay.score_slug)
    if_png_path = output_if_png_path(output_dir, label, overlay.score_slug)
    if_html_png_path = output_if_html_png_path(output_dir, label, overlay.score_slug)
    align_check_path = output_align_check_path(output_dir, label, overlay.score_slug)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    save_if_png(if_bg_full, if_png_path)
    save_if_png(if_bg_html, if_html_png_path)
    plot_score_celltype_align_check(
        adata,
        coords=coords,
        if_bg=if_bg_full,
        overlay=overlay,
        out_path=align_check_path,
        title=f"{overlay.score_label} alignment check",
        highlight_celltype=highlight_celltype,
        color_map=color_map,
        spot_size=args.spot_size,
        highlight_opacity=args.highlight_opacity,
    )
    fig.write_html(html_path, include_plotlyjs="cdn", full_html=True)
    write_spot_table(adata, coords, overlay.score_col, celltype_col, csv_path)

    print(f"Saved: {html_path}")
    print(f"Saved: {if_png_path}")
    print(f"Saved: {if_html_png_path}")
    print(f"Saved: {align_check_path}")
    print(f"Saved: {csv_path}")


if __name__ == "__main__":
    main()
