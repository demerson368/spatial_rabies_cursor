#!/usr/bin/env python
"""Interactive HTML viewer v2: optional gene/signature score with multi-column cell-type highlights.

Labels come from **``adata.obs``** (repeat ``--celltype_col`` up to twice,
e.g. ``rabies_assignment`` and ``cellestimate``). Dropdown includes singles and
all cross-column pairs. ``--no_score`` skips the score layer.

Rabies-style obs columns (name contains ``rabies``): booleans and 0/1 map to
``positive`` / ``negative``; other columns use string labels as-is.

Requires: pip install plotly

Example (h5ad obs lanes, no score)::

    python plot_score_celltype_html_v2.py \\
      --h5ad data/processed/.../spatial_*_rabies_assignment.h5ad \\
      --if_dir ... --if_image ... \\
      --no_score \\
      --celltype_col rabies_assignment --celltype_col cellestimate \\
      --highlight_celltype positive --highlight_celltype Glut \\
      --highlight_opacity 1.0 --highlight_opacity 0.2

Example (signature plus fusion score from obs)::

    python plot_score_celltype_html_v2.py \\
      --h5ad ... --if_dir ... --if_image ... \\
      --signature nrxn_signature --score_vmin 2 \\
      --extra_score_col tltumscore_C11orf95_relafus1 \\
      --extra_score_vmin 0.1 --extra_score_vmax 0.5
"""

from __future__ import annotations

import argparse
import itertools
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
    "nrxn_signature": (
        "nrxn.txt",
        "nrxn_signature",
    ),
    "nlgn_signature": (
        "nlgn.txt",
        "nlgn_signature",
    ),
    "signatures_set1": (
        "signatures_set1.txt",
        "signatures_set1",
    ),
    "signatures_set2": (
        "signatures_set2.txt",
        "signatures_set2",
    ),
    "signatures_synaptic": (
        "signatures_synaptic.txt",
        "signatures_synaptic",
    ),

}
COHORT = "Tuyu_Apr2026"
DEFAULT_IF_DIR = (
    "/research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/"
    "common/Tuyu_Apr2026/IF_images"
)
DEFAULT_CELLTYPE_COL = "cellestimate"
DEFAULT_SCORE_COLORMAP = "Oranges_r" #"Purples"
DEFAULT_EXTRA_SCORE_COLORMAP = "Greens" #"Reds"
DEFAULT_EXTRA_SCORE_OPACITY = 0.45
CELLTYPE_COL_ALIASES = {"cell_estimate": "cellestimate"}
MAX_MENU_ITEMS_WARN = 200
NO_SCORE_SLUG = "celltype_only"
DEFAULT_HIGHLIGHT_OPACITY = 0.45

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
EXTRA_LABEL_COLS: dict[str, str] = {
    "positive": "#00FFFF",
    "negative": "#888888",
}


def hex_color_for_label(label: str) -> str:
    """Resolve highlight hex from UBCLASS, extra (e.g. rabies), or gray fallback."""
    if label in UBCLASS_COLS:
        return UBCLASS_COLS[label]
    lower = label.lower()
    if lower in EXTRA_LABEL_COLS:
        return EXTRA_LABEL_COLS[lower]
    return UBCLASS_FALLBACK_COLOR


def build_ubclass_color_map(
    cell_types: list[str] | set[str],
) -> dict[str, tuple[float, float, float]]:
    """Assign colors to highlight labels for plotting."""
    color_map: dict[str, tuple[float, float, float]] = {}
    for label in sorted(set(cell_types)):
        hex_color = hex_color_for_label(label)
        if label not in UBCLASS_COLS and label.lower() not in EXTRA_LABEL_COLS:
            warnings.warn(
                f"No color entry for {label!r}; using gray",
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
    score_log: str = "none"
    explicit_score_limits: bool = False
    colorscale: str = DEFAULT_SCORE_COLORMAP
    marker_opacity: float = 0.8


@dataclass(frozen=True)
class CelltypeLane:
    """One ``adata.obs`` annotation column used for highlight borders."""

    obs_key: str
    source_col: str
    display_name: str


@dataclass(frozen=True)
class LaneLabel:
    """Cell-type label on a specific annotation lane."""

    lane_idx: int
    cell_type: str


@dataclass(frozen=True)
class HighlightMode:
    """One dropdown option controlling which highlight traces are visible."""

    menu_label: str
    title: str
    active: tuple[LaneLabel, ...]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the v2 interactive HTML viewer."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--h5ad", type=Path, default=None, help="Extended h5ad path")
    parser.add_argument("--label", default=None, help="Sample label for defaults")
    parser.add_argument("--cohort", default=COHORT, help="Cohort under data/processed/")
    parser.add_argument(
        "--celltype_col",
        action="append",
        default=None,
        help="adata.obs column(s) for highlights; repeat for second lane (max 2)",
    )
    parser.add_argument(
        "--highlight_celltype",
        action="append",
        default=None,
        help="Initial highlight class(es); pair with second --celltype_col by order",
    )
    parser.add_argument(
        "--highlight_opacity",
        action="append",
        type=float,
        default=None,
        help=(
            "Highlight border opacity (0-1); pairs with --highlight_celltype "
            "by order; one value without --highlight_celltype applies to all traces"
        ),
    )
    parser.add_argument(
        "--max_pair_options",
        type=int,
        default=80,
        help="Max pair entries in dropdown before error",
    )
    parser.add_argument(
        "--no_score",
        action="store_true",
        help="Skip gene/signature scoring; IF and highlights only",
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
        help="IF channel to display",
    )
    parser.add_argument("--if_flip_y", action="store_true", help="Flip IF vertically")
    parser.add_argument("--if_flip_x", action="store_true", help="Flip IF horizontally")
    parser.add_argument(
        "--if_image_html",
        type=Path,
        default=None,
        help="Optional lower-res IF for HTML only",
    )
    parser.add_argument(
        "--if_downsample",
        type=float,
        default=0.25,
        help="Downsample factor for HTML IF layer",
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
        help="Signature set to score (requires --gene omitted)",
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
        help="Opacity of score markers (0-1)",
    )
    parser.add_argument("--score_vmin", type=float, default=None, help="Score display min")
    parser.add_argument("--score_vmax", type=float, default=None, help="Score display max")
    parser.add_argument(
        "--score_log",
        choices=("none", "log1p", "log10"),
        default="none",
        help="Log transform applied to scores for display, masking, and color scale",
    )
    parser.add_argument(
        "--extra_score_col",
        default=None,
        help=(
            "Additional score from adata.obs (e.g. tltumscore_C11orf95_relafus1); "
            "norm_ prefix added if needed"
        ),
    )
    parser.add_argument(
        "--extra_score_label",
        default=None,
        help="Display label for --extra_score_col (default: column name)",
    )
    parser.add_argument(
        "--extra_score_vmin",
        type=float,
        default=None,
        help="Display min for --extra_score_col",
    )
    parser.add_argument(
        "--extra_score_vmax",
        type=float,
        default=None,
        help="Display max for --extra_score_col",
    )
    parser.add_argument(
        "--extra_score_log",
        choices=("none", "log1p", "log10"),
        default="none",
        help="Log transform for --extra_score_col",
    )
    parser.add_argument(
        "--extra_score_colorscale",
        default=DEFAULT_EXTRA_SCORE_COLORMAP,
        help="Plotly colorscale for --extra_score_col",
    )
    parser.add_argument(
        "--extra_score_opacity",
        type=float,
        default=None,
        help="Marker opacity for extra score (default: 0.45 when stacked)",
    )
    parser.add_argument(
        "--use_preprocess_h5ad",
        action="store_true",
        help="Resolve spatial_{label}.h5ad instead of _ext.h5ad",
    )
    return parser.parse_args()


def requested_celltype_columns(args: argparse.Namespace) -> list[str]:
    """Normalize ``--celltype_col`` to a list of 1–2 column names."""
    cols = args.celltype_col
    if cols is None or len(cols) == 0:
        return [DEFAULT_CELLTYPE_COL]
    if len(cols) > 2:
        raise ValueError(
            f"At most 2 --celltype_col values allowed, got {len(cols)}: {cols}"
        )
    return cols


def resolve_highlight_opacity_by_celltype(
    highlight_values: list[str],
    raw_opacities: list[float] | None,
) -> dict[str, float]:
    """Map cell-type labels to highlight border opacity from CLI repeats.

    Args:
        highlight_values: ``--highlight_celltype`` values in order.
        raw_opacities: ``--highlight_opacity`` values in order, or None.

    Returns:
        dict[str, float]: Label to opacity; empty if defaults apply everywhere.
        Special key ``"__global__"`` means one opacity for all traces when no
        ``--highlight_celltype`` was given.

    Raises:
        ValueError: If opacities are out of range or exceed cell-type count.
    """
    if raw_opacities is None or len(raw_opacities) == 0:
        return {}

    for alpha in raw_opacities:
        if not 0.0 <= alpha <= 1.0:
            raise ValueError(
                f"--highlight_opacity must be in [0, 1], got {alpha}"
            )

    if not highlight_values:
        if len(raw_opacities) > 1:
            raise ValueError(
                "Multiple --highlight_opacity values require matching "
                "--highlight_celltype entries"
            )
        return {"__global__": raw_opacities[0]}

    if len(raw_opacities) > len(highlight_values):
        raise ValueError(
            f"More --highlight_opacity ({len(raw_opacities)}) than "
            f"--highlight_celltype ({len(highlight_values)}) values"
        )

    if len(raw_opacities) < len(highlight_values):
        warnings.warn(
            f"Fewer --highlight_opacity ({len(raw_opacities)}) than "
            f"--highlight_celltype ({len(highlight_values)}); unpaired labels "
            f"use default {DEFAULT_HIGHLIGHT_OPACITY}",
            stacklevel=2,
        )

    result: dict[str, float] = {}
    for label, alpha in zip(highlight_values, raw_opacities, strict=False):
        result[label] = alpha
    return result


def highlight_trace_opacity(
    cell_type: str,
    opacity_by_celltype: dict[str, float],
) -> float:
    """Return border opacity for one highlight cell-type label."""
    if "__global__" in opacity_by_celltype:
        return opacity_by_celltype["__global__"]
    return opacity_by_celltype.get(cell_type, DEFAULT_HIGHLIGHT_OPACITY)


def resolve_obs_column(adata: sc.AnnData, requested: str) -> str:
    """Resolve an ``adata.obs`` column name, honoring aliases.

    Args:
        adata: Spatial AnnData object.
        requested: CLI ``--celltype_col`` value.

    Returns:
        str: Resolved column name present in ``adata.obs``.

    Raises:
        ValueError: If no matching column exists.
    """
    canonical = CELLTYPE_COL_ALIASES.get(requested.lower(), requested)
    columns = adata.obs.columns
    if canonical in columns:
        return str(canonical)
    if requested in columns:
        return str(requested)
    lower_map = {str(c).lower(): str(c) for c in columns}
    if canonical.lower() in lower_map:
        return lower_map[canonical.lower()]
    if requested.lower() in lower_map:
        return lower_map[requested.lower()]
    available = ", ".join(str(c) for c in columns)
    raise ValueError(
        f"obs column {requested!r} not found in h5ad. Available: {available}"
    )


def resolve_obs_score_column(adata: sc.AnnData, requested: str) -> str:
    """Resolve a continuous score column in ``adata.obs``.

    Tries the requested name, ``norm_{name}``, and case-insensitive matches.

    Args:
        adata: Spatial AnnData object.
        requested: CLI ``--extra_score_col`` value.

    Returns:
        str: Resolved column name.

    Raises:
        ValueError: If no matching column exists.
    """
    candidates: list[str] = [requested]
    if not requested.startswith("norm_"):
        candidates.append(f"norm_{requested}")
    else:
        candidates.append(requested[len("norm_") :])
    lower_map = {str(c).lower(): str(c) for c in adata.obs.columns}
    for candidate in candidates:
        if candidate in adata.obs.columns:
            return str(candidate)
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    available = ", ".join(str(c) for c in adata.obs.columns)
    raise ValueError(
        f"Score obs column {requested!r} not found in h5ad. Available: {available}"
    )


def obs_score_slug(column: str) -> str:
    """Filesystem slug from an obs score column name."""
    slug = re.sub(r"^norm_", "", column)
    slug = re.sub(r"[^\w]+", "_", slug).strip("_")
    return slug[:80] if slug else "obs_score"


def overlay_score_summary(overlays: list[ScoreOverlay]) -> str:
    """Short label for titles when one or more score layers are active."""
    if not overlays:
        return "cell types only"
    if len(overlays) == 1:
        return overlays[0].score_label
    return " + ".join(overlay.score_label for overlay in overlays)


def _coerce_rabies_label(value: object) -> str:
    """Map rabies-assignment values to ``positive`` or ``negative`` strings."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    if isinstance(value, (bool, np.bool_)):
        return "positive" if value else "negative"
    if isinstance(value, (int, np.integer)):
        return "positive" if int(value) == 1 else "negative" if int(value) == 0 else str(value)
    if isinstance(value, (float, np.floating)) and not np.isnan(value):
        if value in (0.0, 0):
            return "negative"
        if value in (1.0, 1):
            return "positive"
    text = str(value).strip()
    if text.lower() in ("nan", "none", ""):
        return ""
    lowered = text.lower()
    if lowered in ("true", "1", "yes", "positive"):
        return "positive"
    if lowered in ("false", "0", "no", "negative"):
        return "negative"
    return text


def normalize_lane_labels(series: pd.Series, col_name: str) -> pd.Series:
    """Convert obs values to stable highlight label strings.

    Args:
        series: Raw ``adata.obs`` column.
        col_name: Resolved obs column name.

    Returns:
        pd.Series: Normalized string labels (empty where missing).
    """
    if "rabies" in col_name.lower():
        return series.map(_coerce_rabies_label)
    labels = series.astype(str).replace({"nan": "", "None": "", "NaN": ""})
    return labels


def attach_obs_label_lanes(
    adata: sc.AnnData,
    obs_cols: list[str],
) -> list[CelltypeLane]:
    """Attach 1–2 normalized label lanes from ``adata.obs``.

    Args:
        adata: Spatial AnnData (modified in place for ``_celltype_lane*`` keys).
        obs_cols: Requested obs column names (1–2).

    Returns:
        list[CelltypeLane]: Lane metadata for plotting.
    """
    lanes: list[CelltypeLane] = []
    for idx, requested in enumerate(obs_cols):
        resolved = resolve_obs_column(adata, requested)
        lane = CelltypeLane(
            obs_key=f"_celltype_lane{idx}",
            source_col=resolved,
            display_name=resolved,
        )
        adata.obs[lane.obs_key] = normalize_lane_labels(
            adata.obs[resolved], resolved
        ).astype(str)
        lanes.append(lane)
    return lanes


def label_from_h5ad_path(h5ad_path: Path) -> str:
    """Infer a short sample label from an h5ad filename for outputs and IF lookup.

    Args:
        h5ad_path: Path to the spatial h5ad file.

    Returns:
        str: Sample label (prefer ``..._008um``-style token when present).
    """
    stem = h5ad_path.stem
    if stem.startswith("spatial_"):
        stem = stem[len("spatial_") :]
    match = re.match(r"^(Tuyu_Apr2026_.+?_\d{3}um)", stem)
    if match:
        return match.group(1)
    for suffix in (
        "_merged_nontumor_v1_broad_ext_rabies_assignment",
        "_rabies_assignment",
        "_ext",
    ):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return infer_label_from_h5ad(h5ad_path)


def lane_unique_labels(adata: sc.AnnData, lane: CelltypeLane) -> list[str]:
    """Sorted non-empty unique labels for one lane."""
    series = adata.obs[lane.obs_key].astype(str)
    return sorted(
        {ct for ct in series.unique() if ct and ct.lower() not in ("nan", "none", "")}
    )


def signature_gene_path(set_key: str) -> Path:
    """Resolve the reference gene list path for a signature set key."""
    filename = SIGNATURE_SETS[set_key][0]
    return _REPO_ROOT / "reference" / filename


def transform_scores_for_display(
    scores: np.ndarray,
    score_log: str,
) -> np.ndarray:
    """Apply optional log transform for score display and masking.

    Args:
        scores: Raw per-spot score values.
        score_log: One of ``none``, ``log1p``, ``log10``.

    Returns:
        np.ndarray: Transformed scores (same shape as input).

    Raises:
        ValueError: If ``score_log`` is not supported.
    """
    arr = np.asarray(scores, dtype=float)
    if score_log == "none":
        return arr
    if score_log == "log1p":
        return np.log1p(np.clip(arr, 0.0, None))
    if score_log == "log10":
        return np.log10(np.clip(arr, 1e-10, None))
    raise ValueError(f"Unsupported score_log {score_log!r}")


def resolve_score_color_limits(
    values: np.ndarray,
    score_vmin: float | None,
    score_vmax: float | None,
) -> tuple[float, float, bool]:
    """Resolve display limits on the (possibly log-transformed) score array.

    Args:
        values: Display-scale score values.
        score_vmin: CLI minimum or ``None`` for data min.
        score_vmax: CLI maximum or ``None`` for data max.

    Returns:
        tuple[float, float, bool]: ``vmin``, ``vmax``, and whether either limit
            was set explicitly on the CLI.
    """
    explicit = score_vmin is not None or score_vmax is not None
    vmin = score_vmin if score_vmin is not None else float(np.nanmin(values))
    vmax = score_vmax if score_vmax is not None else float(np.nanmax(values))
    return vmin, vmax, explicit


def resolve_score_overlay_v2(
    adata: sc.AnnData,
    args: argparse.Namespace,
) -> ScoreOverlay | None:
    """Score a gene or signature set, or return None when ``--no_score``."""
    if args.no_score:
        return None

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
        raw = adata.obs[result.obs_col].to_numpy(dtype=float)
        values = transform_scores_for_display(raw, args.score_log)
        vmin, vmax, explicit = resolve_score_color_limits(
            values, args.score_vmin, args.score_vmax
        )
        return ScoreOverlay(
            score_col=result.obs_col,
            score_label=f"{result.matched_var} ({args.expression_layer})",
            score_slug=f"gene_{gene_obs_slug(result.matched_var)}",
            vmin=vmin,
            vmax=vmax,
            score_log=args.score_log,
            explicit_score_limits=explicit,
            colorscale=DEFAULT_SCORE_COLORMAP,
            marker_opacity=args.score_opacity,
        )

    if args.signature:
        set_key = args.signature
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
        values = transform_scores_for_display(
            adata.obs[score_col].to_numpy(dtype=float),
            args.score_log,
        )
        vmin, vmax, explicit = resolve_score_color_limits(
            values, args.score_vmin, args.score_vmax
        )
        return ScoreOverlay(
            score_col=score_col,
            score_label=set_key.replace("_", " "),
            score_slug=set_key,
            vmin=vmin,
            vmax=vmax,
            score_log=args.score_log,
            explicit_score_limits=explicit,
            colorscale=DEFAULT_SCORE_COLORMAP,
            marker_opacity=args.score_opacity,
        )

    raise ValueError(
        "Specify --gene, --signature, --extra_score_col, or --no_score"
    )


def resolve_extra_score_overlay(
    adata: sc.AnnData,
    args: argparse.Namespace,
) -> ScoreOverlay | None:
    """Build a score overlay from an existing ``adata.obs`` column.

    Args:
        adata: Spatial AnnData object.
        args: Parsed CLI namespace.

    Returns:
        ScoreOverlay | None: Overlay metadata, or ``None`` if not requested.
    """
    if not args.extra_score_col:
        return None

    col = resolve_obs_score_column(adata, args.extra_score_col)
    raw = adata.obs[col].to_numpy(dtype=float)
    values = transform_scores_for_display(raw, args.extra_score_log)
    vmin, vmax, explicit = resolve_score_color_limits(
        values,
        args.extra_score_vmin,
        args.extra_score_vmax,
    )
    label = args.extra_score_label or col
    extra_opacity = args.extra_score_opacity
    if extra_opacity is None:
        extra_opacity = DEFAULT_EXTRA_SCORE_OPACITY
    return ScoreOverlay(
        score_col=col,
        score_label=label,
        score_slug=f"extra_{obs_score_slug(col)}",
        vmin=vmin,
        vmax=vmax,
        score_log=args.extra_score_log,
        explicit_score_limits=explicit,
        colorscale=args.extra_score_colorscale,
        marker_opacity=extra_opacity,
    )


def build_score_overlays(
    adata: sc.AnnData,
    args: argparse.Namespace,
) -> list[ScoreOverlay]:
    """Resolve primary (gene/signature) and optional extra obs score overlays.

    Args:
        adata: Spatial AnnData object.
        args: Parsed CLI namespace.

    Returns:
        list[ScoreOverlay]: Zero, one, or two overlays in draw order.

    Raises:
        ValueError: If no score source is requested.
    """
    overlays: list[ScoreOverlay] = []
    genes = parse_gene_args(args.gene)
    wants_primary = not args.no_score and (bool(genes) or args.signature is not None)
    if wants_primary:
        overlays.append(resolve_score_overlay_v2(adata, args))
    extra = resolve_extra_score_overlay(adata, args)
    if extra is not None:
        overlays.append(extra)
    if not overlays:
        if args.no_score and not args.extra_score_col:
            return []
        raise ValueError(
            "Specify --gene, --signature, --extra_score_col, or --no_score"
        )
    return overlays


def score_slug_for_outputs(overlays: list[ScoreOverlay]) -> str:
    """Filesystem slug for output filenames."""
    if not overlays:
        return NO_SCORE_SLUG
    if len(overlays) == 1:
        return overlays[0].score_slug
    return "_".join(overlay.score_slug for overlay in overlays)


def apply_if_orientation(if_bg: np.ndarray, flip_x: bool, flip_y: bool) -> np.ndarray:
    """Apply optional horizontal/vertical flips to an IF background array."""
    oriented = if_bg
    if flip_y:
        oriented = np.flipud(oriented)
    if flip_x:
        oriented = np.fliplr(oriented)
    return oriented


def downsample_if_array(if_bg: np.ndarray, factor: float) -> np.ndarray:
    """Downsample an IF background array for lighter HTML embedding."""
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
    """Load, align, orient, and downsample an IF array for HTML display."""
    if_bg = prepare_if_background(if_path, hires_shape, channel)
    if_bg = apply_if_orientation(if_bg, flip_x=flip_x, flip_y=flip_y)
    return downsample_if_array(if_bg, downsample)


def resolve_if_html_path(
    args: argparse.Namespace, label: str, default_if_path: Path
) -> Path:
    """Resolve the IF image path used for the HTML background layer."""
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
    """Return Plotly ``go.Image`` kwargs scaled to a full coordinate extent."""
    return dict(
        x0=0,
        y0=0,
        dx=full_width / disp_width,
        dy=full_height / disp_height,
    )


def spatial_pixel_coords(adata: sc.AnnData, img_key: str) -> np.ndarray:
    """Return spot coordinates scaled to hires image pixel space."""
    library_id = list(adata.uns["spatial"].keys())[0]
    scale = adata.uns["spatial"][library_id]["scalefactors"][f"tissue_{img_key}_scalef"]
    return np.asarray(adata.obsm["spatial"]) * scale


def rgba_to_plotly(color: tuple[float, ...] | str, alpha: float = 1.0) -> str:
    """Convert a matplotlib RGB or RGBA tuple to a Plotly rgba() string."""
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
    """Convert a matplotlib RGB or RGBA tuple to an RGBA tuple."""
    if len(color) == 3:
        r, g, b = color
    elif len(color) == 4:
        r, g, b, _ = color
    else:
        raise ValueError(f"Expected RGB or RGBA tuple, got length {len(color)}")
    return (float(r), float(g), float(b), float(alpha))


def if_array_to_uint8_rgb(if_bg: np.ndarray) -> np.ndarray:
    """Convert an IF background array to uint8 RGB for image display."""
    if if_bg.ndim == 2:
        gray = (np.clip(if_bg, 0.0, 1.0) * 255).astype(np.uint8)
        return np.stack([gray, gray, gray], axis=-1)
    rgb = (np.clip(if_bg[..., :3], 0.0, 1.0) * 255).astype(np.uint8)
    if rgb.shape[-1] == 1:
        return np.concatenate([rgb, rgb, rgb], axis=-1)
    return rgb


def save_if_png(if_bg: np.ndarray, out_path: Path) -> None:
    """Save the aligned IF background as a PNG file."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(if_array_to_uint8_rgb(if_bg)).save(out_path)


def output_if_png_path(output_dir: Path, label: str, score_slug: str) -> Path:
    return output_dir / f"spatial_{label}_{score_slug}_if_background.png"


def output_if_html_png_path(output_dir: Path, label: str, score_slug: str) -> Path:
    return output_dir / f"spatial_{label}_{score_slug}_if_background_html.png"


def output_align_check_path(output_dir: Path, label: str, score_slug: str) -> Path:
    return output_dir / f"spatial_{label}_{score_slug}_celltype_v2_align_check.png"


def output_html_path(output_dir: Path, label: str, score_slug: str) -> Path:
    return output_dir / f"spatial_{label}_{score_slug}_celltype_v2.html"


def output_csv_path(output_dir: Path, label: str, score_slug: str) -> Path:
    return output_dir / f"spatial_{label}_{score_slug}_celltype_v2_spots.csv"


def format_lane_label(lane: CelltypeLane, cell_type: str) -> str:
    """Short label for menus: ``col: type``."""
    return f"{lane.display_name}: {cell_type}"


def build_highlight_modes(
    lanes: list[CelltypeLane],
    labels_by_lane: list[list[str]],
    max_pair_options: int,
    overlays: list[ScoreOverlay],
) -> list[HighlightMode]:
    """Build None, single, and all-pair highlight dropdown modes."""
    score_part = overlay_score_summary(overlays)
    modes: list[HighlightMode] = [
        HighlightMode(
            menu_label="None",
            title=f"{score_part} — no highlight",
            active=(),
        )
    ]

    for lane_idx, lane in enumerate(lanes):
        for cell_type in labels_by_lane[lane_idx]:
            modes.append(
                HighlightMode(
                    menu_label=format_lane_label(lane, cell_type),
                    title=f"{score_part} — highlight {format_lane_label(lane, cell_type)}",
                    active=(LaneLabel(lane_idx, cell_type),),
                )
            )

    pair_modes: list[HighlightMode] = []
    if len(lanes) == 1:
        types0 = labels_by_lane[0]
        for t1, t2 in itertools.combinations(types0, 2):
            pair_modes.append(
                HighlightMode(
                    menu_label=f"{t1} + {t2}",
                    title=f"{score_part} — highlight {t1} + {t2} ({lanes[0].display_name})",
                    active=(
                        LaneLabel(0, t1),
                        LaneLabel(0, t2),
                    ),
                )
            )
    else:
        types0, types1 = labels_by_lane[0], labels_by_lane[1]
        for t0, t1 in itertools.product(types0, types1):
            pair_modes.append(
                HighlightMode(
                    menu_label=(
                        f"{format_lane_label(lanes[0], t0)} + "
                        f"{format_lane_label(lanes[1], t1)}"
                    ),
                    title=(
                        f"{score_part} — "
                        f"{format_lane_label(lanes[0], t0)} + "
                        f"{format_lane_label(lanes[1], t1)}"
                    ),
                    active=(
                        LaneLabel(0, t0),
                        LaneLabel(1, t1),
                    ),
                )
            )

    if len(pair_modes) > max_pair_options:
        raise ValueError(
            f"Pair dropdown would have {len(pair_modes)} entries (max "
            f"{max_pair_options}). Use fewer types, one column, or increase "
            f"--max_pair_options."
        )

    modes.extend(pair_modes)

    if len(modes) > MAX_MENU_ITEMS_WARN:
        warnings.warn(
            f"Highlight dropdown has {len(modes)} entries; HTML may be slow to use",
            stacklevel=2,
        )
    return modes


def resolve_initial_highlight_mode(
    modes: list[HighlightMode],
    lanes: list[CelltypeLane],
    labels_by_lane: list[list[str]],
    highlight_values: list[str] | None,
) -> HighlightMode:
    """Pick initial dropdown mode from CLI ``--highlight_celltype`` values."""
    if not modes:
        raise ValueError("No highlight modes defined")

    if not highlight_values:
        return modes[0]

    if len(highlight_values) == 1:
        target = highlight_values[0]
        for lane_idx, types in enumerate(labels_by_lane):
            if target in types:
                for mode in modes:
                    if mode.active == (LaneLabel(lane_idx, target),):
                        return mode
        warnings.warn(
            f"--highlight_celltype {target!r} not found in any lane; using None",
            stacklevel=2,
        )
        return modes[0]

    if len(highlight_values) >= 2:
        h0, h1 = highlight_values[0], highlight_values[1]
        if len(lanes) == 1:
            desired = frozenset({LaneLabel(0, h0), LaneLabel(0, h1)})
            for mode in modes:
                if frozenset(mode.active) == desired:
                    return mode
            warnings.warn(
                f"Pair {h0!r} + {h1!r} not in dropdown; using None",
                stacklevel=2,
            )
            return modes[0]

        desired = (LaneLabel(0, h0), LaneLabel(1, h1))
        for mode in modes:
            if mode.active == desired:
                return mode
        warnings.warn(
            f"Cross-column pair {h0!r} + {h1!r} not in dropdown; using None",
            stacklevel=2,
        )
        return modes[0]

    return modes[0]


def build_visibility_for_mode(
    mode: HighlightMode,
    n_highlight_traces: int,
    n_score_traces: int,
    trace_key_to_index: dict[tuple[int, str], int],
) -> list[bool]:
    """Visibility list: [IF, highlight traces..., score trace(s)]."""
    visibility = [True] + [False] * n_highlight_traces
    for lane_label in mode.active:
        key = (lane_label.lane_idx, lane_label.cell_type)
        if key in trace_key_to_index:
            visibility[1 + trace_key_to_index[key]] = True
    visibility.extend([True] * n_score_traces)
    return visibility


def build_spot_hover_lines(
    adata: sc.AnnData,
    lanes: list[CelltypeLane],
    coords: np.ndarray,
    overlays: list[ScoreOverlay],
) -> list[str]:
    """Per-spot hover HTML for Plotly."""
    lines: list[str] = []
    n = adata.n_obs
    lane_series = [adata.obs[lane.obs_key].astype(str) for lane in lanes]
    display_scores: list[np.ndarray] = []
    for overlay in overlays:
        raw = adata.obs[overlay.score_col].to_numpy(dtype=float)
        display_scores.append(transform_scores_for_display(raw, overlay.score_log))

    for i in range(n):
        parts = [f"Barcode: {adata.obs_names[i]}"]
        for overlay, scores in zip(overlays, display_scores, strict=True):
            parts.append(f"{overlay.score_label}: {scores[i]:.4g}")
        for lane, series in zip(lanes, lane_series, strict=True):
            parts.append(f"{lane.display_name}: {series.iloc[i]}")
        parts.append(f"x: {coords[i, 0]:.1f}, y: {coords[i, 1]:.1f}")
        lines.append("<br>".join(parts))
    return lines


def write_spot_table_v2(
    adata: sc.AnnData,
    coords: np.ndarray,
    lanes: list[CelltypeLane],
    overlays: list[ScoreOverlay],
    out_path: Path,
) -> None:
    """Write per-spot table with all lane columns and optional score."""
    data: dict[str, object] = {
        "Barcode": adata.obs_names,
        "x": coords[:, 0],
        "y": coords[:, 1],
    }
    for idx, overlay in enumerate(overlays):
        key = overlay.score_col if len(overlays) == 1 else f"score_{idx}_{overlay.score_col}"
        data[key] = adata.obs[overlay.score_col].to_numpy()
    for lane in lanes:
        data[lane.source_col] = adata.obs[lane.obs_key].astype(str)
    df = pd.DataFrame(data)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)


def plot_score_celltype_align_check_v2(
    adata: sc.AnnData,
    coords: np.ndarray,
    if_bg: np.ndarray,
    overlays: list[ScoreOverlay],
    out_path: Path,
    title: str,
    lanes: list[CelltypeLane],
    initial_mode: HighlightMode,
    color_map: dict[str, tuple[float, ...]],
    spot_size: float,
    highlight_opacity_by_celltype: dict[str, float] | None = None,
    if_alpha: float = 0.85,
    dpi: int = 300,
) -> None:
    """Save alignment PNG with optional score and highlight borders for initial mode."""
    opacity_map = highlight_opacity_by_celltype or {}

    fig, ax = plt.subplots(figsize=(10, 10))
    if if_bg.ndim == 2:
        ax.imshow(if_bg, cmap="gray", alpha=if_alpha)
    else:
        ax.imshow(if_bg, alpha=if_alpha)

    zorder = 3
    for lane_label in initial_mode.active:
        lane = lanes[lane_label.lane_idx]
        celltypes = adata.obs[lane.obs_key].astype(str)
        if lane_label.cell_type not in color_map:
            continue
        highlight_mask = celltypes.to_numpy() == lane_label.cell_type
        alpha = highlight_trace_opacity(lane_label.cell_type, opacity_map)
        edge_color = color_to_matplotlib_rgba(
            color_map[lane_label.cell_type],
            alpha,
        )
        ax.scatter(
            coords[highlight_mask, 0],
            coords[highlight_mask, 1],
            s=spot_size,
            facecolors="none",
            edgecolors=[edge_color],
            linewidths=0.8,
            zorder=zorder,
        )
        zorder += 1

    for overlay in overlays:
        raw_scores = adata.obs[overlay.score_col].to_numpy(dtype=float)
        scores = transform_scores_for_display(raw_scores, overlay.score_log)
        display_mask = score_display_mask(scores, overlay.vmin, overlay.vmax)
        scatter = ax.scatter(
            coords[display_mask, 0],
            coords[display_mask, 1],
            c=scores[display_mask],
            cmap=overlay.colorscale,
            s=spot_size,
            alpha=overlay.marker_opacity,
            linewidths=0,
            zorder=zorder + 1,
            vmin=overlay.vmin,
            vmax=overlay.vmax,
        )
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(f"{overlay.score_label} [{overlay.vmin:g}, {overlay.vmax:g}]")
        zorder += 1

    ax.set_title(title, fontsize=12)
    ax.set_aspect("equal")
    ax.axis("off")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_interactive_figure_v2(
    adata: sc.AnnData,
    coords: np.ndarray,
    if_bg_display: np.ndarray,
    extent_shape: tuple[int, int],
    lanes: list[CelltypeLane],
    labels_by_lane: list[list[str]],
    overlays: list[ScoreOverlay],
    modes: list[HighlightMode],
    initial_mode: HighlightMode,
    spot_size: float,
    highlight_opacity_by_celltype: dict[str, float] | None = None,
) -> "go.Figure":
    """Build Plotly figure with IF, per-lane highlights, optional score, pair dropdown."""
    opacity_map = highlight_opacity_by_celltype or {}

    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise ImportError(
            "plotly is required for HTML output. Install with: pip install plotly"
        ) from exc

    all_labels = [t for types in labels_by_lane for t in types]
    color_map = build_ubclass_color_map(all_labels)

    extent_height, extent_width = extent_shape
    disp_height, disp_width = if_bg_display.shape[:2]
    hover = build_spot_hover_lines(adata, lanes, coords, overlays)

    if_rgb = if_array_to_uint8_rgb(if_bg_display)
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
    trace_key_to_index: dict[tuple[int, str], int] = {}
    trace_idx = 0
    initial_active = frozenset(initial_mode.active)

    for lane_idx, lane in enumerate(lanes):
        celltypes = adata.obs[lane.obs_key].astype(str)
        for cell_type in labels_by_lane[lane_idx]:
            mask = celltypes.to_numpy() == cell_type
            alpha = highlight_trace_opacity(cell_type, opacity_map)
            border_color = rgba_to_plotly(color_map[cell_type], alpha=alpha)
            visible = LaneLabel(lane_idx, cell_type) in initial_active
            trace_key_to_index[(lane_idx, cell_type)] = trace_idx
            fig.add_trace(
                go.Scatter(
                    x=coords[mask, 0],
                    y=coords[mask, 1],
                    mode="markers",
                    name=f"highlight {format_lane_label(lane, cell_type)}",
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
            trace_idx += 1

    n_highlight = trace_idx
    n_score_traces = len(overlays)
    for score_idx, overlay in enumerate(overlays):
        raw_scores = adata.obs[overlay.score_col].to_numpy(dtype=float)
        scores = transform_scores_for_display(raw_scores, overlay.score_log)
        score_mask = score_display_mask(scores, overlay.vmin, overlay.vmax)
        score_hover = [hover[i] for i in np.where(score_mask)[0]]
        colorbar_y = 1.0 - score_idx * 0.38
        marker_kwargs: dict[str, object] = dict(
            size=spot_size,
            color=scores[score_mask],
            colorscale=overlay.colorscale,
            cmin=overlay.vmin,
            cmax=overlay.vmax,
            opacity=overlay.marker_opacity,
            line=dict(width=0),
            colorbar=dict(
                title=overlay.score_label,
                len=0.35,
                yanchor="top",
                y=colorbar_y,
                x=1.02,
            ),
        )
        if overlay.explicit_score_limits:
            marker_kwargs["cauto"] = False
        fig.add_trace(
            go.Scatter(
                x=coords[score_mask, 0],
                y=coords[score_mask, 1],
                mode="markers",
                name=overlay.score_label,
                marker=marker_kwargs,
                text=score_hover,
                hoverinfo="text",
            )
        )

    buttons: list[dict] = []
    for mode in modes:
        visibility = build_visibility_for_mode(
            mode,
            n_highlight_traces=n_highlight,
            n_score_traces=n_score_traces,
            trace_key_to_index=trace_key_to_index,
        )
        buttons.append(
            dict(
                label=mode.menu_label,
                method="update",
                args=[{"visible": visibility}, {"title": mode.title}],
            )
        )

    fig.update_layout(
        title=initial_mode.title,
        width=1000,
        height=1000,
        xaxis=dict(range=[0, extent_width], constrain="domain"),
        yaxis=dict(range=[extent_height, 0], scaleanchor="x", scaleratio=1),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=20, r=120 if n_score_traces > 1 else 20, t=80, b=20),
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
    """Load data, optional score, and write v2 HTML plus outputs."""
    args = parse_args()
    if not 0 < args.if_downsample <= 1.0:
        raise ValueError(
            f"--if_downsample must be in (0, 1], got {args.if_downsample}"
        )
    if not 0.0 <= args.score_opacity <= 1.0:
        raise ValueError(
            f"--score_opacity must be in [0, 1], got {args.score_opacity}"
        )
    if args.score_vmin is not None and args.score_vmax is not None:
        if args.score_vmin >= args.score_vmax:
            raise ValueError(
                f"--score_vmin ({args.score_vmin}) must be less than "
                f"--score_vmax ({args.score_vmax})"
            )
    if args.extra_score_vmin is not None and args.extra_score_vmax is not None:
        if args.extra_score_vmin >= args.extra_score_vmax:
            raise ValueError(
                f"--extra_score_vmin ({args.extra_score_vmin}) must be less than "
                f"--extra_score_vmax ({args.extra_score_vmax})"
            )
    if args.extra_score_opacity is not None and not 0.0 <= args.extra_score_opacity <= 1.0:
        raise ValueError(
            f"--extra_score_opacity must be in [0, 1], got {args.extra_score_opacity}"
        )

    celltype_cols = requested_celltype_columns(args)
    highlight_values = args.highlight_celltype or []
    highlight_opacity_by_celltype = resolve_highlight_opacity_by_celltype(
        highlight_values,
        args.highlight_opacity,
    )

    h5ad_path = resolve_h5ad_path(args)
    label = args.label or label_from_h5ad_path(h5ad_path)
    output_dir = args.output_dir or (_REPO_ROOT / "results" / args.cohort)

    print(f"Loading {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    if "spatial" not in adata.obsm or "spatial" not in adata.uns:
        raise ValueError(f"{h5ad_path} missing spatial coordinates or images")

    lanes = attach_obs_label_lanes(adata, celltype_cols)
    labels_by_lane = [lane_unique_labels(adata, lane) for lane in lanes]
    lane_summary = ", ".join(
        f"{lane.source_col} ({len(types)} types)"
        for lane, types in zip(lanes, labels_by_lane, strict=True)
    )
    print(f"Label lanes from h5ad obs — {len(adata)} spots: {lane_summary}")

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

    score_overlays = build_score_overlays(adata, args)
    if not score_overlays:
        print("Score overlay: disabled (--no_score, no --extra_score_col)")
    else:
        for overlay in score_overlays:
            log_note = f", log={overlay.score_log}" if overlay.score_log != "none" else ""
            print(
                f"Score overlay: {overlay.score_label} "
                f"[{overlay.vmin:g}, {overlay.vmax:g}]{log_note} "
                f"(col={overlay.score_col})"
            )

    score_slug = score_slug_for_outputs(score_overlays)
    modes = build_highlight_modes(
        lanes,
        labels_by_lane,
        max_pair_options=args.max_pair_options,
        overlays=score_overlays,
    )
    initial_mode = resolve_initial_highlight_mode(
        modes,
        lanes,
        labels_by_lane,
        highlight_values,
    )
    print(
        f"Highlight dropdown: {len(modes)} options "
        f"(initial: {initial_mode.menu_label})"
    )

    coords = spatial_pixel_coords(adata, args.img_key)
    all_labels = [t for types in labels_by_lane for t in types]
    color_map = build_ubclass_color_map(all_labels)

    fig = build_interactive_figure_v2(
        adata,
        coords=coords,
        if_bg_display=if_bg_html,
        extent_shape=if_bg_full.shape[:2],
        lanes=lanes,
        labels_by_lane=labels_by_lane,
        overlays=score_overlays,
        modes=modes,
        initial_mode=initial_mode,
        spot_size=args.spot_size,
        highlight_opacity_by_celltype=highlight_opacity_by_celltype,
    )

    html_path = output_html_path(output_dir, label, score_slug)
    csv_path = output_csv_path(output_dir, label, score_slug)
    if_png_path = output_if_png_path(output_dir, label, score_slug)
    if_html_png_path = output_if_html_png_path(output_dir, label, score_slug)
    align_check_path = output_align_check_path(output_dir, label, score_slug)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    save_if_png(if_bg_full, if_png_path)
    save_if_png(if_bg_html, if_html_png_path)

    align_title = initial_mode.title
    if score_overlays:
        align_title = (
            f"{overlay_score_summary(score_overlays)} alignment check — "
            f"{initial_mode.menu_label}"
        )
    else:
        align_title = f"Cell type alignment — {initial_mode.menu_label}"

    plot_score_celltype_align_check_v2(
        adata,
        coords=coords,
        if_bg=if_bg_full,
        overlays=score_overlays,
        out_path=align_check_path,
        title=align_title,
        lanes=lanes,
        initial_mode=initial_mode,
        color_map=color_map,
        spot_size=args.spot_size,
        highlight_opacity_by_celltype=highlight_opacity_by_celltype,
    )
    fig.write_html(html_path, include_plotlyjs="cdn", full_html=True)
    write_spot_table_v2(adata, coords, lanes, score_overlays, csv_path)

    print(f"Saved: {html_path}")
    print(f"Saved: {if_png_path}")
    print(f"Saved: {if_html_png_path}")
    print(f"Saved: {align_check_path}")
    print(f"Saved: {csv_path}")


if __name__ == "__main__":
    main()
