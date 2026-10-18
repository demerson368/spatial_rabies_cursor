#!/usr/bin/env python
"""Score Visium HD spots for Allen/Yao signature gene sets or individual genes.

Loads an extended spatial h5ad from ``data/processed/{cohort}/`` and writes
spatial PNGs with spot values overlaid on a registered IF image.

Default mode: scanpy gene-set scores for GABA/Glut intersection, GABA-unique,
Glut-unique, and zr_fusion reference lists.

Exclusive single-gene mode: pass ``--gene`` to plot per-spot expression for
one or more genes instead of signature sets.
"""

from __future__ import annotations

import argparse
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image
from skimage import transform

COHORT = "Tuyu_Apr2026"
DEFAULT_IF_DIR = (
    "/research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/"
    "common/Tuyu_Apr2026/IF_images"
)
MIN_GENES = 2 #3 

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

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "scripts"))

from celltype_biomarker_enrichment import (  # noqa: E402
    build_gene_map,
    get_count_matrix,
    get_expr_matrix,
    load_gene_set,
)


@dataclass(frozen=True)
class SignatureScoreResult:
    """Outcome of scoring one reference gene set."""

    set_key: str
    score_col: str
    reference_path: Path
    reference_genes: list[str]
    matched_genes: list[str]
    missing_genes: list[str]
    scored: bool


@dataclass(frozen=True)
class GeneExpressionResult:
    """Outcome of extracting per-spot expression for one gene."""

    requested_gene: str
    matched_var: str | None
    obs_col: str
    expression_layer: str
    min_value: float | None
    max_value: float | None
    plotted: bool


def parse_gene_args(gene_args: list[str] | None) -> list[str]:
    """Flatten and deduplicate gene symbols from repeated ``--gene`` flags.

    Args:
        gene_args: Values collected from ``parser.add_argument("--gene", action="append")``.

    Returns:
        list[str]: Deduplicated gene symbols in request order.

    Raises:
        ValueError: If ``gene_args`` is non-empty but yields no symbols after parsing.
    """
    if not gene_args:
        return []

    genes: list[str] = []
    seen: set[str] = set()
    for token in gene_args:
        for part in re.split(r"[,\s]+", token):
            gene = part.strip()
            if gene and gene not in seen:
                seen.add(gene)
                genes.append(gene)

    if gene_args and not genes:
        raise ValueError("No gene symbols found in --gene arguments")
    return genes


def gene_obs_slug(gene_symbol: str) -> str:
    """Build a filesystem-safe slug from a gene symbol.

    Args:
        gene_symbol: Gene symbol as stored in ``adata.var_names``.

    Returns:
        str: Slug with non-alphanumeric runs replaced by ``_``.
    """
    slug = re.sub(r"[^A-Za-z0-9]+", "_", gene_symbol).strip("_")
    return slug or "gene"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for GABA/Glut signature scoring.

    Returns:
        argparse.Namespace: Parsed paths, plot options, and output flags.
    """
    repo_root = _REPO_ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=None,
        help="Extended h5ad under data/processed/{cohort}/ (spatial_{label}_ext.h5ad)",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Sample label for output filenames and default h5ad lookup",
    )
    parser.add_argument(
        "--cohort",
        default=COHORT,
        help="Cohort subdirectory under data/processed/",
    )
    parser.add_argument(
        "--if_image",
        type=Path,
        default=None,
        help=(
            "IF TIFF/PNG path. May be a full path or a filename resolved under "
            "--if_dir. Omit to auto-pick from --if_dir using config/samples.tsv."
        ),
    )
    parser.add_argument(
        "--if_dir",
        type=Path,
        default=Path(DEFAULT_IF_DIR),
        help=(
            "Directory of cohort IF images. Used for auto-discovery when --if_image "
            "is omitted, and as fallback when --if_image is a filename only."
        ),
    )
    parser.add_argument(
        "--if_channel",
        choices=("red", "green", "gray", "rgb"),
        default="green",
        help="Channel to display from RGB IF stacks (default green for GFP)",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=None,
        help="Directory for PNG and coverage CSV (default results/{cohort})",
    )
    parser.add_argument(
        "--gaba_glut_intersection_genes",
        type=Path,
        default=repo_root / "reference" / SIGNATURE_SETS["gaba_glut_intersection"][0],
        help="GABA ∩ Glut shared gene list",
    )
    parser.add_argument(
        "--gaba_unique_genes",
        type=Path,
        default=repo_root / "reference" / SIGNATURE_SETS["gaba_unique"][0],
        help="GABA-exclusive gene list",
    )
    parser.add_argument(
        "--glut_unique_genes",
        type=Path,
        default=repo_root / "reference" / SIGNATURE_SETS["glut_unique"][0],
        help="Glut-exclusive gene list",
    )
    parser.add_argument(
        "--zr_fusion_genes",
        type=Path,
        default=repo_root / "reference" / SIGNATURE_SETS["zr_fusion"][0],
        help="zr_fusion gene list",
    )
    parser.add_argument(
        "--img_key",
        default="hires",
        help="Spatial image key in adata.uns['spatial'] for alignment dimensions",
    )
    parser.add_argument(
        "--spot_size",
        type=float,
        default=4.0,
        help="Matplotlib scatter marker size for Visium HD bins",
    )
    parser.add_argument(
        "--if_alpha",
        type=float,
        default=0.85,
        help="Alpha for IF background image",
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
        "--dpi",
        type=int,
        default=300,
        help="PNG resolution",
    )
    parser.add_argument(
        "--write_h5ad",
        type=Path,
        default=None,
        help="Optional output h5ad path to persist score columns",
    )
    parser.add_argument(
        "--use_preprocess_h5ad",
        action="store_true",
        help="Resolve spatial_{label}.h5ad (preprocess) instead of _ext.h5ad",
    )
    parser.add_argument(
        "--gene",
        action="append",
        default=None,
        help=(
            "Gene symbol to plot (repeatable; comma-separated ok). "
            "Exclusive: skips signature sets when set."
        ),
    )
    parser.add_argument(
        "--expression_layer",
        choices=("normalized", "counts"),
        default="normalized",
        help="Expression layer for --gene mode (default: log-normalized adata.X)",
    )
    parser.epilog = (
        "Examples:\n"
        "  # Signature sets (default):\n"
        "  python score_gaba_glut_signatures.py \\\n"
        "    --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um\n\n"
        "  # Single gene (exclusive — no signature sets):\n"
        "  python score_gaba_glut_signatures.py --label ... --gene Gad1\n\n"
        "  # Multiple genes:\n"
        "  python score_gaba_glut_signatures.py --label ... --gene Gad1 --gene Slc17a7\n\n"
        "  # Gene counts with fixed colormap:\n"
        "  python score_gaba_glut_signatures.py --label ... --gene Gad1 \\\n"
        "    --expression_layer counts --score_vmin 0 --score_vmax 5\n\n"
        "  # Explicit IF directory (omit --if_image; picks file via samples.tsv):\n"
        "  python score_gaba_glut_signatures.py --label ... \\\n"
        "    --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/"
        "common/Tuyu_Apr2026/IF_images\n\n"
        "  # IF filename only (resolved under --if_dir, not repo root):\n"
        "  python score_gaba_glut_signatures.py --label ... --if_dir .../IF_images \\\n"
        "    --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif\n\n"
        "  # Fixed colormap range across plots:\n"
        "  python score_gaba_glut_signatures.py --label ... \\\n"
        "    --score_vmin -0.2 --score_vmax 0.8"
    )
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    return parser.parse_args()


def infer_label_from_h5ad(h5ad_path: Path) -> str:
    """Infer sample label from an extended h5ad filename.

    Args:
        h5ad_path: Path like ``spatial_{label}_ext.h5ad`` or ref-tagged variant.

    Returns:
        str: Sample label stem ending in bin size (e.g. ``..._008um``).
    """
    stem = h5ad_path.stem
    if stem.startswith("spatial_"):
        stem = stem[len("spatial_") :]
    if stem.endswith("_ext"):
        stem = stem[: -len("_ext")]

    match = re.match(r"^(Tuyu_Apr2026_.+?_\d{3}um)", stem)
    if match:
        return match.group(1)
    return stem


def resolve_h5ad_path(args: argparse.Namespace) -> Path:
    """Resolve the extended h5ad input path from CLI arguments.

    Args:
        args: Parsed CLI namespace with ``h5ad``, ``label``, and ``cohort``.

    Returns:
        Path: Existing extended h5ad file.

    Raises:
        FileNotFoundError: If no matching h5ad can be resolved.
        ValueError: If neither ``--h5ad`` nor ``--label`` is provided.
    """
    if args.h5ad is not None:
        h5ad = args.h5ad.expanduser().resolve()
        if not h5ad.is_file():
            raise FileNotFoundError(f"h5ad not found: {h5ad}")
        return h5ad

    if not args.label:
        raise ValueError("Provide --h5ad or --label to locate the extended h5ad")

    processed_dir = _REPO_ROOT / "data" / "processed" / args.cohort
    if args.use_preprocess_h5ad:
        direct = processed_dir / f"spatial_{args.label}.h5ad"
        if direct.is_file():
            return direct
        raise FileNotFoundError(
            f"No preprocess h5ad found for label {args.label!r} under {processed_dir}"
        )

    direct = processed_dir / f"spatial_{args.label}_ext.h5ad"
    if direct.is_file():
        return direct

    matches = sorted(processed_dir.glob(f"spatial_{args.label}_*_ext.h5ad"))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise FileNotFoundError(
            f"Multiple ref-tagged extended h5ad files match label {args.label!r}: "
            f"{', '.join(p.name for p in matches)}. Pass --h5ad explicitly."
        )
    
    print(f"spatial_{args.label}_*_ext.h5ad")
    raise FileNotFoundError(
        f"No extended h5ad found for label {args.label!r} under {processed_dir}"
    )


def sr_sample_dir_from_label(label: str) -> str | None:
    """Look up Space Ranger sample dir token for a pipeline label.

    Args:
        label: Sample label matching ``config/samples.tsv``.

    Returns:
        str | None: ``sr_sample_dir`` value (e.g. ``846531_HD_..._D90-5-A1-GFP_R1``)
            or ``None`` if not found.
    """
    samples_tsv = _REPO_ROOT / "config" / "samples.tsv"
    if not samples_tsv.is_file():
        return None
    df = pd.read_csv(samples_tsv, sep="\t", comment="#")
    df.columns = df.columns.str.strip()
    matches = df.loc[df["label"].astype(str) == label, "sr_sample_dir"]
    if len(matches) == 1:
        return str(matches.iloc[0])
    return None


def sample_id_from_label(label: str) -> str:
    """Extract the Space Ranger sample id token from a pipeline label.

    Args:
        label: Sample label like ``Tuyu_Apr2026_846531_HD_BA_WM_..._008um``.

    Returns:
        str: Sample id portion before the bin suffix (e.g. ``846531_HD_..._GFP_R1``).
    """
    parts = label.split("_")
    if len(parts) >= 3 and parts[0] == "Tuyu" and parts[1] == "Apr2026":
        core = parts[2:]
    else:
        core = parts
    if core and core[-1].endswith("um"):
        core = core[:-1]
    return "_".join(core)


def resolve_if_image_path(args: argparse.Namespace, label: str) -> Path:
    """Resolve the IF image path from CLI arguments.

    Resolution order when ``--if_image`` is set:

    1. Path as given (absolute or relative to the current working directory)
    2. ``{if_dir}/{if_image}`` when step 1 is missing
    3. ``{if_dir}/{basename}`` when ``if_image`` is a filename only

    When ``--if_image`` is omitted, auto-discover under ``--if_dir`` using
    ``sr_sample_dir`` from ``config/samples.tsv``.

    Args:
        args: Parsed CLI namespace with ``if_image`` and ``if_dir``.
        label: Sample label used for auto-discovery in ``if_dir``.

    Returns:
        Path: Existing IF image file.

    Raises:
        FileNotFoundError: If no IF image can be resolved.
    """
    if args.if_image is not None:
        raw = args.if_image.expanduser()
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
            f"IF image not found: {args.if_image}. Tried: {tried}. "
            f"Use a full path or pass only the filename with --if_dir set."
        )

    if_dir = args.if_dir.expanduser().resolve()
    if not if_dir.is_dir():
        raise FileNotFoundError(
            f"--if_image not provided and IF directory missing: {if_dir}"
        )

    sample_id = sr_sample_dir_from_label(label) or sample_id_from_label(label)
    patterns = (f"*{sample_id}*",)
    candidates: list[Path] = []
    suffixes = (".tif", ".tiff", ".png", ".jpg", ".jpeg")
    for pattern in patterns:
        for match in if_dir.glob(pattern):
            if match.suffix.lower() in suffixes:
                candidates.append(match)
    unique = sorted({p.resolve() for p in candidates if p.is_file()})
    if len(unique) == 1:
        return unique[0]
    if len(unique) > 1:
        raise FileNotFoundError(
            f"Multiple IF images match sample {label!r} in {if_dir}: "
            f"{', '.join(p.name for p in unique)}. Pass --if_image explicitly."
        )

    raise FileNotFoundError(
        f"No IF image found for sample {label!r} under {if_dir}. "
        f"Pass --if_image or check config/samples.tsv sr_sample_dir."
    )


def load_if_image_array(if_path: Path) -> np.ndarray:
    """Load an IF image as a float array in ``[0, 1]``.

    Args:
        if_path: Path to a TIFF or PNG IF image.

    Returns:
        np.ndarray: Image array shaped ``(H, W)`` or ``(H, W, C)``.
    """
    with Image.open(if_path) as img:
        arr = np.asarray(img)
    if arr.dtype == np.uint8:
        arr = arr.astype(np.float32) / 255.0
    elif arr.dtype == np.uint16:
        arr = arr.astype(np.float32) / 65535.0
    else:
        arr = arr.astype(np.float32)
        if arr.max() > 1.0:
            arr = arr / arr.max()
    return arr


def prepare_if_background(
    if_path: Path,
    target_shape: tuple[int, ...],
    channel: str,
) -> np.ndarray:
    """Load and resize an IF image to match the spatial hires canvas.

    Args:
        if_path: Path to the IF image file.
        target_shape: Shape of ``adata.uns['spatial'][library_id]['images']['hires']``.
        channel: One of ``red``, ``green``, ``gray``, or ``rgb``.

    Returns:
        np.ndarray: Background array for ``imshow`` — 2D grayscale or 3-channel RGB.
    """
    if_arr = load_if_image_array(if_path)
    target_hw = (target_shape[0], target_shape[1])
    resized = transform.resize(if_arr, target_hw, preserve_range=True, anti_aliasing=True)
    resized = np.clip(resized, 0.0, 1.0)

    if resized.ndim == 2:
        if channel == "rgb":
            return np.stack([resized, resized, resized], axis=-1)
        return resized

    if channel == "rgb":
        if resized.shape[-1] >= 3:
            return resized[..., :3]
        return np.stack([resized[..., 0]] * 3, axis=-1)

    if channel == "gray":
        if resized.shape[-1] >= 3:
            return resized[..., :3].mean(axis=-1)
        return resized[..., 0]

    channel_idx = {"red": 0, "green": 1}[channel]
    if resized.shape[-1] <= channel_idx:
        return resized[..., 0]
    return resized[..., channel_idx]


def score_signature_set(
    adata: sc.AnnData,
    set_key: str,
    score_col: str,
    gene_path: Path,
    min_genes: int = MIN_GENES,
) -> SignatureScoreResult:
    """Score spots for one Allen/Yao signature gene list.

    Args:
        adata: Spatial AnnData with log-normalized expression in ``.X``.
        set_key: Short key for outputs (e.g. ``gaba_unique``).
        score_col: ``adata.obs`` column name for the score.
        gene_path: Path to the reference gene list file.
        min_genes: Minimum matched genes required to run ``score_genes``.

    Returns:
        SignatureScoreResult: Coverage metadata and whether scoring succeeded.
    """
    reference_genes = load_gene_set(gene_path)
    gene_map = build_gene_map(reference_genes, adata.var_names.tolist())
    matched = [gene_map[g] for g in reference_genes if g in gene_map]
    missing = [g for g in reference_genes if g not in gene_map]

    if len(matched) < min_genes:
        warnings.warn(
            f"{set_key}: only {len(matched)}/{len(reference_genes)} genes matched "
            f"(need >= {min_genes}); skipping score.",
            stacklevel=2,
        )
        return SignatureScoreResult(
            set_key=set_key,
            score_col=score_col,
            reference_path=gene_path,
            reference_genes=reference_genes,
            matched_genes=matched,
            missing_genes=missing,
            scored=False,
        )

    sc.tl.score_genes(
        adata,
        gene_list=matched,
        score_name=score_col,
        use_raw=False,
    )
    return SignatureScoreResult(
        set_key=set_key,
        score_col=score_col,
        reference_path=gene_path,
        reference_genes=reference_genes,
        matched_genes=matched,
        missing_genes=missing,
        scored=True,
    )


def score_single_gene(
    adata: sc.AnnData,
    gene_symbol: str,
    expression_layer: str,
) -> GeneExpressionResult:
    """Extract per-spot expression for one gene into ``adata.obs``.

    Args:
        adata: Spatial AnnData with expression in ``.X`` or ``layers['counts']``.
        gene_symbol: Requested gene symbol (case-insensitive match).
        expression_layer: ``normalized`` for ``adata.X`` or ``counts`` for raw counts.

    Returns:
        GeneExpressionResult: Metadata including ``obs_col`` and whether the gene matched.
    """
    gene_map = build_gene_map([gene_symbol], adata.var_names.tolist())
    if gene_symbol not in gene_map:
        warnings.warn(
            f"Gene {gene_symbol!r} not found in adata.var_names; skipping.",
            stacklevel=2,
        )
        return GeneExpressionResult(
            requested_gene=gene_symbol,
            matched_var=None,
            obs_col="",
            expression_layer=expression_layer,
            min_value=None,
            max_value=None,
            plotted=False,
        )

    matched_var = gene_map[gene_symbol]
    obs_col = f"expr_gene_{gene_obs_slug(matched_var)}"

    if expression_layer == "counts":
        matrix = get_count_matrix(adata)
    else:
        matrix = get_expr_matrix(adata)

    var_index = {g: i for i, g in enumerate(adata.var_names)}
    values = np.asarray(matrix[:, var_index[matched_var]]).ravel()
    adata.obs[obs_col] = values

    return GeneExpressionResult(
        requested_gene=gene_symbol,
        matched_var=matched_var,
        obs_col=obs_col,
        expression_layer=expression_layer,
        min_value=float(np.nanmin(values)),
        max_value=float(np.nanmax(values)),
        plotted=True,
    )


def score_display_mask(
    scores: np.ndarray,
    vmin: float | None,
    vmax: float | None,
) -> np.ndarray:
    """Return a boolean mask for spots within the display score range.

    Args:
        scores: Per-spot score values.
        vmin: Minimum display threshold; ``None`` uses the data minimum.
        vmax: Maximum display threshold; ``None`` uses the data maximum.

    Returns:
        np.ndarray: Boolean mask shaped like ``scores``; ``False`` for out-of-range
            and NaN values.
    """
    lo = vmin if vmin is not None else float(np.nanmin(scores))
    hi = vmax if vmax is not None else float(np.nanmax(scores))
    return (scores >= lo) & (scores <= hi)


def plot_score_on_if(
    adata: sc.AnnData,
    score_col: str,
    if_bg: np.ndarray,
    out_path: Path,
    title: str,
    img_key: str = "hires",
    spot_size: float = 4.0,
    if_alpha: float = 0.85,
    score_vmin: float | None = None,
    score_vmax: float | None = None,
    dpi: int = 300,
) -> None:
    """Save a spatial plot of spot scores overlaid on an IF background.

    Args:
        adata: AnnData with ``obsm['spatial']`` and spatial metadata.
        score_col: ``adata.obs`` column with continuous scores.
        if_bg: IF background array aligned to hires dimensions.
        out_path: Output PNG path.
        title: Figure title.
        img_key: Image resolution key in ``adata.uns['spatial']``.
        spot_size: Matplotlib scatter marker size.
        if_alpha: Alpha for the IF background.
        score_vmin: Minimum display threshold; spots below are hidden. ``None`` uses data min.
        score_vmax: Maximum display threshold; spots above are hidden. ``None`` uses data max.
        dpi: PNG resolution.

    Returns:
        None

    Raises:
        ValueError: If both limits are set and ``score_vmin >= score_vmax``.
    """
    if score_vmin is not None and score_vmax is not None and score_vmin >= score_vmax:
        raise ValueError(
            f"score_vmin ({score_vmin}) must be less than score_vmax ({score_vmax})"
        )

    library_id = list(adata.uns["spatial"].keys())[0]
    scale = adata.uns["spatial"][library_id]["scalefactors"][f"tissue_{img_key}_scalef"]
    coords = adata.obsm["spatial"] * scale
    scores = adata.obs[score_col].to_numpy()
    display_mask = score_display_mask(scores, score_vmin, score_vmax)

    fig, ax = plt.subplots(figsize=(10, 10))
    if if_bg.ndim == 2:
        ax.imshow(if_bg, cmap="gray", alpha=if_alpha)
    else:
        ax.imshow(if_bg, alpha=if_alpha)

    scatter = ax.scatter(
        coords[display_mask, 0],
        coords[display_mask, 1],
        c=scores[display_mask],
        cmap="magma",
        s=spot_size,
        alpha=0.9,
        linewidths=0,
        zorder=3,
        vmin=score_vmin,
        vmax=score_vmax,
    )
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    if score_vmin is not None or score_vmax is not None:
        lo = score_vmin if score_vmin is not None else float(np.nanmin(scores))
        hi = score_vmax if score_vmax is not None else float(np.nanmax(scores))
        cbar.set_label(f"{score_col} [{lo:g}, {hi:g}]")
    else:
        cbar.set_label(score_col)

    ax.set_title(title, fontsize=12)
    ax.set_aspect("equal")
    ax.axis("off")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_gene_coverage_csv(results: list[SignatureScoreResult], out_path: Path) -> None:
    """Write per-set gene coverage summary to CSV.

    Args:
        results: Scoring results for each signature set.
        out_path: Output CSV path.

    Returns:
        None
    """
    rows: list[dict[str, object]] = []
    for result in results:
        rows.append(
            {
                "set_key": result.set_key,
                "score_col": result.score_col,
                "reference_path": str(result.reference_path),
                "n_reference": len(result.reference_genes),
                "n_matched": len(result.matched_genes),
                "n_missing": len(result.missing_genes),
                "scored": result.scored,
                "matched_genes": ";".join(result.matched_genes),
                "missing_genes": ";".join(result.missing_genes),
            }
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, index=False)


def png_output_path(output_dir: Path, label: str, set_key: str) -> Path:
    """Build the PNG output path for one signature set.

    Args:
        output_dir: Results directory.
        label: Sample label.
        set_key: Signature set key (e.g. ``gaba_unique``).

    Returns:
        Path: ``spatial_{label}_{set_key}_allen_score.png``.
    """
    return output_dir / f"spatial_{label}_{set_key}_allen_score.png"


def gene_png_output_path(output_dir: Path, label: str, gene_symbol: str) -> Path:
    """Build the PNG output path for one gene expression plot.

    Args:
        output_dir: Results directory.
        label: Sample label.
        gene_symbol: Matched gene symbol in ``adata.var_names``.

    Returns:
        Path: ``spatial_{label}_gene_{gene}_allen_score.png``.
    """
    return output_dir / f"spatial_{label}_gene_{gene_obs_slug(gene_symbol)}_allen_score.png"


def write_gene_expression_summary_csv(
    results: list[GeneExpressionResult],
    out_path: Path,
) -> None:
    """Write per-gene expression summary for single-gene mode.

    Args:
        results: Gene expression extraction results.
        out_path: Output CSV path.

    Returns:
        None
    """
    rows: list[dict[str, object]] = []
    for result in results:
        rows.append(
            {
                "requested_gene": result.requested_gene,
                "matched_var": result.matched_var or "",
                "obs_col": result.obs_col,
                "expression_layer": result.expression_layer,
                "min_value": result.min_value,
                "max_value": result.max_value,
                "plotted": result.plotted,
            }
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, index=False)


def run_signature_scoring(
    adata: sc.AnnData,
    args: argparse.Namespace,
    label: str,
    output_dir: Path,
    if_bg: np.ndarray,
    gene_paths: dict[str, Path],
) -> list[SignatureScoreResult]:
    """Score and plot all signature gene sets.

    Args:
        adata: Spatial AnnData object.
        args: Parsed CLI namespace.
        label: Sample label for output filenames.
        output_dir: Directory for PNG and CSV outputs.
        if_bg: IF background image array.
        gene_paths: Mapping from set key to reference gene list path.

    Returns:
        list[SignatureScoreResult]: Results for each signature set.
    """
    results: list[SignatureScoreResult] = []
    for set_key, (_, score_col) in SIGNATURE_SETS.items():
        result = score_signature_set(
            adata,
            set_key=set_key,
            score_col=score_col,
            gene_path=gene_paths[set_key],
        )
        results.append(result)
        print(
            f"{set_key}: {len(result.matched_genes)}/{len(result.reference_genes)} "
            f"genes matched"
            + (" — scored" if result.scored else " — skipped")
        )

        if not result.scored:
            continue

        title = (
            f"{set_key.replace('_', ' ')} "
            f"({len(result.matched_genes)}/{len(result.reference_genes)} genes)"
        )
        out_png = png_output_path(output_dir, label, set_key)
        plot_score_on_if(
            adata,
            score_col=score_col,
            if_bg=if_bg,
            out_path=out_png,
            title=title,
            img_key=args.img_key,
            spot_size=args.spot_size,
            if_alpha=args.if_alpha,
            score_vmin=args.score_vmin,
            score_vmax=args.score_vmax,
            dpi=args.dpi,
        )
        print(f"Saved: {out_png}")

    coverage_csv = output_dir / f"spatial_{label}_gaba_glut_fusion_signature_gene_coverage.csv"
    write_gene_coverage_csv(results, coverage_csv)
    print(f"Saved: {coverage_csv}")
    return results


def run_gene_scoring(
    adata: sc.AnnData,
    args: argparse.Namespace,
    label: str,
    output_dir: Path,
    if_bg: np.ndarray,
    genes: list[str],
) -> list[GeneExpressionResult]:
    """Plot per-spot expression for requested genes (exclusive mode).

    Args:
        adata: Spatial AnnData object.
        args: Parsed CLI namespace.
        label: Sample label for output filenames.
        output_dir: Directory for PNG and CSV outputs.
        if_bg: IF background image array.
        genes: Requested gene symbols.

    Returns:
        list[GeneExpressionResult]: Results for each requested gene.

    Raises:
        ValueError: If no requested genes match ``adata.var_names``.
    """
    results: list[GeneExpressionResult] = []
    for gene in genes:
        result = score_single_gene(adata, gene, args.expression_layer)
        results.append(result)
        if not result.plotted or result.matched_var is None:
            continue

        print(
            f"{gene} -> {result.matched_var} "
            f"({result.expression_layer}: {result.min_value:g}–{result.max_value:g})"
        )
        out_png = gene_png_output_path(output_dir, label, result.matched_var)
        plot_score_on_if(
            adata,
            score_col=result.obs_col,
            if_bg=if_bg,
            out_path=out_png,
            title=f"{result.matched_var} ({args.expression_layer})",
            img_key=args.img_key,
            spot_size=args.spot_size,
            if_alpha=args.if_alpha,
            score_vmin=args.score_vmin,
            score_vmax=args.score_vmax,
            dpi=args.dpi,
        )
        print(f"Saved: {out_png}")

    if not any(r.plotted for r in results):
        unmatched = [r.requested_gene for r in results]
        raise ValueError(
            f"No requested genes matched adata.var_names: {', '.join(unmatched)}"
        )

    summary_csv = output_dir / f"spatial_{label}_gene_expression_summary.csv"
    write_gene_expression_summary_csv(results, summary_csv)
    print(f"Saved: {summary_csv}")
    return results


def main() -> None:
    """Load h5ad, score signature sets or single genes, and write IF-overlay PNGs.

    Returns:
        None
    """
    args = parse_args()
    genes = parse_gene_args(args.gene)
    if args.score_vmin is not None and args.score_vmax is not None:
        if args.score_vmin >= args.score_vmax:
            raise ValueError(
                f"--score_vmin ({args.score_vmin}) must be less than "
                f"--score_vmax ({args.score_vmax})"
            )

    h5ad_path = resolve_h5ad_path(args)
    label = args.label or infer_label_from_h5ad(h5ad_path)
    output_dir = args.output_dir or (_REPO_ROOT / "results" / args.cohort)

    print(f"Loading {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    if "spatial" not in adata.obsm or "spatial" not in adata.uns:
        raise ValueError(f"{h5ad_path} missing spatial coordinates or images")

    if_path = resolve_if_image_path(args, label)
    print(f"Using IF image: {if_path}")

    library_id = list(adata.uns["spatial"].keys())[0]
    hires_shape = adata.uns["spatial"][library_id]["images"][args.img_key].shape
    if_bg = prepare_if_background(if_path, hires_shape, args.if_channel)

    if genes:
        run_gene_scoring(adata, args, label, output_dir, if_bg, genes)
    else:
        gene_paths = {
            "gaba_glut_intersection": args.gaba_glut_intersection_genes,
            "gaba_unique": args.gaba_unique_genes,
            "glut_unique": args.glut_unique_genes,
            "zr_fusion": args.zr_fusion_genes,
        }
        run_signature_scoring(adata, args, label, output_dir, if_bg, gene_paths)

    if args.write_h5ad is not None:
        args.write_h5ad.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(args.write_h5ad)
        print(f"Saved: {args.write_h5ad}")


if __name__ == "__main__":
    main()
