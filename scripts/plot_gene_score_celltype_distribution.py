#!/usr/bin/env python
"""Per-bin gene or signature score distributions by cell-type column.

Uses the same scoring as ``plot_score_celltype_html_v2.py``: ``--gene`` (``score_single_gene``)
or ``--signature`` (``score_signature_set`` / ``sc.tl.score_genes``), optional log transform.
Writes violin PNGs with Mann-Whitney U (two groups) or Kruskal-Wallis (three or more).

Example (one or more ``--h5ad``; bins are pooled across inputs)::

    python plot_gene_score_celltype_distribution.py \\
      --h5ad ../data/processed/Tuyu_Apr2026/spatial_sample_a_rabies_assignment.h5ad \\
      --gene Slc17a7 \\
      --celltype Glut

    python plot_gene_score_celltype_distribution.py \\
      --h5ad ../data/processed/Tuyu_Apr2026/spatial_*_cest_w_tumor.h5ad \\
      --gene Slc17a7 \\
      --all_celltypes \\
      --celltype_col cest_w_tumor \\
      --include_celltype GABA \\
      --include_celltype Glut \\
      --include_celltype "Microglia NN" \\
      --include_celltype tumor

    python plot_gene_score_celltype_distribution.py \\
      --h5ad ../data/processed/.../spatial_*_cest_w_tumor.h5ad \\
      --signature nrxn_signature \\
      --all_celltypes \\
      --celltype_col cest_w_tumor \\
      --include_celltype GABA --include_celltype Glut --include_celltype tumor
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
from scipy.stats import kruskal, mannwhitneyu

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "scripts"))

from gfp_neighbor_cell2location import build_celltype_color_map  # noqa: E402
from score_gaba_glut_signatures import (  # noqa: E402
    gene_obs_slug,
    infer_label_from_h5ad,
    parse_gene_args,
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
DEFAULT_CELLTYPE_COL = "cellestimate"
CELLTYPE_COL_ALIASES = {"cell_estimate": "cellestimate"}
MIN_BINS_WARN = 30


@dataclass(frozen=True)
class ScoreTarget:
    """Gene or signature set used for per-bin score values."""

    kind: str
    label: str
    file_slug: str
    set_key: str | None = None
    score_col: str | None = None


@dataclass(frozen=True)
class NonparametricTestResult:
    """Result of a Mann-Whitney U or Kruskal-Wallis test."""

    test: str
    statistic: float
    pvalue: float
    n_groups: int
    comparison: str
    group_names: tuple[str, ...]
    group_sizes: tuple[int, ...]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed CLI options.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--h5ad",
        type=Path,
        action="append",
        required=True,
        help="Spatial h5ad path (repeat to combine bins across multiple objects)",
    )
    parser.add_argument(
        "--gene",
        action="append",
        default=None,
        help="Gene symbol (exclusive with --signature; first value used if repeated)",
    )
    parser.add_argument(
        "--signature",
        choices=tuple(SIGNATURE_SETS.keys()),
        default=None,
        help="Signature set key to score (exclusive with --gene)",
    )
    parser.add_argument(
        "--celltype",
        default=None,
        help="Label in ``--celltype_col`` for single-type mode (exact string match)",
    )
    parser.add_argument(
        "--all_celltypes",
        action="store_true",
        help="One violin per cell-type label in ``--celltype_col`` (pooled h5ads)",
    )
    parser.add_argument(
        "--include_celltype",
        action="append",
        default=None,
        help=(
            "With --all_celltypes: plot only these labels (repeat flag). "
            "Exact match on --celltype_col."
        ),
    )
    parser.add_argument(
        "--celltype_col",
        default=DEFAULT_CELLTYPE_COL,
        help=f"adata.obs column for cell-type labels (default {DEFAULT_CELLTYPE_COL})",
    )
    parser.add_argument(
        "--expression_layer",
        choices=("normalized", "counts"),
        default="normalized",
        help="Expression layer for gene scoring (matches plot_score_celltype_html_v2)",
    )
    parser.add_argument(
        "--score_log",
        choices=("none", "log1p", "log10"),
        default="none",
        help="Log transform applied to scores before plotting",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Output filename stem (default: inferred from one h5ad, or combined_* for many)",
    )
    parser.add_argument(
        "--cohort",
        default=COHORT,
        help="Cohort subdirectory under results/ when --output_dir omitted",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=None,
        help="Output directory (default results/{cohort})",
    )
    return parser.parse_args()


def validate_analysis_mode(args: argparse.Namespace) -> None:
    """Ensure exactly one of single-type or all-celltypes mode is selected.

    Args:
        args: Parsed CLI namespace.

    Raises:
        ValueError: If mode flags are missing or conflicting.
    """
    if args.all_celltypes and args.celltype is not None:
        raise ValueError("Pass either --celltype or --all_celltypes, not both")
    if not args.all_celltypes and args.celltype is None:
        raise ValueError("Pass --celltype for single-type mode or --all_celltypes")
    if args.include_celltype and not args.all_celltypes:
        raise ValueError("--include_celltype requires --all_celltypes")


def signature_gene_path(set_key: str) -> Path:
    """Resolve the reference gene list path for a signature set key."""
    filename = SIGNATURE_SETS[set_key][0]
    return _REPO_ROOT / "reference" / filename


def validate_score_target(args: argparse.Namespace) -> tuple[ScoreTarget, str | None]:
    """Resolve gene or signature score target from CLI flags.

    Args:
        args: Parsed CLI namespace.

    Returns:
        tuple[ScoreTarget, str | None]: Target metadata and requested gene symbol
            (``None`` for signature mode).

    Raises:
        ValueError: If both or neither of ``--gene`` and ``--signature`` are set.
    """
    genes = parse_gene_args(args.gene) if args.gene else []
    if genes and args.signature:
        raise ValueError("Pass either --gene or --signature, not both")
    if args.signature:
        set_key = args.signature
        _, score_col = SIGNATURE_SETS[set_key]
        return (
            ScoreTarget(
                kind="signature",
                label=set_key.replace("_", " "),
                file_slug=f"signature_{label_slug(set_key)}",
                set_key=set_key,
                score_col=score_col,
            ),
            None,
        )
    if not genes:
        raise ValueError("Pass --gene or --signature")
    if len(genes) > 1:
        warnings.warn(
            f"Multiple --gene values provided; using first only: {genes[0]}",
            stacklevel=2,
        )
    symbol = genes[0]
    return (
        ScoreTarget(
            kind="gene",
            label=symbol,
            file_slug=f"gene_{gene_obs_slug(symbol)}",
        ),
        symbol,
    )


def parse_include_celltypes(args: argparse.Namespace) -> list[str] | None:
    """Parse repeatable ``--include_celltype`` values, deduped in user order.

    Args:
        args: Parsed CLI namespace.

    Returns:
        list[str] | None: Included labels, or ``None`` when unset.
    """
    raw = args.include_celltype
    if not raw:
        return None
    seen: set[str] = set()
    ordered: list[str] = []
    for label in raw:
        if label not in seen:
            seen.add(label)
            ordered.append(label)
    return ordered


def resolve_celltype_col(name: str) -> str:
    """Normalize cell-type column aliases."""
    return CELLTYPE_COL_ALIASES.get(name, name)


def label_slug(text: str) -> str:
    """Build a filesystem-safe slug from a free-text label."""
    slug = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return slug or "label"


def transform_scores_for_display(scores: np.ndarray, score_log: str) -> np.ndarray:
    """Apply optional log transform (same semantics as plot_score_celltype_html_v2)."""
    arr = np.asarray(scores, dtype=float)
    if score_log == "none":
        return arr
    if score_log == "log1p":
        return np.log1p(np.clip(arr, 0.0, None))
    if score_log == "log10":
        return np.log10(np.clip(arr, 1e-10, None))
    raise ValueError(f"Unsupported score_log {score_log!r}")


def score_axis_label(matched_var: str, expression_layer: str, score_log: str) -> str:
    """Build y-axis label for violin plots."""
    base = f"{matched_var} ({expression_layer})"
    if score_log == "log1p":
        return f"log1p({base})"
    if score_log == "log10":
        return f"log10({base})"
    return base


def ylabel_for_target(
    target: ScoreTarget,
    expression_layer: str,
    score_log: str,
) -> str:
    """Build y-axis label for gene or signature score distributions."""
    if target.kind == "gene":
        return score_axis_label(target.label, expression_layer, score_log)
    base = f"{target.label} (signature)"
    if score_log == "log1p":
        return f"log1p({base})"
    if score_log == "log10":
        return f"log10({base})"
    return base


def score_metadata_rows(
    target: ScoreTarget,
    args: argparse.Namespace,
    n_h5ad_inputs: int,
) -> dict[str, object]:
    """Common summary/test metadata for gene or signature runs."""
    layer = args.expression_layer if target.kind == "gene" else "signature"
    return {
        "score_kind": target.kind,
        "score_target": target.label,
        "gene": target.label,
        "n_h5ad_inputs": n_h5ad_inputs,
        "expression_layer": layer,
        "score_log": args.score_log,
    }


def format_pvalue(pvalue: float) -> str:
    """Format a p-value for plot labels and console output."""
    if not np.isfinite(pvalue):
        return "NA"
    if pvalue < 1e-4:
        return f"{pvalue:.2e}"
    return f"{pvalue:.4g}"


def finite_group(values: np.ndarray) -> np.ndarray:
    """Return finite score values for one group."""
    return np.asarray(values, dtype=float)[np.isfinite(values)]


def warn_small_groups(group_names: tuple[str, ...], group_sizes: tuple[int, ...]) -> None:
    """Warn when any plotted group has fewer than ``MIN_BINS_WARN`` bins."""
    for name, size in zip(group_names, group_sizes, strict=True):
        if size < MIN_BINS_WARN:
            warnings.warn(
                f"Group {name!r} has only {size} bins; distribution may be noisy.",
                stacklevel=2,
            )


def run_nonparametric_test(
    group_names: tuple[str, ...],
    groups: tuple[np.ndarray, ...],
    comparison: str,
) -> NonparametricTestResult:
    """Run Mann-Whitney U (two groups) or Kruskal-Wallis (three or more).

    Args:
        group_names: Labels for each group (same order as ``groups``).
        groups: Per-group score arrays.
        comparison: Short description stored in output CSV.

    Returns:
        NonparametricTestResult: Test name, statistic, p-value, and group metadata.

    Raises:
        ValueError: If fewer than two groups or any group is empty.
    """
    if len(group_names) != len(groups):
        raise ValueError("group_names and groups must have the same length")
    if len(groups) < 2:
        raise ValueError("At least two groups are required for a non-parametric test")

    cleaned = tuple(finite_group(g) for g in groups)
    sizes = tuple(int(g.size) for g in cleaned)
    if any(size == 0 for size in sizes):
        raise ValueError(
            "All groups must have at least one finite score; "
            f"sizes={dict(zip(group_names, sizes, strict=True))}"
        )

    if len(cleaned) == 2:
        result = mannwhitneyu(cleaned[0], cleaned[1], alternative="two-sided")
        return NonparametricTestResult(
            test="mannwhitneyu_two_sided",
            statistic=float(result.statistic),
            pvalue=float(result.pvalue),
            n_groups=2,
            comparison=comparison,
            group_names=group_names,
            group_sizes=sizes,
        )

    result = kruskal(*cleaned)
    if not np.isfinite(result.pvalue):
        raise ValueError(
            "Kruskal-Wallis returned a non-finite p-value; "
            "check for constant scores within groups."
        )
    return NonparametricTestResult(
        test="kruskalwallis",
        statistic=float(result.statistic),
        pvalue=float(result.pvalue),
        n_groups=len(cleaned),
        comparison=comparison,
        group_names=group_names,
        group_sizes=sizes,
    )


def test_annotation_text(result: NonparametricTestResult) -> str:
    """Build axes annotation text for a non-parametric test result."""
    if result.test == "mannwhitneyu_two_sided":
        header = "Mann-Whitney U"
    else:
        header = "Kruskal-Wallis H"
    return f"{header}\np = {format_pvalue(result.pvalue)}"


def annotate_test(ax: plt.Axes, result: NonparametricTestResult) -> None:
    """Add non-parametric test p-value annotation to an axes."""
    ax.text(
        0.02,
        0.98,
        test_annotation_text(result),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
    )


def style_violin_bodies(parts: dict, colors: tuple[str, ...]) -> None:
    """Apply face colors to matplotlib ``violinplot`` bodies."""
    for body, color in zip(parts["bodies"], colors, strict=True):
        body.set_facecolor(color)
        body.set_edgecolor("black")
        body.set_alpha(0.75)
    if parts.get("cmedians") is not None:
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(1.5)


def distribution_summary(values: np.ndarray) -> dict[str, float | int]:
    """Summary statistics for one score vector."""
    finite = finite_group(values)
    if finite.size == 0:
        return {
            "n": 0,
            "mean": np.nan,
            "median": np.nan,
            "std": np.nan,
            "q25": np.nan,
            "q75": np.nan,
        }
    return {
        "n": int(finite.size),
        "mean": float(np.mean(finite)),
        "median": float(np.median(finite)),
        "std": float(np.std(finite, ddof=1)) if finite.size > 1 else 0.0,
        "q25": float(np.percentile(finite, 25)),
        "q75": float(np.percentile(finite, 75)),
    }


def plot_pair_violins(
    ax: plt.Axes,
    other: np.ndarray,
    subset: np.ndarray,
    ylabel: str,
    celltype: str,
    result: NonparametricTestResult,
) -> None:
    """Violin plot: selected cell type vs bins outside that label."""
    color_other = "steelblue"
    color_sub = "coral"
    parts = ax.violinplot(
        [other, subset],
        positions=[1, 2],
        vert=True,
        showmeans=False,
        showmedians=True,
        showextrema=True,
    )
    style_violin_bodies(parts, (color_other, color_sub))
    ax.set_xticks([1, 2])
    ax.set_xticklabels(
        [
            f"Other bins\n(n={other.size})",
            f"{celltype}\n(n={subset.size})",
        ]
    )
    ax.set_ylabel(ylabel)
    ax.set_title(f"{ylabel}\n{celltype} vs other bins")
    annotate_test(ax, result)


def plot_all_celltypes_violins(
    ax: plt.Axes,
    labeled_groups: list[tuple[str, np.ndarray]],
    ylabel: str,
    celltype_col: str,
    result: NonparametricTestResult,
) -> None:
    """One violin per cell-type label."""
    labels = [name for name, _ in labeled_groups]
    arrays = [scores for _, scores in labeled_groups]
    color_map = build_celltype_color_map(labels)
    colors = tuple(color_map[label] for label in labels)
    positions = list(range(1, len(labels) + 1))

    parts = ax.violinplot(
        arrays,
        positions=positions,
        vert=True,
        showmeans=False,
        showmedians=True,
        showextrema=True,
    )
    style_violin_bodies(parts, colors)
    ax.set_xticks(positions)
    ax.set_xticklabels(
        [f"{label}\n(n={arr.size})" for label, arr in zip(labels, arrays, strict=True)],
        rotation=45,
        ha="right",
    )
    ax.set_ylabel(ylabel)
    ax.set_title(f"{ylabel}\nby {celltype_col}")
    annotate_test(ax, result)


def groups_by_celltype(
    table: pd.DataFrame,
    celltype_col: str,
    include: list[str] | None = None,
) -> list[tuple[str, np.ndarray]]:
    """Collect (label, scores) pairs from a long bin table.

    Args:
        table: Long table with ``celltype_col`` and ``score``.
        celltype_col: Observation column with cell-type labels.
        include: If set, only these labels (in this order) are returned.

    Returns:
        list[tuple[str, np.ndarray]]: Non-empty groups for plotting.

    Raises:
        ValueError: If included labels are missing or fewer than two groups remain.
    """
    available = set(table[celltype_col].astype(str).unique())
    if include is not None:
        missing = [label for label in include if label not in available]
        if missing:
            preview = ", ".join(sorted(available)[:20])
            suffix = "..." if len(available) > 20 else ""
            raise ValueError(
                f"Unknown --include_celltype label(s): {missing}. "
                f"Available in {celltype_col!r}: {preview}{suffix}"
            )
        label_order = include
    else:
        label_order = sorted(available)

    groups: list[tuple[str, np.ndarray]] = []
    for label in label_order:
        scores = table.loc[table[celltype_col] == label, "score"].to_numpy(dtype=float)
        scores = finite_group(scores)
        if scores.size == 0:
            warnings.warn(f"Skipping empty group {label!r}", stacklevel=2)
            continue
        groups.append((label, scores))
    if len(groups) < 2:
        raise ValueError(
            f"Need at least two non-empty {celltype_col} groups; found {len(groups)}"
        )
    return groups


def test_result_to_csv_row(
    result: NonparametricTestResult,
    *,
    celltype_col: str,
    target: ScoreTarget,
    n_h5ad_inputs: int,
    expression_layer: str,
    score_log: str,
    celltype: str | None = None,
    included_celltypes: list[str] | None = None,
) -> dict[str, object]:
    """Serialize a test result for CSV output."""
    row: dict[str, object] = {
        "test": result.test,
        "comparison": result.comparison,
        "statistic": result.statistic,
        "pvalue": result.pvalue,
        "n_groups": result.n_groups,
        "group_names": ";".join(result.group_names),
        "group_sizes": ";".join(str(n) for n in result.group_sizes),
        "celltype_col": celltype_col,
        "score_kind": target.kind,
        "score_target": target.label,
        "gene": target.label,
        "n_h5ad_inputs": n_h5ad_inputs,
        "expression_layer": expression_layer,
        "score_log": score_log,
    }
    if celltype is not None:
        row["celltype"] = celltype
    if included_celltypes is not None:
        row["included_celltypes"] = ";".join(included_celltypes)
    return row


def resolve_h5ad_paths(raw_paths: list[Path]) -> list[Path]:
    """Expand, resolve, and validate one or more h5ad inputs."""
    if not raw_paths:
        raise ValueError("At least one --h5ad path is required")
    resolved: list[Path] = []
    for raw in raw_paths:
        path = raw.expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"h5ad not found: {path}")
        resolved.append(path)
    return resolved


def output_label_from_h5ads(h5ad_paths: list[Path], explicit: str | None) -> str:
    """Build output filename stem for single or pooled h5ad inputs."""
    if explicit:
        return explicit
    if len(h5ad_paths) == 1:
        return infer_label_from_h5ad(h5ad_paths[0])
    sample_labels = [infer_label_from_h5ad(p) for p in h5ad_paths]
    slugs = [label_slug(s) for s in sample_labels]
    if len(slugs) <= 2:
        return "combined_" + "_".join(slugs)
    return f"combined_{len(slugs)}samples_{slugs[0]}"


def score_h5ad_bins(
    h5ad_path: Path,
    target: ScoreTarget,
    gene_symbol: str | None,
    celltype_col: str,
    expression_layer: str,
    score_log: str,
) -> tuple[pd.DataFrame, str]:
    """Score one h5ad and return a per-bin table with sample metadata."""
    sample_label = infer_label_from_h5ad(h5ad_path)
    adata = sc.read_h5ad(h5ad_path)
    if celltype_col not in adata.obs.columns:
        raise ValueError(
            f"Column {celltype_col!r} missing from {h5ad_path.name} adata.obs. "
            f"Available: {', '.join(sorted(adata.obs.columns.astype(str))[:20])}..."
        )

    score_key = target.label
    if target.kind == "signature":
        assert target.set_key is not None and target.score_col is not None
        sig_result = score_signature_set(
            adata,
            set_key=target.set_key,
            score_col=target.score_col,
            gene_path=signature_gene_path(target.set_key),
        )
        if not sig_result.scored:
            raise ValueError(
                f"Signature {target.set_key!r} could not be scored for {h5ad_path.name} "
                f"({len(sig_result.matched_genes)} genes matched)"
            )
        raw = adata.obs[target.score_col].to_numpy(dtype=float)
        score_key = target.set_key
    else:
        assert gene_symbol is not None
        result = score_single_gene(adata, gene_symbol, expression_layer)
        if not result.plotted or result.matched_var is None:
            raise ValueError(
                f"Gene {gene_symbol!r} not found in adata.var_names for {h5ad_path.name}"
            )
        raw = adata.obs[result.obs_col].to_numpy(dtype=float)
        score_key = result.matched_var

    scores = transform_scores_for_display(raw, score_log)
    obs_names = adata.obs_names.astype(str)
    prefixed_obs = [f"{sample_label}::{name}" for name in obs_names]

    table = pd.DataFrame(
        {
            "obs_name": prefixed_obs,
            "sample_label": sample_label,
            "h5ad": str(h5ad_path),
            celltype_col: adata.obs[celltype_col].astype(str).to_numpy(),
            "score": scores,
        }
    )
    return table.dropna(subset=["score"]), score_key


def combine_h5ad_tables(
    h5ad_paths: list[Path],
    target: ScoreTarget,
    gene_symbol: str | None,
    celltype_col: str,
    expression_layer: str,
    score_log: str,
) -> tuple[pd.DataFrame, ScoreTarget]:
    """Load and score each h5ad; concatenate bins into one table."""
    frames: list[pd.DataFrame] = []
    matched_keys: list[str] = []
    for path in h5ad_paths:
        frame, matched = score_h5ad_bins(
            path,
            target,
            gene_symbol,
            celltype_col,
            expression_layer,
            score_log,
        )
        frames.append(frame)
        matched_keys.append(matched)
        kind_label = "signature" if target.kind == "signature" else "gene"
        print(
            f"{path.name}: {frame.shape[0]} bins, "
            f"{kind_label}={matched}, sample={frame['sample_label'].iloc[0]}"
        )

    if target.kind == "gene":
        unique = sorted(set(matched_keys))
        if len(unique) > 1:
            warnings.warn(
                f"Matched var_names differ across h5ad inputs: {unique}; "
                f"using {unique[0]!r} for outputs.",
                stacklevel=2,
            )
        matched = unique[0]
        target = ScoreTarget(
            kind="gene",
            label=matched,
            file_slug=f"gene_{gene_obs_slug(matched)}",
        )
    return pd.concat(frames, ignore_index=True), target


def run_single_celltype_analysis(
    table: pd.DataFrame,
    args: argparse.Namespace,
    celltype_col: str,
    target: ScoreTarget,
    label: str,
    col_slug: str,
    output_dir: Path,
    h5ad_paths: list[Path],
) -> None:
    """Write violin, summary, and test outputs for one selected cell type."""
    mask = table[celltype_col] == args.celltype
    subset = finite_group(table.loc[mask, "score"].to_numpy(dtype=float))
    if subset.size == 0:
        available = sorted(table[celltype_col].unique().tolist())
        preview = ", ".join(available[:15])
        suffix = "..." if len(available) > 15 else ""
        raise ValueError(
            f"No bins with {celltype_col}={args.celltype!r}. "
            f"Available labels include: {preview}{suffix}"
        )

    other = finite_group(table.loc[~mask, "score"].to_numpy(dtype=float))
    group_names = ("other_bins", args.celltype)
    test_result = run_nonparametric_test(
        group_names,
        (other, subset),
        comparison=f"{args.celltype}_vs_other_bins",
    )
    warn_small_groups(test_result.group_names, test_result.group_sizes)

    ylabel = ylabel_for_target(target, args.expression_layer, args.score_log)
    ct_slug = label_slug(args.celltype)
    meta = score_metadata_rows(target, args, len(h5ad_paths))

    fig, ax = plt.subplots(figsize=(7, 5))
    plot_pair_violins(ax, other, subset, ylabel, args.celltype, test_result)
    plt.tight_layout()
    png_path = (
        output_dir
        / f"spatial_{label}_{target.file_slug}_{ct_slug}_{col_slug}_distribution.png"
    )
    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    summary_rows = [
        {
            "group": args.celltype,
            "celltype_col": celltype_col,
            **meta,
            **distribution_summary(subset),
        },
        {
            "group": "other_bins",
            "celltype_col": celltype_col,
            **meta,
            **distribution_summary(other),
        },
    ]
    summary_csv = (
        output_dir
        / f"spatial_{label}_{target.file_slug}_{ct_slug}_distribution_summary.csv"
    )
    pd.DataFrame(summary_rows).to_csv(summary_csv, index=False)

    test_csv = (
        output_dir
        / f"spatial_{label}_{target.file_slug}_{ct_slug}_nonparametric_test.csv"
    )
    pd.DataFrame(
        [
            test_result_to_csv_row(
                test_result,
                celltype_col=celltype_col,
                target=target,
                n_h5ad_inputs=len(h5ad_paths),
                expression_layer=str(meta["expression_layer"]),
                score_log=args.score_log,
                celltype=args.celltype,
            )
        ]
    ).to_csv(test_csv, index=False)

    print(f"Wrote {summary_csv}")
    print(f"Wrote {test_csv}")
    print(f"Wrote {png_path}")
    print(
        f"{args.celltype}: n={subset.size}, median={np.median(subset):.4f}, "
        f"mean={np.mean(subset):.4f}"
    )
    print(
        f"{test_result.test}: stat={test_result.statistic:.4f}, "
        f"p={format_pvalue(test_result.pvalue)}"
    )


def run_all_celltypes_analysis(
    table: pd.DataFrame,
    args: argparse.Namespace,
    celltype_col: str,
    target: ScoreTarget,
    label: str,
    col_slug: str,
    output_dir: Path,
    h5ad_paths: list[Path],
) -> None:
    """Write multi-violin, summary, and test outputs for all cell types."""
    include = parse_include_celltypes(args)
    labeled_groups = groups_by_celltype(table, celltype_col, include=include)
    group_names = tuple(name for name, _ in labeled_groups)
    group_arrays = tuple(scores for _, scores in labeled_groups)
    if include is not None:
        comparison = f"selected_celltypes_in_{celltype_col}"
    else:
        comparison = f"all_celltypes_in_{celltype_col}"
    test_result = run_nonparametric_test(group_names, group_arrays, comparison)
    warn_small_groups(test_result.group_names, test_result.group_sizes)

    ylabel = ylabel_for_target(target, args.expression_layer, args.score_log)
    meta = score_metadata_rows(target, args, len(h5ad_paths))
    fig_width = max(8.0, 1.2 * len(labeled_groups))
    fig, ax = plt.subplots(figsize=(fig_width, 5))
    plot_all_celltypes_violins(ax, labeled_groups, ylabel, celltype_col, test_result)
    plt.tight_layout()
    suffix = "_all_celltypes_selected" if include is not None else "_all_celltypes"
    png_path = (
        output_dir
        / f"spatial_{label}_{target.file_slug}_by_{col_slug}{suffix}_distribution.png"
    )
    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    summary_rows = [
        {
            "group": name,
            "celltype_col": celltype_col,
            **meta,
            **distribution_summary(scores),
        }
        for name, scores in labeled_groups
    ]
    summary_csv = (
        output_dir
        / f"spatial_{label}_{target.file_slug}_by_{col_slug}{suffix}_summary.csv"
    )
    pd.DataFrame(summary_rows).to_csv(summary_csv, index=False)

    test_csv = (
        output_dir
        / f"spatial_{label}_{target.file_slug}_by_{col_slug}{suffix}_nonparametric_test.csv"
    )
    pd.DataFrame(
        [
            test_result_to_csv_row(
                test_result,
                celltype_col=celltype_col,
                target=target,
                n_h5ad_inputs=len(h5ad_paths),
                expression_layer=str(meta["expression_layer"]),
                score_log=args.score_log,
                included_celltypes=include,
            )
        ]
    ).to_csv(test_csv, index=False)

    print(f"Wrote {summary_csv}")
    print(f"Wrote {test_csv}")
    print(f"Wrote {png_path}")
    print(
        f"{test_result.test} ({test_result.n_groups} groups): "
        f"stat={test_result.statistic:.4f}, p={format_pvalue(test_result.pvalue)}"
    )


def main() -> None:
    """Load h5ad(s), score gene or signature, write CSVs and distribution PNG."""
    args = parse_args()
    validate_analysis_mode(args)
    target, gene_symbol = validate_score_target(args)
    h5ad_paths = resolve_h5ad_paths(args.h5ad)

    celltype_col = resolve_celltype_col(args.celltype_col)
    label = output_label_from_h5ads(h5ad_paths, args.label)
    output_dir = args.output_dir or (_REPO_ROOT / "results" / args.cohort)
    output_dir.mkdir(parents=True, exist_ok=True)

    if len(h5ad_paths) > 1:
        print(f"Pooling {len(h5ad_paths)} h5ad inputs -> output stem {label!r}")

    table, target = combine_h5ad_tables(
        h5ad_paths,
        target,
        gene_symbol,
        celltype_col,
        args.expression_layer,
        args.score_log,
    )

    col_slug = label_slug(celltype_col)
    long_csv = output_dir / f"spatial_{label}_{target.file_slug}_by_{col_slug}.csv"
    table.to_csv(long_csv, index=False)
    print(f"Wrote {long_csv}")

    if args.all_celltypes:
        run_all_celltypes_analysis(
            table,
            args,
            celltype_col,
            target,
            label,
            col_slug,
            output_dir,
            h5ad_paths,
        )
    else:
        run_single_celltype_analysis(
            table,
            args,
            celltype_col,
            target,
            label,
            col_slug,
            output_dir,
            h5ad_paths,
        )


if __name__ == "__main__":
    main()
