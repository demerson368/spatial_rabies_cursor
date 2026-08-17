#!/usr/bin/env python
"""Rank signature genes in cell2location-labeled bins and plot spatially.

Loads a Loupe export CSV and matching extended h5ad, identifies spots with a
target cell2location label (default P GABA), scores reference signature genes
by mean count or normalized expression within those bins, writes a ranked CSV,
a multi-panel spatial expression figure, and a separate bar chart of the
ranked signature genes.
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
import squidpy as sq

COHORT = "Tuyu_Apr2026"
DEFAULT_CELLTYPE_COL = "cell2location"
DEFAULT_CELLTYPE_VALUE = "P GABA"
DEFAULT_BARCODE_COL = "Barcode"
DEFAULT_GENE_SET = "p_gaba_biomarkers.txt"

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "scripts"))

from celltype_biomarker_enrichment import (  # noqa: E402
    GROUP_OBS_COL,
    build_gene_map,
    celltype_slug,
    get_count_matrix,
    get_expr_matrix,
    infer_label,
    load_gene_set,
    load_merged_adata,
    output_stem,
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for cell-type spatial gene ranking.

    Returns:
        argparse.Namespace: Parsed paths, label filters, and plot options.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--loupe_csv",
        type=Path,
        required=True,
        help="Loupe export CSV with cell2location labels",
    )
    parser.add_argument(
        "--h5ad",
        type=Path,
        required=True,
        help="Extended h5ad with expression (same barcodes as Loupe CSV)",
    )
    parser.add_argument(
        "--gene_set_file",
        type=Path,
        default=_REPO_ROOT / "reference" / DEFAULT_GENE_SET,
        help="Reference signature gene list (comma- or newline-separated)",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=_REPO_ROOT / "results" / COHORT,
        help="Directory for CSV and PNG outputs",
    )
    parser.add_argument(
        "--celltype_col",
        default=DEFAULT_CELLTYPE_COL,
        help="Loupe column with cell2location assignment",
    )
    parser.add_argument(
        "--celltype_value",
        default=DEFAULT_CELLTYPE_VALUE,
        help='Value in celltype_col defining target spots (e.g. "P GABA")',
    )
    parser.add_argument(
        "--celltype_slug",
        default=None,
        help="Slug for output filenames (default: derived from celltype_value)",
    )
    parser.add_argument(
        "--barcode_col",
        default=DEFAULT_BARCODE_COL,
        help="Barcode column shared by Loupe CSV and h5ad.obs_names",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Sample label for output filenames (default: inferred from loupe_csv)",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=12,
        help="Number of top signature genes to plot",
    )
    parser.add_argument(
        "--layer",
        choices=("counts", "normalized"),
        default="counts",
        help="Expression layer for ranking and spatial coloring",
    )
    parser.add_argument(
        "--min_pct",
        type=float,
        default=0.0,
        help="Minimum fraction of target bins expressing a gene (plot filter)",
    )
    parser.add_argument(
        "--ncols",
        type=int,
        default=3,
        help="Number of columns in the spatial multi-panel figure",
    )
    parser.add_argument(
        "--all_genes",
        action="store_true",
        help="Rank genome-wide instead of the reference signature gene set",
    )
    return parser.parse_args()


def _target_expression_stats(
    adata: sc.AnnData,
    target_value: str,
) -> tuple[int, np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    """Compute per-gene expression statistics within target-labeled bins.

    Args:
        adata: AnnData with ``obs[GROUP_OBS_COL]`` and count/expression matrices.
        target_value: Label for target spots in ``celltype_group``.

    Returns:
        tuple: ``n_target``, ``mean_count``, ``pct_expressed``, ``mean_expr``,
            and a ``var_name`` -> column index map.

    Raises:
        ValueError: If no spots are labeled with ``target_value``.
    """
    is_target = adata.obs[GROUP_OBS_COL] == target_value
    n_target = int(is_target.sum())
    if n_target == 0:
        raise ValueError(f"No spots labeled {target_value!r}")

    counts = get_count_matrix(adata)
    expr = get_expr_matrix(adata)
    target_counts = counts[is_target]
    target_expr = expr[is_target]

    mean_count = np.asarray(target_counts.mean(axis=0)).ravel()
    pct_expressed = np.asarray((target_counts > 0).mean(axis=0)).ravel()
    mean_expr = np.asarray(target_expr.mean(axis=0)).ravel()
    var_index = {g: i for i, g in enumerate(adata.var_names)}
    return n_target, mean_count, pct_expressed, mean_expr, var_index


def rank_signature_genes_in_target(
    adata: sc.AnnData,
    reference_genes: list[str],
    gene_map: dict[str, str],
    target_value: str,
    layer: str,
    min_pct: float,
) -> pd.DataFrame:
    """Rank reference signature genes by mean expression in target bins.

    Args:
        adata: AnnData with ``obs[GROUP_OBS_COL]`` and count/expression matrices.
        reference_genes: Signature gene symbols from the reference file.
        gene_map: Reference-to-panel gene symbol mapping.
        target_value: Label for target spots in ``celltype_group``.
        layer: ``counts`` or ``normalized`` — primary ranking metric.
        min_pct: Minimum fraction of target bins with count > 0 for plotting.

    Returns:
        pd.DataFrame: One row per reference gene with panel match and stats.
    """
    n_target, mean_count, pct_expressed, mean_expr, var_index = _target_expression_stats(
        adata, target_value
    )

    rows: list[dict] = []
    for ref_gene in reference_genes:
        matched = gene_map.get(ref_gene)
        if matched is None:
            rows.append(
                {
                    "gene": ref_gene,
                    "matched_var_name": None,
                    "in_panel": False,
                    "mean_count": np.nan,
                    "pct_expressed": np.nan,
                    "mean_expr": np.nan,
                    "n_target_bins": n_target,
                }
            )
            continue
        idx = var_index[matched]
        rows.append(
            {
                "gene": ref_gene,
                "matched_var_name": matched,
                "in_panel": True,
                "mean_count": float(mean_count[idx]),
                "pct_expressed": float(pct_expressed[idx]),
                "mean_expr": float(mean_expr[idx]),
                "n_target_bins": n_target,
            }
        )

    df = pd.DataFrame(rows)
    rank_col = "mean_count" if layer == "counts" else "mean_expr"
    ranked = df[df["in_panel"]].copy()
    ranked = ranked.sort_values(rank_col, ascending=False, na_position="last")
    ranked["rank"] = range(1, len(ranked) + 1)
    df = df.merge(ranked[["gene", "rank"]], on="gene", how="left")
    df["plot_eligible"] = (
        df["in_panel"]
        & df["pct_expressed"].ge(min_pct)
        & df[rank_col].gt(0)
    )
    plot_order = df[df["plot_eligible"]].sort_values(rank_col, ascending=False)
    df["_plot_order"] = df["gene"].map(
        {gene: i for i, gene in enumerate(plot_order["gene"], start=1)}
    )
    return df.sort_values(
        by=["_plot_order", "gene"],
        ascending=[True, True],
        na_position="last",
    ).drop(columns="_plot_order")


def rank_genes_mean_in_target(
    adata: sc.AnnData,
    target_value: str,
    layer: str,
    min_pct: float,
) -> pd.DataFrame:
    """Rank all genes by mean expression within target-labeled bins only.

    Args:
        adata: AnnData with ``obs[GROUP_OBS_COL]`` and count/expression matrices.
        target_value: Label for target spots in ``celltype_group``.
        layer: ``counts`` or ``normalized`` — primary ranking metric.
        min_pct: Minimum fraction of target bins with count > 0.

    Returns:
        pd.DataFrame: Ranked genes with mean_count, pct_expressed, and mean_expr.
    """
    n_target, mean_count, pct_expressed, mean_expr, _ = _target_expression_stats(
        adata, target_value
    )

    df = pd.DataFrame(
        {
            "gene": adata.var_names,
            "matched_var_name": adata.var_names,
            "in_panel": True,
            "mean_count": mean_count,
            "pct_expressed": pct_expressed,
            "mean_expr": mean_expr,
            "n_target_bins": n_target,
            "plot_eligible": (pct_expressed >= min_pct) & (mean_count > 0),
        }
    )
    rank_col = "mean_count" if layer == "counts" else "mean_expr"
    df = df[df["plot_eligible"]].copy()
    df = df.sort_values(rank_col, ascending=False).reset_index(drop=True)
    df.insert(0, "rank", df.index + 1)
    return df


def genes_to_plot(ranked: pd.DataFrame, top_n: int, layer: str) -> list[str]:
    """Select matched var names for spatial plotting from a ranked table.

    Args:
        ranked: Ranked gene table with ``plot_eligible`` and ``matched_var_name``.
        top_n: Maximum number of genes to return.
        layer: ``counts`` or ``normalized`` — used to sort when rank is missing.

    Returns:
        list[str]: Gene symbols present in ``adata.var_names`` for plotting.
    """
    rank_col = "mean_count" if layer == "counts" else "mean_expr"
    if "plot_eligible" in ranked.columns:
        candidates = ranked[ranked["plot_eligible"]].copy()
    else:
        candidates = ranked.copy()
    if "rank" in candidates.columns:
        candidates = candidates.sort_values("rank", na_position="last")
    else:
        candidates = candidates.sort_values(rank_col, ascending=False)
    genes: list[str] = []
    for matched in candidates["matched_var_name"].fillna(candidates["gene"]):
        if pd.isna(matched):
            continue
        gene = str(matched)
        if gene not in genes:
            genes.append(gene)
        if len(genes) >= top_n:
            break
    return genes


def _plot_layer(adata: sc.AnnData, layer: str) -> sc.AnnData | tuple[sc.AnnData, sc.AnnData]:
    """Return AnnData suitable for spatial plotting with the requested layer.

    Args:
        adata: Source AnnData object.
        layer: ``counts`` or ``normalized``.

    Returns:
        sc.AnnData | tuple[sc.AnnData, sc.AnnData]: Plotting AnnData, or a tuple
            of plotting AnnData and the original when a temporary copy was made.
    """
    if layer == "normalized":
        return adata
    if "counts" not in adata.layers:
        return adata
    plot_adata = adata.copy()
    plot_adata.X = adata.layers["counts"]
    return plot_adata, adata


def plot_spatial_top_genes(
    adata: sc.AnnData,
    genes: list[str],
    layer: str,
    celltype_value: str,
    out_path: Path,
    ncols: int,
    *,
    signature_mode: bool,
) -> None:
    """Save a multi-panel spatial scatter plot for ranked genes.

    Args:
        adata: AnnData with spatial coordinates and expression.
        genes: Gene symbols to plot (must exist in ``adata.var_names``).
        layer: ``counts`` or ``normalized`` for spot coloring.
        celltype_value: Target label shown in the figure title.
        out_path: Output PNG path.
        ncols: Number of columns in the panel grid.
        signature_mode: If True, title notes signature-gene evaluation.

    Returns:
        None
    """
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        raise ValueError("No requested genes found in adata.var_names")

    plot_source = _plot_layer(adata, layer)
    if isinstance(plot_source, tuple):
        plot_adata, _ = plot_source
    else:
        plot_adata = plot_source

    layer_label = "counts" if layer == "counts" else "normalized expression"
    sq.pl.spatial_scatter(
        plot_adata,
        color=genes,
        size=1.5,
        ncols=ncols,
        title=[f"{g} ({layer_label})" for g in genes],
    )
    fig = plt.gcf()
    gene_label = "signature genes" if signature_mode else "genes"
    fig.suptitle(
        f"Top {gene_label} in {celltype_value} bins — {layer_label}",
        y=1.02,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_gene_ranking_bar(
    ranked: pd.DataFrame,
    top_n: int,
    layer: str,
    celltype_value: str,
    out_path: Path,
    *,
    signature_mode: bool,
) -> None:
    """Save a horizontal bar chart of ranked genes in target bins.

    Args:
        ranked: Ranked gene table.
        top_n: Number of top genes to include.
        layer: ``counts`` or ``normalized`` — metric shown on the x-axis.
        celltype_value: Target label shown in the figure title.
        out_path: Output PNG path.
        signature_mode: If True, title notes signature-gene evaluation.

    Returns:
        None
    """
    rank_col = "mean_count" if layer == "counts" else "mean_expr"
    x_label = "Mean count" if layer == "counts" else "Mean normalized expression"
    plot_genes = genes_to_plot(ranked, top_n, layer)
    if not plot_genes:
        raise ValueError("No plottable genes after applying rank and panel filters")

    label_col = "gene"
    top = ranked[ranked["matched_var_name"].fillna(ranked["gene"]).isin(plot_genes)].copy()
    top = top.sort_values(rank_col, ascending=True)

    fig, ax = plt.subplots(figsize=(8, max(4, 0.35 * len(top))))
    bars = ax.barh(top[label_col], top[rank_col], color="#2166ac")
    ax.set_xlabel(x_label)
    ax.set_ylabel("Gene")
    title_kind = "signature genes" if signature_mode else "genes"
    ax.set_title(
        f"Top {len(top)} {title_kind} by {x_label.lower()} in {celltype_value} bins"
    )

    for bar, pct in zip(bars, top["pct_expressed"], strict=True):
        width = bar.get_width()
        pct_label = "n/a" if pd.isna(pct) else f"{pct:.0%}"
        ax.text(
            width,
            bar.get_y() + bar.get_height() / 2,
            f"  {pct_label}",
            va="center",
            ha="left",
            fontsize=8,
        )

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    """Load data, rank signature genes in target bins, and write CSV plus plots.

    Returns:
        None
    """
    args = parse_args()
    args.celltype_slug = args.celltype_slug or celltype_slug(args.celltype_value)
    label = args.label or infer_label(args.loupe_csv)
    suffix = "signature_genes" if not args.all_genes else "top_genes"
    file_stem = f"{output_stem(label, args.celltype_slug)}_{suffix}"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_merged_adata(
        args.loupe_csv,
        args.h5ad,
        args.barcode_col,
        args.celltype_col,
        args.celltype_value,
    )
    n_target = int((adata.obs[GROUP_OBS_COL] == args.celltype_value).sum())
    print(
        f"Merged {adata.n_obs:,} spots ({n_target:,} {args.celltype_value}, "
        f"{adata.n_obs - n_target:,} other)"
    )

    if args.all_genes:
        ranked = rank_genes_mean_in_target(
            adata, args.celltype_value, args.layer, args.min_pct
        )
        signature_mode = False
    else:
        reference_genes = load_gene_set(args.gene_set_file)
        gene_map = build_gene_map(reference_genes, adata.var_names.tolist())
        print(
            f"Loaded {len(reference_genes)} signature genes from {args.gene_set_file}; "
            f"matched {len(gene_map)}/{len(reference_genes)} to Visium panel"
        )
        ranked = rank_signature_genes_in_target(
            adata,
            reference_genes,
            gene_map,
            args.celltype_value,
            args.layer,
            args.min_pct,
        )
        signature_mode = True

    csv_path = args.output_dir / f"{file_stem}_ranked.csv"
    ranked.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")

    plot_table = ranked[ranked.get("plot_eligible", True)].copy()
    if "rank" in plot_table.columns:
        preview = plot_table.dropna(subset=["rank"]).sort_values("rank").head(10)
    else:
        preview = plot_table.head(10)
    print(preview.to_string(index=False))

    top_genes = genes_to_plot(ranked, args.top_n, args.layer)
    spatial_path = args.output_dir / f"{file_stem}_spatial.png"
    plot_spatial_top_genes(
        adata,
        top_genes,
        args.layer,
        args.celltype_value,
        spatial_path,
        args.ncols,
        signature_mode=signature_mode,
    )
    print(f"Saved: {spatial_path}")

    ranking_path = args.output_dir / f"{file_stem}_ranking.png"
    plot_gene_ranking_bar(
        ranked,
        args.top_n,
        args.layer,
        args.celltype_value,
        ranking_path,
        signature_mode=signature_mode,
    )
    print(f"Saved: {ranking_path}")


if __name__ == "__main__":
    main()
