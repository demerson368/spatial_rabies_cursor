#!/usr/bin/env python
"""Quantify Allen marker enrichment at cell2location-labeled spots.

Loads a Loupe export CSV (spot annotations) and matching extended h5ad
(expression), splits spots into a target cell2location label vs all other
labels, and reports per-marker detection, differential expression,
hypergeometric ORA, custom GSEA on a reference gene set, and GO enrichment
on upregulated genes.
"""

from __future__ import annotations

import argparse
import re
import shutil
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.sparse import issparse

COHORT = "Tuyu_Apr2026"
GROUP_OBS_COL = "celltype_group"
DEFAULT_CELLTYPE_COL = "cell2location"
DEFAULT_CELLTYPE_VALUE = "P GABA"
DEFAULT_BARCODE_COL = "Barcode"
DEFAULT_GENE_SET = "p_gaba_genes_allen_yao_et_al.txt"
GO_LIBRARIES = (
    "GO_Biological_Process_2023",
    "GO_Cellular_Component_2023",
    "GO_Molecular_Function_2023",
)


def celltype_slug(celltype_value: str) -> str:
    """Convert a cell2location label to a filesystem-safe slug.

    Args:
        celltype_value: Cell type label (e.g. ``P GABA``, ``MB Dopa``).

    Returns:
        str: Lowercase slug with non-alphanumeric runs replaced by ``_``
            (e.g. ``p_gaba``, ``mb_dopa``).
    """
    slug = celltype_value.lower().strip()
    slug = re.sub(r"[^a-z0-9]+", "_", slug)
    return slug.strip("_") or "celltype"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for cell-type biomarker enrichment.

    Returns:
        argparse.Namespace: Parsed paths, column names, thresholds, and flags.
    """
    repo_root = Path(__file__).resolve().parents[1]
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
        default=repo_root / "reference" / DEFAULT_GENE_SET,
        help="Reference marker gene list (comma- or newline-separated)",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=repo_root / "results" / COHORT,
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
        help='Value in celltype_col defining target spots (e.g. "P GABA", "MB Dopa")',
    )
    parser.add_argument(
        "--celltype_slug",
        default=None,
        help="Slug for output filenames (default: derived from celltype_value)",
    )
    parser.add_argument(
        "--gsea_gene_set_name",
        default=None,
        help="Prefix for gseapy prerank gene-set name (default: {celltype_slug}_markers)",
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
        "--gsea_permutations",
        type=int,
        default=1000,
        help="Permutations for gseapy prerank GSEA",
    )
    parser.add_argument(
        "--ora_top_n",
        type=int,
        default=200,
        help="Top-N upregulated genes for hypergeometric ORA of marker set",
    )
    parser.add_argument(
        "--go_fdr",
        type=float,
        default=0.05,
        help="Adjusted p-value cutoff for GO input gene selection",
    )
    parser.add_argument(
        "--go_log2fc",
        type=float,
        default=0.25,
        help="Minimum log2 fold-change for GO input gene selection",
    )
    parser.add_argument(
        "--skip_gsea",
        action="store_true",
        help="Skip gseapy GSEA and GO steps (biomarker table + DE only)",
    )
    args = parser.parse_args()
    args.celltype_slug = args.celltype_slug or celltype_slug(args.celltype_value)
    args.gsea_gene_set_name = args.gsea_gene_set_name or f"{args.celltype_slug}_markers"
    return args


def infer_label(loupe_path: Path) -> str:
    """Infer sample label from a Loupe CSV filename.

    Args:
        loupe_path: Path like ``spatial_{label}_for_loupe.csv``.

    Returns:
        str: Sample label stem used in output filenames.
    """
    stem = loupe_path.stem
    if stem.startswith("spatial_"):
        stem = stem[len("spatial_") :]
    if stem.endswith("_for_loupe"):
        stem = stem[: -len("_for_loupe")]
    return stem


def output_stem(label: str, celltype_slug_value: str) -> str:
    """Build the shared filename stem for one sample and cell type.

    Args:
        label: Sample label.
        celltype_slug_value: Filesystem slug for the target cell type.

    Returns:
        str: Stem like ``spatial_{label}_{celltype_slug}``.
    """
    return f"spatial_{label}_{celltype_slug_value}"


def load_gene_set(path: Path) -> list[str]:
    """Load marker genes from a reference file.

    Accepts comma-separated and/or newline-separated symbols; empty tokens
    and duplicates are removed while preserving order.

    Args:
        path: Path to the gene list file.

    Returns:
        list[str]: Deduplicated gene symbols in file order.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
    """
    text = path.read_text()
    genes = [g.strip() for g in re.split(r"[,\s]+", text) if g.strip()]
    seen: set[str] = set()
    ordered: list[str] = []
    for gene in genes:
        if gene not in seen:
            seen.add(gene)
            ordered.append(gene)
    return ordered


def build_gene_map(reference_genes: list[str], var_names: list[str]) -> dict[str, str]:
    """Map reference gene symbols to ``adata.var_names`` case-insensitively.

    Args:
        reference_genes: Symbols from the reference list.
        var_names: Gene symbols in the expression matrix.

    Returns:
        dict[str, str]: Mapping from reference symbol to matched ``var_names`` entry.
    """
    exact = set(var_names)
    lower_to_actual = {g.lower(): g for g in var_names}
    mapping: dict[str, str] = {}
    for ref in reference_genes:
        if ref in exact:
            mapping[ref] = ref
        elif ref.lower() in lower_to_actual:
            mapping[ref] = lower_to_actual[ref.lower()]
    return mapping


def get_count_matrix(adata: sc.AnnData) -> np.ndarray:
    """Return a dense count matrix for expression detection stats.

    Args:
        adata: AnnData with raw counts in ``layers['counts']`` or ``.X``.

    Returns:
        np.ndarray: Count matrix shaped (n_spots, n_genes).
    """
    if "counts" in adata.layers:
        counts = adata.layers["counts"]
    else:
        counts = adata.X
    if issparse(counts):
        return counts.toarray()
    return np.asarray(counts)


def get_expr_matrix(adata: sc.AnnData) -> np.ndarray:
    """Return a dense normalized expression matrix from ``adata.X``.

    Args:
        adata: AnnData with log-normalized expression in ``.X``.

    Returns:
        np.ndarray: Expression matrix shaped (n_spots, n_genes).
    """
    if issparse(adata.X):
        return adata.X.toarray()
    return np.asarray(adata.X)


def load_merged_adata(
    loupe_csv: Path,
    h5ad_path: Path,
    barcode_col: str,
    celltype_col: str,
    celltype_value: str,
) -> sc.AnnData:
    """Load h5ad and attach Loupe labels as ``celltype_group``.

    Args:
        loupe_csv: Loupe export CSV path.
        h5ad_path: Extended h5ad path.
        barcode_col: Shared barcode column / obs index name.
        celltype_col: Column with cell2location dominant type.
        celltype_value: Value defining target test spots.

    Returns:
        sc.AnnData: Subset-aligned AnnData with ``obs[GROUP_OBS_COL]`` set.

    Raises:
        ValueError: If no barcodes overlap between Loupe CSV and h5ad.
    """
    loupe = pd.read_csv(loupe_csv)
    adata = sc.read_h5ad(h5ad_path)

    if barcode_col in loupe.columns:
        loupe = loupe.set_index(barcode_col)
    else:
        raise ValueError(f"Column {barcode_col!r} not found in {loupe_csv}")

    overlap = adata.obs_names.intersection(loupe.index)
    if len(overlap) == 0:
        raise ValueError(f"No overlapping barcodes between {loupe_csv} and {h5ad_path}")

    adata = adata[overlap].copy()
    labels = loupe.loc[adata.obs_names, celltype_col].astype(str)
    adata.obs[GROUP_OBS_COL] = np.where(
        labels == celltype_value,
        celltype_value,
        "other",
    )
    adata.obs[GROUP_OBS_COL] = adata.obs[GROUP_OBS_COL].astype("category")
    return adata


def biomarker_enrichment_table(
    adata: sc.AnnData,
    reference_genes: list[str],
    gene_map: dict[str, str],
    target_value: str,
) -> pd.DataFrame:
    """Compute per-marker detection and expression enrichment statistics.

    Args:
        adata: AnnData with ``obs[GROUP_OBS_COL]`` and expression data.
        reference_genes: Full reference gene list.
        gene_map: Reference-to-panel gene symbol mapping.
        target_value: Label for target spots in ``celltype_group``.

    Returns:
        pd.DataFrame: One row per reference gene with detection and DE stats.
    """
    is_target = adata.obs[GROUP_OBS_COL] == target_value
    counts = get_count_matrix(adata)
    expr = get_expr_matrix(adata)
    var_index = {g: i for i, g in enumerate(adata.var_names)}

    rows: list[dict] = []
    for ref_gene in reference_genes:
        in_panel = ref_gene in gene_map
        matched = gene_map.get(ref_gene)
        if in_panel and matched is not None:
            idx = var_index[matched]
            target_counts = counts[is_target, idx]
            other_counts = counts[~is_target, idx]
            target_expr = expr[is_target, idx]
            other_expr = expr[~is_target, idx]
            pct_target = float(np.mean(target_counts > 0))
            pct_other = float(np.mean(other_counts > 0))
            mean_target = float(np.mean(target_expr))
            mean_other = float(np.mean(other_expr))
            log2fc = float(np.log2((mean_target + 1e-9) / (mean_other + 1e-9)))
            combined = np.concatenate([target_expr, other_expr])
            if combined.size > 0 and np.unique(combined).size > 1:
                _, pval = stats.mannwhitneyu(
                    target_expr, other_expr, alternative="two-sided"
                )
            else:
                pval = 1.0
        else:
            matched = None
            pct_target = np.nan
            pct_other = np.nan
            mean_target = np.nan
            mean_other = np.nan
            log2fc = np.nan
            pval = np.nan

        rows.append(
            {
                "gene": ref_gene,
                "matched_var_name": matched,
                "in_panel": in_panel,
                "pct_expressed_target": pct_target,
                "pct_expressed_other": pct_other,
                "mean_expr_target": mean_target,
                "mean_expr_other": mean_other,
                "log2fc": log2fc,
                "wilcoxon_p": pval,
            }
        )

    df = pd.DataFrame(rows)
    valid = np.isfinite(df["wilcoxon_p"].to_numpy())
    fdr = np.full(len(df), np.nan)
    if valid.any():
        from statsmodels.stats.multitest import multipletests

        _, adj, _, _ = multipletests(df.loc[valid, "wilcoxon_p"], method="fdr_bh")
        fdr[valid] = adj
    df["fdr"] = fdr
    return df


def hypergeometric_ora(
    ranked_genes: pd.DataFrame,
    marker_genes_in_panel: list[str],
    top_n: int,
) -> pd.DataFrame:
    """Test over-representation of marker genes among top upregulated genes.

    Args:
        ranked_genes: Genome-wide DE table with ``gene`` and ``log2fc`` columns,
            sorted by descending ``log2fc``.
        marker_genes_in_panel: Panel-matched marker symbols to test.
        top_n: Number of top upregulated genes defining the hit set.

    Returns:
        pd.DataFrame: Single-row ORA summary with hypergeometric p-value.
    """
    universe = ranked_genes["gene"].tolist()
    universe_set = set(universe)
    markers = [g for g in marker_genes_in_panel if g in universe_set]
    top_genes = set(ranked_genes.nlargest(top_n, "log2fc")["gene"])
    hits = len(top_genes.intersection(markers))
    m_total = len(universe)
    n_markers = len(markers)
    p_value = (
        float(stats.hypergeom.sf(hits - 1, m_total, n_markers, top_n))
        if n_markers > 0
        else np.nan
    )
    return pd.DataFrame(
        [
            {
                "universe_genes": m_total,
                "marker_genes_in_panel": n_markers,
                "top_n": top_n,
                "markers_in_top_n": hits,
                "expected_markers_in_top_n": n_markers * top_n / m_total
                if m_total > 0
                else np.nan,
                "hypergeom_p": p_value,
            }
        ]
    )


def rank_genes_target_vs_other(adata: sc.AnnData, target_value: str) -> pd.DataFrame:
    """Run genome-wide Wilcoxon DE between target and other spots.

    Args:
        adata: AnnData with ``obs[GROUP_OBS_COL]``.
        target_value: Target group label.

    Returns:
        pd.DataFrame: Ranked genes with log2FC and Wilcoxon statistics.
    """
    sc.tl.rank_genes_groups(
        adata,
        groupby=GROUP_OBS_COL,
        groups=[target_value],
        reference="other",
        method="wilcoxon",
        use_raw=False,
    )
    result = adata.uns["rank_genes_groups"]
    genes = result["names"][target_value]
    logfc = result["logfoldchanges"][target_value]
    pvals = result["pvals"][target_value]
    pvals_adj = result["pvals_adj"][target_value]
    scores = result["scores"][target_value]

    df = pd.DataFrame(
        {
            "gene": genes,
            "log2fc": logfc,
            "score": scores,
            "wilcoxon_p": pvals,
            "fdr": pvals_adj,
        }
    )
    return df.sort_values("log2fc", ascending=False).reset_index(drop=True)


def run_gsea_prerank(
    ranked_genes: pd.DataFrame,
    matched_markers: list[str],
    output_prefix: Path,
    label: str,
    celltype_slug_value: str,
    gsea_gene_set_name: str,
    permutations: int,
) -> tuple[pd.DataFrame | None, Path | None, Path | None]:
    """Run prerank GSEA on a custom marker gene set.

    Args:
        ranked_genes: Genome-wide ranked DE table.
        matched_markers: Panel-matched marker symbols included in the GMT.
        output_prefix: Path prefix for GSEA result files (no suffix).
        label: Sample label used in output paths.
        celltype_slug_value: Slug for the target cell type.
        gsea_gene_set_name: Base name for the gseapy gene set.
        permutations: Number of permutations for ``gseapy.prerank``.

    Returns:
        tuple[pd.DataFrame | None, Path | None, Path | None]: GSEA result table,
            nested gseapy PDF path, and flat report PDF copy (either may be None).
    """
    if len(matched_markers) < 5:
        warnings.warn(
            f"Only {len(matched_markers)} markers in panel; skipping GSEA (need >= 5)."
        )
        return None, None, None

    import gseapy as gp

    gene_set_name = f"{gsea_gene_set_name}_{label}"
    gsea_outdir = output_prefix.parent / "prerank" / label / celltype_slug_value
    gsea_outdir.mkdir(parents=True, exist_ok=True)

    rnk = ranked_genes[["gene", "log2fc"]].dropna()
    gene_sets = {gene_set_name: matched_markers}
    pre_res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        outdir=str(gsea_outdir),
        prefix=output_prefix.name,
        permutation_num=permutations,
        seed=0,
        verbose=False,
    )

    nested_pdf = gsea_outdir / "prerank" / f"{gene_set_name}.pdf"
    report_pdf = output_prefix.parent / f"{output_prefix.name}_report.pdf"
    if nested_pdf.is_file():
        shutil.copy2(nested_pdf, report_pdf)
    else:
        warnings.warn(f"Expected gseapy PDF not found: {nested_pdf}")
        report_pdf = None

    return pre_res.res2d, nested_pdf if nested_pdf.is_file() else None, report_pdf


def run_go_enrichment(
    ranked_genes: pd.DataFrame,
    go_fdr: float,
    go_log2fc: float,
    output_dir: Path,
    file_stem: str,
) -> dict[str, pd.DataFrame]:
    """Run GO ORA via Enrichr on upregulated genes in target spots.

    Args:
        ranked_genes: Genome-wide DE table.
        go_fdr: Adjusted p-value cutoff for input genes.
        go_log2fc: Minimum log2FC for input genes.
        output_dir: Directory for GO CSV outputs.
        file_stem: Filename stem (``spatial_{label}_{celltype_slug}``).

    Returns:
        dict[str, pd.DataFrame]: GO results keyed by library short name.
    """
    import gseapy as gp

    sig = ranked_genes[
        (ranked_genes["fdr"] < go_fdr) & (ranked_genes["log2fc"] > go_log2fc)
    ]
    if len(sig) < 10:
        sig = ranked_genes.nlargest(200, "log2fc")
    gene_list = sig["gene"].dropna().astype(str).tolist()
    if len(gene_list) < 10:
        warnings.warn("Fewer than 10 genes for GO enrichment; skipping GO.")
        return {}

    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets=list(GO_LIBRARIES),
        organism="mouse",
        outdir=None,
    )
    results: dict[str, pd.DataFrame] = {}
    short_names = {"bp": GO_LIBRARIES[0], "cc": GO_LIBRARIES[1], "mf": GO_LIBRARIES[2]}
    for short, library in short_names.items():
        df = enr.results[enr.results["Gene_set"] == library].copy()
        out_path = output_dir / f"{file_stem}_go_{short}.csv"
        df.to_csv(out_path, index=False)
        results[short] = df
    return results


def plot_biomarker_dotplot(
    df: pd.DataFrame,
    out_path: Path,
    label: str,
    celltype_value: str,
) -> None:
    """Save a dotplot of marker log2FC vs detection in target spots.

    Args:
        df: Biomarker enrichment table.
        out_path: Output PNG path.
        label: Sample label for the title.
        celltype_value: Target cell type display name.

    Returns:
        None
    """
    plot_df = df[df["in_panel"]].copy().sort_values("log2fc", ascending=True)
    if plot_df.empty:
        return

    fig, ax = plt.subplots(figsize=(8, max(4, 0.18 * len(plot_df))))
    y = np.arange(len(plot_df))
    colors = np.where(plot_df["log2fc"] >= 0, "#b2182b", "#2166ac")
    ax.barh(y, plot_df["log2fc"], color=colors, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["gene"], fontsize=8)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel(f"log2 fold-change ({celltype_value} vs other)")
    ax.set_title(f"{celltype_value} markers — {label}")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_go_dotplot(
    go_df: pd.DataFrame,
    out_path: Path,
    label: str,
    celltype_value: str,
    top_n: int = 15,
) -> None:
    """Save a horizontal bar plot of top GO biological process terms.

    Args:
        go_df: GO enrichment results for one library.
        out_path: Output PNG path.
        label: Sample label for the title.
        celltype_value: Target cell type display name.
        top_n: Maximum terms to display.

    Returns:
        None
    """
    if go_df.empty:
        return
    plot_df = go_df.nsmallest(top_n, "Adjusted P-value").copy()
    plot_df["neglog10p"] = -np.log10(plot_df["Adjusted P-value"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * len(plot_df))))
    y = np.arange(len(plot_df))
    ax.barh(y, plot_df["neglog10p"], color="#4575b4")
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["Term"], fontsize=8)
    ax.set_xlabel("-log10(adjusted p-value)")
    ax.set_title(f"GO BP enrichment — {celltype_value} spots — {label}")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_gsea_summary(
    gsea_df: pd.DataFrame,
    out_path: Path,
    label: str,
    celltype_value: str,
) -> None:
    """Save a simple bar plot of GSEA NES for the custom marker gene set.

    Args:
        gsea_df: GSEA prerank results (``res2d``).
        out_path: Output PNG path.
        label: Sample label for the title.
        celltype_value: Target cell type display name.

    Returns:
        None
    """
    if gsea_df is None or gsea_df.empty:
        return
    row = gsea_df.iloc[0]
    nes = float(row.get("NES", np.nan))
    fdr = float(row.get("FDR q-val", row.get("fdr", np.nan)))
    fig, ax = plt.subplots(figsize=(4, 3))
    color = "#b2182b" if nes >= 0 else "#2166ac"
    ax.bar([f"{celltype_value} markers"], [nes], color=color)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Normalized enrichment score (NES)")
    ax.set_title(f"GSEA — {label}\nFDR={fdr:.3g}")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    """Orchestrate cell-type biomarker enrichment, GSEA, and GO analysis.

    Loads Loupe labels and expression, computes per-marker and genome-wide
    statistics, optionally runs gseapy GSEA/GO, and writes CSV/PNG outputs.

    Returns:
        None
    """
    args = parse_args()
    label = args.label or infer_label(args.loupe_csv)
    file_stem = output_stem(label, args.celltype_slug)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    reference_genes = load_gene_set(args.gene_set_file)
    print(
        f"Loaded {len(reference_genes)} reference genes for "
        f"{args.celltype_value!r} ({args.celltype_slug})"
    )

    adata = load_merged_adata(
        args.loupe_csv,
        args.h5ad,
        args.barcode_col,
        args.celltype_col,
        args.celltype_value,
    )
    n_target = int((adata.obs[GROUP_OBS_COL] == args.celltype_value).sum())
    n_other = int((adata.obs[GROUP_OBS_COL] == "other").sum())
    print(
        f"Merged {adata.n_obs:,} spots ({n_target:,} {args.celltype_value}, "
        f"{n_other:,} other)"
    )

    gene_map = build_gene_map(reference_genes, adata.var_names.tolist())
    print(f"Matched {len(gene_map)}/{len(reference_genes)} markers to Visium panel")

    biomarker_df = biomarker_enrichment_table(
        adata, reference_genes, gene_map, args.celltype_value
    )
    biomarker_path = args.output_dir / f"{file_stem}_biomarker_enrichment.csv"
    biomarker_df.to_csv(biomarker_path, index=False)
    print(f"Saved {biomarker_path}")

    ranked_genes = rank_genes_target_vs_other(adata, args.celltype_value)
    de_path = args.output_dir / f"{file_stem}_de_ranked.csv"
    ranked_genes.to_csv(de_path, index=False)
    print(f"Saved {de_path}")

    matched_panel = [gene_map[g] for g in reference_genes if g in gene_map]
    ora_df = hypergeometric_ora(ranked_genes, matched_panel, args.ora_top_n)
    ora_path = args.output_dir / f"{file_stem}_marker_ora.csv"
    ora_df.to_csv(ora_path, index=False)
    print(f"Saved {ora_path}")

    dotplot_path = args.output_dir / f"{file_stem}_biomarker_dotplot.png"
    plot_biomarker_dotplot(biomarker_df, dotplot_path, label, args.celltype_value)
    print(f"Saved {dotplot_path}")

    if args.skip_gsea:
        print("Skipping GSEA/GO (--skip_gsea)")
        return

    gsea_prefix = args.output_dir / f"{file_stem}_gsea"
    gsea_df, nested_pdf, report_pdf = run_gsea_prerank(
        ranked_genes,
        matched_panel,
        gsea_prefix,
        label,
        args.celltype_slug,
        args.gsea_gene_set_name,
        args.gsea_permutations,
    )
    if gsea_df is not None:
        gsea_csv = args.output_dir / f"{file_stem}_gsea.csv"
        gsea_df.to_csv(gsea_csv, index=False)
        print(f"Saved {gsea_csv}")
        gsea_plot = args.output_dir / f"{file_stem}_gsea.png"
        plot_gsea_summary(gsea_df, gsea_plot, label, args.celltype_value)
        print(f"Saved {gsea_plot}")
        if nested_pdf is not None:
            print(f"Saved {nested_pdf}")
        if report_pdf is not None:
            print(f"Saved {report_pdf}")

    go_results = run_go_enrichment(
        ranked_genes, args.go_fdr, args.go_log2fc, args.output_dir, file_stem
    )
    if "bp" in go_results:
        go_plot = args.output_dir / f"{file_stem}_go_dotplot.png"
        plot_go_dotplot(go_results["bp"], go_plot, label, args.celltype_value)
        print(f"Saved {go_plot}")

    print(f"{args.celltype_value} biomarker enrichment complete")


if __name__ == "__main__":
    main()
