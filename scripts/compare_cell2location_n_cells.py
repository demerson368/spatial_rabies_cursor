#!/usr/bin/env python
"""Compare cell2location runs with different N_cells_per_location settings.

Loads two ``sp.h5ad`` outputs (typically N=2 vs N=1), reports dominant-label
agreement, secondary-type shifts on GABA-dominant bins, and optional cellestimate
concordance from a Loupe export CSV.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

COHORT = "Tuyu_Apr2026"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for N=2 vs N=1 cell2location comparison.

    Returns:
        argparse.Namespace: Parsed paths, sample label, and output directory.
    """
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", type=str, required=True, help="Sample label stem")
    parser.add_argument(
        "--ref_tag",
        type=str,
        default="merged_nontumor_v1_broad",
        help="Ref tag used in cell2location output directory names",
    )
    parser.add_argument(
        "--n2_h5ad",
        type=Path,
        default=None,
        help="N=2 sp.h5ad path (default: results_cell2location/c2l_spatial_{label}_{ref_tag}/sp.h5ad)",
    )
    parser.add_argument(
        "--n1_h5ad",
        type=Path,
        default=None,
        help="N=1 sp.h5ad path (default: .../c2l_spatial_{label}_{ref_tag}_n1/sp.h5ad)",
    )
    parser.add_argument(
        "--loupe_csv",
        type=Path,
        default=None,
        help="Optional Loupe CSV for cellestimate cross-check",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=repo_root / "results" / COHORT / "cell2location_n1_vs_n2",
        help="Directory for comparison tables",
    )
    return parser.parse_args()


def load_abundances(adata: sc.AnnData) -> pd.DataFrame:
    """Build a spot x cell-type abundance table from cell2location posterior export.

    Args:
        adata: Cell2location spatial AnnData with ``means_cell_abundance_w_sf``.

    Returns:
        pd.DataFrame: Abundances indexed by spot barcode.

    Raises:
        KeyError: If posterior abundances or factor names are missing.
    """
    if "means_cell_abundance_w_sf" not in adata.obsm:
        raise KeyError("obsm['means_cell_abundance_w_sf'] not found")
    factor_names = list(adata.uns["mod"]["factor_names"])
    matrix = adata.obsm["means_cell_abundance_w_sf"]
    if isinstance(matrix, pd.DataFrame):
        abund = matrix.copy()
        if not abund.index.equals(adata.obs_names):
            abund.index = adata.obs_names
        return abund
    return pd.DataFrame(
        np.asarray(matrix),
        index=adata.obs_names,
        columns=factor_names,
    )


def rank_summary(abund: pd.DataFrame) -> pd.DataFrame:
    """Summarize dominant and secondary cell types from abundance matrix.

    Args:
        abund: Spot x cell-type posterior mean abundances.

    Returns:
        pd.DataFrame: Columns dominant_celltype, dominant_abundance,
            secondary_celltype, secondary_abundance, margin.
    """
    vals = abund.to_numpy()
    order = np.argsort(-vals, axis=1)
    cols = abund.columns.to_numpy()
    rank1 = cols[order[:, 0]]
    rank2 = cols[order[:, 1]]
    rank1_ab = vals[np.arange(len(vals)), order[:, 0]]
    rank2_ab = vals[np.arange(len(vals)), order[:, 1]]
    return pd.DataFrame(
        {
            "dominant_celltype": rank1,
            "dominant_abundance": rank1_ab,
            "secondary_celltype": rank2,
            "secondary_abundance": rank2_ab,
            "margin": rank1_ab - rank2_ab,
        },
        index=abund.index,
    )


def secondary_counts_on_gaba(summary: pd.DataFrame) -> pd.Series:
    """Count secondary cell types among GABA-dominant spots.

    Args:
        summary: Output of :func:`rank_summary`.

    Returns:
        pd.Series: Value counts of secondary_celltype on GABA-dominant bins.
    """
    gaba = summary[summary["dominant_celltype"] == "GABA"]
    return gaba["secondary_celltype"].value_counts()


def main() -> None:
    """Load N=2 and N=1 cell2location outputs and write comparison tables."""
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]

    n2_path = args.n2_h5ad or (
        repo_root
        / f"results_cell2location/c2l_spatial_{args.label}_{args.ref_tag}/sp.h5ad"
    )
    n1_path = args.n1_h5ad or (
        repo_root
        / f"results_cell2location/c2l_spatial_{args.label}_{args.ref_tag}_n1/sp.h5ad"
    )

    adata_n2 = sc.read_h5ad(n2_path)
    adata_n1 = sc.read_h5ad(n1_path)

    abund_n2 = load_abundances(adata_n2)
    abund_n1 = load_abundances(adata_n1)
    summary_n2 = rank_summary(abund_n2)
    summary_n1 = rank_summary(abund_n1)

    common = summary_n2.index.intersection(summary_n1.index)
    merged = summary_n2.loc[common].add_suffix("_n2").join(
        summary_n1.loc[common].add_suffix("_n1")
    )
    merged["dominant_changed"] = (
        merged["dominant_celltype_n2"] != merged["dominant_celltype_n1"]
    )

    out_dir = args.output_dir / args.label
    out_dir.mkdir(parents=True, exist_ok=True)

    merged.to_csv(out_dir / "spot_summary_n2_vs_n1.csv")

    pct_changed = merged["dominant_changed"].mean() * 100.0
    crosstab = pd.crosstab(
        merged["dominant_celltype_n2"],
        merged["dominant_celltype_n1"],
        margins=True,
    )
    crosstab.to_csv(out_dir / "dominant_crosstab_n2_vs_n1.csv")

    sec_n2 = secondary_counts_on_gaba(summary_n2)
    sec_n1 = secondary_counts_on_gaba(summary_n1)
    sec_n2.to_csv(out_dir / "gaba_secondary_n2.csv", header=["count"])
    sec_n1.to_csv(out_dir / "gaba_secondary_n1.csv", header=["count"])

    report_lines = [
        f"label: {args.label}",
        f"ref_tag: {args.ref_tag}",
        f"n2_h5ad: {n2_path}",
        f"n1_h5ad: {n1_path}",
        f"common_spots: {len(common)}",
        f"dominant_label_changed_pct: {pct_changed:.2f}",
        "",
        "GABA-dominant secondary (N=2):",
        sec_n2.head(10).to_string(),
        "",
        "GABA-dominant secondary (N=1):",
        sec_n1.head(10).to_string(),
    ]

    if args.loupe_csv and args.loupe_csv.is_file():
        loupe = pd.read_csv(args.loupe_csv).set_index("Barcode")
        for suffix, summary in [("n2", summary_n2), ("n1", summary_n1)]:
            sub = summary.loc[common].copy()
            sub["cellestimate"] = loupe.reindex(sub.index)["cellestimate"]
            agree = (sub["dominant_celltype"] == sub["cellestimate"]).mean() * 100.0
            report_lines.extend(
                [
                    "",
                    f"cellestimate agreement ({suffix}, all spots): {agree:.2f}%",
                ]
            )
            gaba = sub[sub["dominant_celltype"] == "GABA"]
            if len(gaba):
                gaba_sec = gaba["cellestimate"].value_counts().head(10)
                gaba.to_csv(out_dir / f"gaba_spots_with_cellestimate_{suffix}.csv")
                report_lines.extend(
                    [
                        f"cellestimate on GABA-dominant ({suffix}):",
                        gaba_sec.to_string(),
                    ]
                )

    report_path = out_dir / "comparison_report.txt"
    report_path.write_text("\n".join(report_lines) + "\n")
    print(report_path.read_text())

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sec_n2.head(8).plot(kind="bar", ax=axes[0], title="GABA secondary (N=2)")
    sec_n1.head(8).plot(kind="bar", ax=axes[1], title="GABA secondary (N=1)")
    for ax in axes:
        ax.set_xlabel("secondary_celltype")
        ax.set_ylabel("count")
        ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(out_dir / "gaba_secondary_n2_vs_n1.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
