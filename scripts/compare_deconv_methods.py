#!/usr/bin/env python
"""Compare cell2location and RCTD deconvolution outputs for the same sample."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def load_dominant(adata: sc.AnnData, col: str = "dominant_celltype") -> pd.Series:
    if col not in adata.obs.columns:
        raise KeyError(f"Missing obs column: {col}")
    return adata.obs[col].astype(str)


def load_weights(adata: sc.AnnData, obsm_key: str) -> pd.DataFrame | None:
    if obsm_key not in adata.obsm:
        return None
    names = adata.uns.get("rctd_run", {}).get("cell_type_names")
    if names is None and "cell2location_run" in adata.uns:
        names = adata.uns.get("mod", {}).get("factor_names")
    matrix = np.asarray(adata.obsm[obsm_key])
    if names is None:
        names = [f"type_{i}" for i in range(matrix.shape[1])]
    return pd.DataFrame(matrix, index=adata.obs_names, columns=list(names))


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", type=str, required=True)
    parser.add_argument("--cohort", type=str, default="Tuyu_Apr2026")
    parser.add_argument("--ref_tag", type=str, default="merged_nontumor_v1_broad")
    parser.add_argument("--c2l_h5ad", type=Path, default=None)
    parser.add_argument("--rctd_h5ad", type=Path, default=None)
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=repo_root / "results" / "deconv_comparison",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]

    c2l_path = args.c2l_h5ad or (
        repo_root
        / f"results_cell2location/{args.cohort}/c2l_spatial_{args.label}_{args.ref_tag}/sp.h5ad"
    )
    rctd_path = args.rctd_h5ad or (
        repo_root
        / f"results_rctd/{args.cohort}/rctd_spatial_{args.label}_{args.ref_tag}/sp.h5ad"
    )

    adata_c2l = sc.read_h5ad(c2l_path)
    adata_rctd = sc.read_h5ad(rctd_path)

    c2l_dom = load_dominant(adata_c2l)
    rctd_dom = load_dominant(adata_rctd)
    common = c2l_dom.index.intersection(rctd_dom.index)
    merged = pd.DataFrame(
        {
            "cell2location": c2l_dom.loc[common],
            "rctd": rctd_dom.loc[common],
        }
    )
    merged["dominant_changed"] = merged["cell2location"] != merged["rctd"]

    out_dir = args.output_dir / args.label
    out_dir.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_dir / "dominant_labels_c2l_vs_rctd.csv")

    crosstab = pd.crosstab(merged["cell2location"], merged["rctd"], margins=True)
    crosstab.to_csv(out_dir / "dominant_crosstab_c2l_vs_rctd.csv")

    pct_changed = merged["dominant_changed"].mean() * 100.0
    report_lines = [
        f"label: {args.label}",
        f"cohort: {args.cohort}",
        f"ref_tag: {args.ref_tag}",
        f"c2l_h5ad: {c2l_path}",
        f"rctd_h5ad: {rctd_path}",
        f"common_spots: {len(common)}",
        f"dominant_label_changed_pct: {pct_changed:.2f}",
    ]

    c2l_weights = load_weights(adata_c2l, "means_cell_abundance_w_sf")
    rctd_weights = load_weights(adata_rctd, "rctd_weights")
    if c2l_weights is not None and rctd_weights is not None:
        shared_types = c2l_weights.columns.intersection(rctd_weights.columns)
        if len(shared_types):
            corr = (
                c2l_weights.loc[common, shared_types]
                .corrwith(rctd_weights.loc[common, shared_types], axis=0)
                .sort_values(ascending=False)
            )
            corr.to_csv(out_dir / "per_celltype_weight_correlation.csv", header=["pearson_r"])
            report_lines.extend(["", "Per-cell-type weight correlation (top 10):", corr.head(10).to_string()])

    report_path = out_dir / "comparison_report.txt"
    report_path.write_text("\n".join(report_lines) + "\n")
    print(report_path.read_text())

    top_types = merged["cell2location"].value_counts().head(8).index
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, col in zip(axes, ["cell2location", "rctd"]):
        counts = merged[col].value_counts().reindex(top_types, fill_value=0)
        counts.plot(kind="bar", ax=ax, title=f"{col} dominant counts")
        ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(out_dir / "dominant_counts_c2l_vs_rctd.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
