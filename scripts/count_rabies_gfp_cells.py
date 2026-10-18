#!/usr/bin/env python3
"""Count cells with >=1 raw UMI on rabies/GFP reporter features (RABV_GFP_ALIASES, exact match).

Positive if any count > 0 on genes_used (case-insensitive exact alias match to var_names).

Example:
  python count_rabies_gfp_cells.py --h5ad integrated_object.h5ad
  python count_rabies_gfp_cells.py --h5ad object.h5ad --batch-key batch
  python count_rabies_gfp_cells.py --h5ad object.h5ad --annotate-h5ad
"""

from __future__ import annotations

import argparse
import re
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

RABV_GFP_ALIASES = [
    "Gfp",
    "GFP",
    "EGFP",
    "eGFP",
    "gfp",
    "egfp",
    "EGFP-L10a",
    "Gfp-L10a",
    "CVSN2c-N",
    "CVSN2c-P",
    "CVSN2c-M",
    "CVSN2c-L",
    "CVSN2c_N",
    "CVSN2c_P",
    "CVSN2c_M",
    "CVSN2c_L",
    "Nucleoprotein",
    "Phosphoprotein",
    "MatrixProtein",
    "Polymerase",
]

HINT_VAR_PATTERN = re.compile(
    r"(?i)gfp|rabies|cvsn|nucleoprotein|phosphoprotein|matrixprotein|polymerase"
)


def load_gene_list(path: Path) -> list[str]:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        genes = pd.read_excel(path).iloc[:, 0].dropna().astype(str).tolist()
    else:
        genes = path.read_text().splitlines()
    return [g.strip() for g in genes if g.strip()]


def resolve_genes_from_aliases(
    var_names: pd.Index,
    alias_list: list[str],
) -> tuple[list[str], list[str]]:
    """Match aliases to var_names by case-insensitive exact string equality."""
    lower_to_var: dict[str, str] = {}
    for v in var_names.astype(str):
        key = v.lower()
        if key not in lower_to_var:
            lower_to_var[key] = v

    found: set[str] = set()
    missing: list[str] = []
    for alias in alias_list:
        hit = lower_to_var.get(alias.lower())
        if hit is not None:
            found.add(hit)
        else:
            missing.append(alias)
    return sorted(found), missing


def resolve_count_layer(adata: sc.AnnData, layer: str | None) -> tuple[str | None, bool]:
    """Return (layer_key, used_x_fallback). layer_key None means use adata.X."""
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer not found: {layer!r}. Available: {list(adata.layers.keys())}")
        return layer, False
    for candidate in ("counts", "counts_raw"):
        if candidate in adata.layers:
            return candidate, False
    warnings.warn(
        "No layers['counts'] or layers['counts_raw']; using adata.X (may be normalized).",
        stacklevel=2,
    )
    return None, True


def _matrix_max_sample(X, n: int = 5000) -> float:
    if issparse(X):
        if X.shape[0] <= n:
            return float(X.max()) if X.nnz else 0.0
        idx = np.random.default_rng(0).choice(X.shape[0], size=n, replace=False)
        sub = X[idx]
        return float(sub.max()) if sub.nnz else 0.0
    flat = np.asarray(X).ravel()
    if flat.size > n:
        flat = np.random.default_rng(0).choice(flat, size=n, replace=False)
    return float(flat.max()) if flat.size else 0.0


def warn_if_not_raw_counts(X, used_x_fallback: bool) -> None:
    if not used_x_fallback:
        return
    mx = _matrix_max_sample(X)
    if mx <= 20:
        warnings.warn(
            f"adata.X max sample value is {mx:.4g}; counts may not be raw UMIs.",
            stacklevel=2,
        )


def get_gene_subset_matrix(adata: sc.AnnData, genes: list[str], layer: str | None):
    sub = adata[:, genes]
    if layer is not None:
        return sub.layers[layer]
    return sub.X


def cells_with_any_umi(X) -> np.ndarray:
    if issparse(X):
        return np.asarray((X > 0).getnnz(axis=1)).ravel() > 0
    return np.asarray((X > 0).any(axis=1)).ravel()


def hint_var_names(var_names: pd.Index) -> str:
    hints = [v for v in var_names.astype(str) if HINT_VAR_PATTERN.search(v)]
    if not hints:
        return "(no var_names matching rabies/GFP hint patterns)"
    shown = hints[:40]
    suffix = f" ... and {len(hints) - 40} more" if len(hints) > 40 else ""
    return ", ".join(shown) + suffix


def assign_rabies_obs(positive: np.ndarray, obs_names: pd.Index) -> pd.Series:
    labels = np.where(positive, "positive", "negative")
    return pd.Series(
        pd.Categorical(labels, categories=["negative", "positive"], ordered=True),
        index=obs_names,
        name="rabies_assignment",
    )


def build_summary(
    adata: sc.AnnData,
    positive: np.ndarray,
    batch_key: str | None,
) -> pd.DataFrame:
    n_cells = adata.n_obs
    n_pos = int(positive.sum())
    rows: list[dict] = [
        {
            "batch": "all",
            "n_cells": n_cells,
            "n_positive": n_pos,
            "pct_positive": 100.0 * n_pos / n_cells if n_cells else 0.0,
        }
    ]

    if batch_key is not None:
        batch = adata.obs[batch_key].astype(str)
        for b in sorted(batch.unique()):
            mask = batch.values == b
            n_b = int(mask.sum())
            n_p = int(positive[mask].sum())
            rows.append(
                {
                    "batch": b,
                    "n_cells": n_b,
                    "n_positive": n_p,
                    "pct_positive": 100.0 * n_p / n_b if n_b else 0.0,
                }
            )
    return pd.DataFrame(rows)


def _aggregate_positive_groups(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    out = (
        df.groupby(group_cols, observed=True)
        .agg(n_cells=("positive", "size"), n_positive=("positive", "sum"))
        .reset_index()
    )
    out["n_positive"] = out["n_positive"].astype(int)
    out["pct_positive"] = 100.0 * out["n_positive"] / out["n_cells"]
    return out


def build_summary_by_celltype(
    positive: np.ndarray,
    cell_type: pd.Series,
    batch: pd.Series | None,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Return (global batch=all by cell_type, per batch x cell_type or None)."""
    pos = np.asarray(positive, dtype=bool)
    ct = cell_type.astype(str).to_numpy()
    by_celltype_frame = pd.DataFrame(
        {
            "batch": "all",
            "cell_type": ct,
            "positive": pos,
        }
    )
    by_celltype = _aggregate_positive_groups(by_celltype_frame, ["batch", "cell_type"])
    by_celltype = by_celltype.sort_values("n_positive", ascending=False)

    if batch is None:
        return by_celltype, None

    frame = pd.DataFrame(
        {
            "batch": batch.astype(str).to_numpy(),
            "cell_type": ct,
            "positive": pos,
        }
    )
    by_batch_celltype = _aggregate_positive_groups(frame, ["batch", "cell_type"])
    by_batch_celltype = by_batch_celltype.sort_values(
        ["batch", "n_positive"],
        ascending=[True, False],
    )
    return by_celltype, by_batch_celltype


def resolve_batch_key(args: argparse.Namespace, adata: sc.AnnData) -> str | None:
    if args.no_batch and args.batch_key is not None:
        raise ValueError("Use either --batch-key or --no-batch, not both.")
    if args.no_batch:
        return None
    if args.batch_key is None:
        return None
    if args.batch_key not in adata.obs.columns:
        warnings.warn(
            f"batch-key {args.batch_key!r} not in obs; skipping batch breakdown.",
            stacklevel=2,
        )
        return None
    return args.batch_key


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count cells with >=1 UMI on rabies/GFP reporter genes (alias exact match)."
    )
    parser.add_argument("--h5ad", type=Path, required=True, help="Input .h5ad")
    parser.add_argument(
        "--batch-key",
        default=None,
        metavar="COL",
        help="obs column for per-batch breakdown (optional; omit for totals only)",
    )
    parser.add_argument(
        "--no-batch",
        action="store_true",
        help="Do not stratify by batch (same as omitting --batch-key)",
    )
    parser.add_argument(
        "--cell-type-col",
        default="cell_type",
        help="obs column for cell type (default: cell_type)",
    )
    parser.add_argument(
        "--no-celltype",
        action="store_true",
        help="Skip cell-type summary outputs",
    )
    parser.add_argument(
        "--celltype-out",
        type=Path,
        default=None,
        help="Global by-celltype CSV (default: {h5ad_stem}_rabies_gfp_by_celltype.csv)",
    )
    parser.add_argument(
        "--layer",
        default=None,
        help="Count layer (default: counts, then counts_raw, else X)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output CSV (default: {h5ad_stem}_rabies_gfp_cell_counts.csv)",
    )
    parser.add_argument(
        "--genes-out",
        type=Path,
        default=None,
        help="Gene audit file (default: {h5ad_stem}_rabies_gfp_genes_used.txt)",
    )
    parser.add_argument(
        "--genes",
        nargs="*",
        default=None,
        help=(
            f"Gene aliases to match (default: {len(RABV_GFP_ALIASES)} built-in names). "
            "Replaces defaults when provided."
        ),
    )
    parser.add_argument(
        "--genes-file",
        type=Path,
        default=None,
        help="Optional file of extra gene aliases (one per line), appended to --genes or defaults",
    )
    parser.add_argument(
        "--annotate-h5ad",
        action="store_true",
        help="Add obs rabies_assignment and write annotated h5ad",
    )
    parser.add_argument(
        "--write-h5ad",
        type=Path,
        default=None,
        help="Annotated h5ad path (default: {stem}_rabies_assignment.h5ad)",
    )
    parser.add_argument(
        "--overwrite-inplace",
        action="store_true",
        help="With --annotate-h5ad, overwrite the input --h5ad path",
    )
    parser.add_argument(
        "--backed",
        action="store_true",
        help="Open h5ad in backed mode (incompatible with --annotate-h5ad)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.annotate_h5ad and args.backed:
        raise ValueError("--annotate-h5ad cannot be used with --backed (reload h5ad in memory to write obs).")
    if args.overwrite_inplace and not args.annotate_h5ad:
        raise ValueError("--overwrite-inplace requires --annotate-h5ad.")

    h5ad_path = args.h5ad.resolve()
    if not h5ad_path.is_file():
        raise FileNotFoundError(h5ad_path)

    stem = h5ad_path.with_suffix("")
    out_path = args.out or stem.parent / f"{stem.name}_rabies_gfp_cell_counts.csv"
    genes_out_path = args.genes_out or stem.parent / f"{stem.name}_rabies_gfp_genes_used.txt"
    celltype_out_path = args.celltype_out or stem.parent / f"{stem.name}_rabies_gfp_by_celltype.csv"
    batch_celltype_out_path = (
        stem.parent / f"{stem.name}_rabies_gfp_by_batch_celltype.csv"
    )

    if args.genes is not None:
        alias_list = list(args.genes)
    else:
        alias_list = list(RABV_GFP_ALIASES)
    if args.genes_file is not None:
        alias_list.extend(load_gene_list(args.genes_file))
    seen_alias: set[str] = set()
    deduped: list[str] = []
    for a in alias_list:
        if a not in seen_alias:
            seen_alias.add(a)
            deduped.append(a)
    alias_list = deduped

    use_backed = args.backed and not args.annotate_h5ad
    adata = sc.read_h5ad(h5ad_path, backed="r" if use_backed else None)
    var_names = adata.var_names

    genes_used, aliases_missing = resolve_genes_from_aliases(var_names, alias_list)
    if not genes_used:
        raise ValueError(
            "No rabies/GFP reporter genes matched in var_names (case-insensitive exact). "
            f"Hints: {hint_var_names(var_names)}"
        )

    layer_key, used_x = resolve_count_layer(adata, args.layer)
    X = get_gene_subset_matrix(adata, genes_used, layer_key)
    warn_if_not_raw_counts(X, used_x)
    positive = cells_with_any_umi(X)

    batch_col = resolve_batch_key(args, adata)
    summary = build_summary(adata, positive, batch_col)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_path, index=False)

    audit_lines = [
        f"h5ad: {h5ad_path}",
        f"count_layer: {layer_key if layer_key is not None else 'X'}",
        f"match_rule: case_insensitive_exact",
        f"n_obs: {adata.n_obs}",
        "",
        f"genes_used ({len(genes_used)}):",
        *[f"  {g}" for g in genes_used],
        "",
        f"aliases_not_in_object ({len(aliases_missing)}):",
        *[f"  {g}" for g in aliases_missing],
    ]
    genes_out_path.write_text("\n".join(audit_lines) + "\n")

    print(f"Count layer: {layer_key if layer_key is not None else 'X'}")
    print(f"Genes used ({len(genes_used)}): {', '.join(genes_used)}")
    if aliases_missing:
        print(f"Aliases not in object ({len(aliases_missing)}): {', '.join(aliases_missing)}")
    print()
    print(summary.to_string(index=False))

    if not args.no_celltype:
        if args.cell_type_col not in adata.obs.columns:
            raise KeyError(
                f"Missing obs column {args.cell_type_col!r}. "
                f"Available: {list(adata.obs.columns[:30])}..."
            )
        cell_type_series = adata.obs[args.cell_type_col].astype(str)
        batch_series = adata.obs[batch_col].astype(str) if batch_col else None
        by_celltype, by_batch_celltype = build_summary_by_celltype(
            positive,
            cell_type_series,
            batch_series,
        )
        by_celltype.to_csv(celltype_out_path, index=False)
        print()
        print(f"Top cell types by n_positive (batch=all, col={args.cell_type_col}):")
        print(
            by_celltype.head(15)[["cell_type", "n_cells", "n_positive", "pct_positive"]].to_string(
                index=False
            )
        )
        print()
        print(f"Wrote: {celltype_out_path}")
        if by_batch_celltype is not None:
            by_batch_celltype.to_csv(batch_celltype_out_path, index=False)
            print(f"Wrote: {batch_celltype_out_path}")
        else:
            print("Skipped batch x cell-type CSV (no batch column).")

    if args.annotate_h5ad:
        if "rabies_assignment" in adata.obs.columns:
            warnings.warn("Overwriting existing obs column 'rabies_assignment'.", stacklevel=2)
        adata.obs["rabies_assignment"] = assign_rabies_obs(positive, adata.obs_names)
        adata.uns["rabies_assignment"] = {
            "genes_used": genes_used,
            "alias_list": alias_list,
            "match_rule": "case_insensitive_exact",
            "count_layer": layer_key if layer_key is not None else "X",
            "positive_rule": ">0 on any genes_used",
            "created_utc": datetime.now(timezone.utc).isoformat(),
        }
        if args.overwrite_inplace:
            write_h5ad_path = h5ad_path
        else:
            write_h5ad_path = args.write_h5ad or (stem.parent / f"{stem.name}_rabies_assignment.h5ad")
        write_h5ad_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(write_h5ad_path)
        n_pos_label = int((adata.obs["rabies_assignment"] == "positive").sum())
        print()
        print(f"rabies_assignment: {n_pos_label} positive / {adata.n_obs} cells")
        print(f"Wrote: {write_h5ad_path}")

    print()
    print(f"Wrote: {out_path}")
    print(f"Wrote: {genes_out_path}")

    if use_backed and adata.file is not None:
        adata.file.close()


if __name__ == "__main__":
    main()
