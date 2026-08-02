#!/usr/bin/env python
"""Permutation test for cell2location enrichment in GFP+ neighborhoods.

Builds on gfp_neighbor_cell2location.py. For each cell type and radius, compares
the observed pooled neighbor fraction around true GFP+ seeds to a null
distribution obtained by repeatedly sampling the same number of random seed
locations from a candidate spot pool (spatial permutation / label swap).

One-sided p-values:
  - p_enrichment: fraction of permutations with null fraction >= observed
  - p_depletion:  fraction of permutations with null fraction <= observed

Benjamini-Hochberg FDR is computed across cell types within each sample × radius.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import gfp_neighbor_cell2location as gnc


DEFAULT_RADII_UM = (8, 16, 24, 48, 80, 120)
COHORT = "Tuyu_Apr2026"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the permutation enrichment script.

    Returns:
        argparse.Namespace: Parsed arguments including input/output paths,
            radii, permutation count, random seed, candidate pool, and
            minimum observed neighbor count for testing.
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
        help="Directory for permutation result tables",
    )
    parser.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="Optional sample labels. Default: all spatial_*_for_loupe.csv files.",
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
        help="Bin size in microns. Default: inferred from filename.",
    )
    parser.add_argument(
        "--gfp_col",
        default=gnc.DEFAULT_GFP_COL,
        help="Column used to identify GFP+ seed spots",
    )
    parser.add_argument(
        "--gfp_value",
        default=gnc.DEFAULT_GFP_VALUE,
        help="Value in gfp_col that marks GFP+ seed spots",
    )
    parser.add_argument(
        "--c2l_col",
        default=gnc.DEFAULT_C2L_COL,
        help="cell2location assignment column",
    )
    parser.add_argument(
        "--n_perm",
        type=int,
        default=999,
        help="Number of permutations (default: 999).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for reproducible permutations.",
    )
    parser.add_argument(
        "--candidate_pool",
        choices=("all", "non_gfp"),
        default="all",
        help="Spot pool for random pseudo-seeds. "
        "'all' = any spot; 'non_gfp' = exclude GFP+ seeds from null.",
    )
    parser.add_argument(
        "--min_observed_count",
        type=int,
        default=5,
        help="Skip cell types with fewer than this many observed neighbor spots.",
    )
    return parser.parse_args()


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Compute Benjamini-Hochberg FDR-adjusted q-values.

    Args:
        p_values: One-dimensional array of raw p-values in [0, 1].

    Returns:
        np.ndarray: BH-adjusted q-values in the same order as ``p_values``.
            Returns an empty array unchanged if ``p_values`` is empty.
    """
    n = len(p_values)
    if n == 0:
        return p_values
    order = np.argsort(p_values)
    ranked = p_values[order]
    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        q[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0.0, 1.0)
    return out


def precompute_neighbor_counts(
    xs: np.ndarray,
    ys: np.ndarray,
    c2l: np.ndarray,
    grid: dict,
    radius_um: float,
    bin_um: float,
) -> tuple[np.ndarray, list[str]]:
    """Precompute neighbor cell2location counts for every spot at one radius.

    For each spot, finds all neighbors within ``radius_um`` and counts how
    many belong to each cell2location label. Used to accelerate permutation
    testing by avoiding repeated spatial lookups.

    Args:
        xs: x grid coordinates for all spots.
        ys: y grid coordinates for all spots.
        c2l: cell2location label for each spot.
        grid: Spatial index from :func:`gnc.build_spatial_index`.
        radius_um: Neighborhood radius in microns.
        bin_um: Bin size in microns.

    Returns:
        tuple containing:
            - np.ndarray: Integer matrix of shape ``(n_spots, n_cell_types)``
              with neighbor counts per spot and cell type.
            - list[str]: Cell type names corresponding to matrix columns.
    """
    cell_types = sorted(set(c2l.tolist()))
    ct_to_i = {ct: i for i, ct in enumerate(cell_types)}
    n_spots = len(xs)
    n_types = len(cell_types)
    counts = np.zeros((n_spots, n_types), dtype=np.int32)

    t0 = time.time()
    for spot_idx in range(n_spots):
        if spot_idx and spot_idx % 10000 == 0:
            elapsed = time.time() - t0
            print(f"    precompute: {spot_idx:,}/{n_spots:,} spots ({elapsed:.1f}s)")
        nbr_idxs = gnc.neighbor_indices(xs, ys, grid, spot_idx, radius_um, bin_um)
        for j in nbr_idxs:
            counts[spot_idx, ct_to_i[c2l[j]]] += 1
    print(f"    precompute done: {n_spots:,} spots in {time.time() - t0:.1f}s")
    return counts, cell_types


def pooled_fractions_from_counts(counts: np.ndarray) -> np.ndarray:
    """Convert neighbor count vector(s) to fractions summing to 1.

    Args:
        counts: Non-negative count array of any shape. Summed over all elements.

    Returns:
        np.ndarray: Same shape as ``counts`` with values divided by the total
            count. Returns zeros if the total is zero.
    """
    total = counts.sum()
    if total == 0:
        return np.zeros_like(counts, dtype=float)
    return counts / total


def permutation_test_enrichment(
    precomputed_counts: np.ndarray,
    cell_types: list[str],
    seed_indices: np.ndarray,
    candidate_indices: np.ndarray,
    global_fractions: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
    min_observed_count: int,
) -> pd.DataFrame:
    """Run spatial permutation test for cell2location neighborhood enrichment.

    Compares the observed pooled neighbor fraction for each cell type around
    true GFP+ seeds to a null distribution built by randomly sampling the same
    number of pseudo-seeds from ``candidate_indices`` on each permutation.

    Args:
        precomputed_counts: Matrix from :func:`precompute_neighbor_counts` with
            shape ``(n_spots, n_cell_types)``.
        cell_types: Cell type names corresponding to matrix columns.
        seed_indices: Indices of observed GFP+ seed spots.
        candidate_indices: Spot indices eligible for random pseudo-seeds.
        global_fractions: Tissue-wide fraction per cell type (for fold-change).
        n_perm: Number of permutations.
        rng: NumPy random generator for reproducible sampling.
        min_observed_count: Minimum observed neighbor count required to test
            a cell type.

    Returns:
        pd.DataFrame: One row per tested cell type with observed/null
            statistics, p-values, FDR q-values, and significance flags. Sorted
            by ``p_enrichment``. Empty if no cell types pass
            ``min_observed_count``.

    Raises:
        ValueError: If there are no seed spots or too few candidate spots.
    """
    n_seeds = len(seed_indices)
    if n_seeds == 0:
        raise ValueError("No GFP+ seed spots found.")
    if n_seeds > len(candidate_indices):
        raise ValueError(
            f"Need at least {n_seeds} candidate spots for permutation, "
            f"but pool has {len(candidate_indices)}."
        )

    obs_counts = precomputed_counts[seed_indices].sum(axis=0)
    obs_total = obs_counts.sum()
    obs_frac = pooled_fractions_from_counts(obs_counts.astype(float))

    null_fracs = np.zeros((n_perm, len(cell_types)), dtype=float)
    for i in range(n_perm):
        null_seeds = rng.choice(candidate_indices, size=n_seeds, replace=False)
        null_counts = precomputed_counts[null_seeds].sum(axis=0)
        null_fracs[i] = pooled_fractions_from_counts(null_counts.astype(float))

    rows: list[dict] = []
    p_enrich_list: list[float] = []
    p_deplete_list: list[float] = []
    tested_indices: list[int] = []

    for j, ct in enumerate(cell_types):
        observed_count = int(obs_counts[j])
        if observed_count < min_observed_count:
            continue

        null_col = null_fracs[:, j]
        p_enrich = (np.sum(null_col >= obs_frac[j]) + 1) / (n_perm + 1)
        p_deplete = (np.sum(null_col <= obs_frac[j]) + 1) / (n_perm + 1)
        g_frac = float(global_fractions[j])
        fold_change = obs_frac[j] / g_frac if g_frac > 0 else float("inf")

        rows.append(
            {
                "cell2location": ct,
                "observed_neighbor_count": observed_count,
                "observed_neighbor_fraction": obs_frac[j],
                "global_fraction": g_frac,
                "fold_change_vs_global": fold_change,
                "null_mean_fraction": float(null_col.mean()),
                "null_std_fraction": float(null_col.std(ddof=0)),
                "null_q025_fraction": float(np.quantile(null_col, 0.025)),
                "null_q975_fraction": float(np.quantile(null_col, 0.975)),
                "p_enrichment": p_enrich,
                "p_depletion": p_deplete,
            }
        )
        p_enrich_list.append(p_enrich)
        p_deplete_list.append(p_deplete)
        tested_indices.append(j)

    if not rows:
        return pd.DataFrame()

    q_enrich = benjamini_hochberg(np.array(p_enrich_list))
    q_deplete = benjamini_hochberg(np.array(p_deplete_list))
    for i, row in enumerate(rows):
        row["fdr_enrichment"] = q_enrich[i]
        row["fdr_depletion"] = q_deplete[i]
        row["significant_enrichment_005"] = q_enrich[i] <= 0.05
        row["significant_depletion_005"] = q_deplete[i] <= 0.05

    out = pd.DataFrame(rows)
    return out.sort_values("p_enrichment")


def analyze_sample_permutation(
    path: Path,
    output_dir: Path,
    radii_um: list[float],
    bin_um_override: float | None,
    gfp_col: str,
    gfp_value: str,
    c2l_col: str,
    n_perm: int,
    seed: int,
    candidate_pool: str,
    min_observed_count: int,
) -> pd.DataFrame:
    """Run permutation enrichment analysis for one sample across radii.

    Loads Loupe data, precomputes neighbor counts at each radius, runs
    :func:`permutation_test_enrichment`, and writes per-radius and combined
    result CSVs.

    Args:
        path: Path to a ``spatial_*_for_loupe.csv`` file.
        output_dir: Directory for permutation result tables.
        radii_um: Neighborhood radii in microns to test.
        bin_um_override: Optional bin size override; inferred if ``None``.
        gfp_col: Column marking GFP probe status.
        gfp_value: Value in ``gfp_col`` that defines GFP+ seeds.
        c2l_col: cell2location assignment column.
        n_perm: Number of permutations per radius.
        seed: Random seed for reproducible pseudo-seed sampling.
        candidate_pool: ``"all"`` uses any spot as null; ``"non_gfp"`` excludes
            GFP+ seeds from the null pool.
        min_observed_count: Minimum observed neighbor count to test a cell type.

    Returns:
        pd.DataFrame: Combined results across all radii for this sample, or an
            empty DataFrame if no cell types were testable.
    """
    label = gnc.infer_label(path)
    bin_um = gnc.infer_bin_um(path, bin_um_override)
    print(f"\nSample: {label} ({path.name}, bin={bin_um:g} µm)")

    df, xs, ys, c2l, is_seed = gnc.load_loupe_csv(
        path, bin_um, gfp_col, gfp_value, c2l_col, gnc.DEFAULT_BARCODE_COL
    )
    grid = gnc.build_spatial_index(xs, ys)
    seed_indices = np.flatnonzero(is_seed)
    n_spots = len(df)
    gfp_fraction = len(seed_indices) / n_spots
    print(f"  Spots: {n_spots:,}; GFP+ seeds: {len(seed_indices):,} ({100 * gfp_fraction:.2f}%)")

    if gfp_fraction > 0.5:
        print(
            "  Warning: >50% of spots are GFP+ seeds. "
            "Permutation null may be weakly informative; consider --candidate_pool non_gfp."
        )

    if candidate_pool == "non_gfp":
        candidate_indices = np.flatnonzero(~is_seed)
    else:
        candidate_indices = np.arange(n_spots, dtype=int)

    global_df = gnc.global_distribution(c2l)
    global_df = global_df.set_index("cell2location")

    rng = np.random.default_rng(seed)
    all_results: list[pd.DataFrame] = []

    for radius_um in radii_um:
        print(f"  radius={radius_um:g} µm: precomputing neighbor counts ...")
        precomputed, cell_types = precompute_neighbor_counts(
            xs, ys, c2l, grid, radius_um, bin_um
        )
        global_fractions = np.array(
            [global_df.loc[ct, "fraction"] if ct in global_df.index else 0.0 for ct in cell_types]
        )

        print(f"  radius={radius_um:g} µm: running {n_perm} permutations ...")
        t0 = time.time()
        result = permutation_test_enrichment(
            precomputed,
            cell_types,
            seed_indices,
            candidate_indices,
            global_fractions,
            n_perm,
            rng,
            min_observed_count,
        )
        print(f"  radius={radius_um:g} µm: done in {time.time() - t0:.1f}s")

        if result.empty:
            print(f"  radius={radius_um:g} µm: no cell types passed min_observed_count")
            continue

        n_sig = int(result["significant_enrichment_005"].sum())
        top = result.sort_values("p_enrichment").head(3)
        top_str = ", ".join(
            f"{r.cell2location} (FC={r.fold_change_vs_global:.2f}, p={r.p_enrichment:.3f})"
            for r in top.itertuples()
        )
        print(f"  radius={radius_um:g} µm: {n_sig} enriched at FDR<=0.05; top: {top_str}")

        result.insert(0, "radius_um", radius_um)
        result.insert(0, "sample", label)
        result.insert(2, "n_spots", n_spots)
        result.insert(3, "n_gfp_seeds", len(seed_indices))
        result.insert(4, "n_perm", n_perm)
        result.insert(5, "candidate_pool", candidate_pool)
        all_results.append(result)

        out_path = output_dir / (
            f"spatial_{label}_gfp_neighbor_c2l_perm_{int(radius_um)}um.csv"
        )
        gnc.write_csv(result, out_path)
        print(f"  Saved {out_path.name}")

    if not all_results:
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)
    combined_path = output_dir / f"spatial_{label}_gfp_neighbor_c2l_perm_all_radii.csv"
    gnc.write_csv(combined, combined_path)
    print(f"  Saved {combined_path.name}")
    return combined


def main() -> None:
    """Entry point: run permutation enrichment tests for all requested samples.

    Parses CLI arguments, discovers Loupe CSVs, calls
    :func:`analyze_sample_permutation` for each, and writes a combined summary
    CSV when results exist.

    Returns:
        None
    """
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    input_files = gnc.discover_input_files(args.input_dir, args.samples)
    if not input_files:
        raise FileNotFoundError(f"No spatial_*_for_loupe.csv files found in {args.input_dir}")

    print(f"Found {len(input_files)} Loupe CSV(s) in {args.input_dir}")
    print(
        f"Permutation settings: n_perm={args.n_perm}, seed={args.seed}, "
        f"candidate_pool={args.candidate_pool}"
    )

    all_combined: list[pd.DataFrame] = []
    for path in input_files:
        combined = analyze_sample_permutation(
            path,
            args.output_dir,
            args.radii_um,
            args.bin_um,
            args.gfp_col,
            args.gfp_value,
            args.c2l_col,
            args.n_perm,
            args.seed,
            args.candidate_pool,
            args.min_observed_count,
        )
        if not combined.empty:
            all_combined.append(combined)

    if all_combined:
        summary = pd.concat(all_combined, ignore_index=True)
        summary_path = args.output_dir / "all_samples_gfp_neighbor_c2l_perm_summary.csv"
        gnc.write_csv(summary, summary_path)
        print(f"\nSaved combined summary: {summary_path}")


if __name__ == "__main__":
    main()
