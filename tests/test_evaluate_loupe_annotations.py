"""Tests for label_onehot_autocorr in evaluate_loupe_annotations.py."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_loupe_annotations import (  # noqa: E402
    build_spatial_adata,
    filter_label_categories,
    label_onehot_autocorr,
)


def _make_clustered_grid(
    n_side: int = 10,
    label_col: str = "celltype",
) -> ad.AnnData:
    """Two spatial blocks: TypeA on left, TypeB on right (strong clustering)."""
    xs, ys = np.meshgrid(np.arange(n_side), np.arange(n_side), indexing="ij")
    xs = xs.ravel()
    ys = ys.ravel()
    labels = np.where(xs < n_side // 2, "TypeA", "TypeB")
    obs = pd.DataFrame({label_col: labels})
    return build_spatial_adata(obs, xs, ys, label_col, n_neighs=6)


def _make_random_grid(
    n_side: int = 10,
    label_col: str = "celltype",
    seed: int = 0,
) -> ad.AnnData:
    """Same grid, labels shuffled (weak spatial structure)."""
    adata = _make_clustered_grid(n_side=n_side, label_col=label_col)
    rng = np.random.default_rng(seed)
    shuffled = rng.permutation(adata.obs[label_col].to_numpy())
    adata.obs[label_col] = pd.Categorical(shuffled)
    return adata


@pytest.fixture(scope="module")
def n_perms() -> int:
    # squidpy requires n_perms > 0; 1 is enough to exercise the code path.
    return 1


def test_label_onehot_autocorr_output_shape(n_perms: int) -> None:
    adata = _make_clustered_grid()
    moran, geary = label_onehot_autocorr(
        adata,
        label_col="celltype",
        min_spots=1,
        max_labels=10,
        n_perms=n_perms,
        n_jobs=1,
    )
    assert set(moran["label_category"]) == {"TypeA", "TypeB"}
    assert set(geary["label_category"]) == {"TypeA", "TypeB"}
    assert "I" in moran.columns
    assert "C" in geary.columns


def test_clustered_labels_have_higher_moran_than_shuffled(n_perms: int) -> None:
    clustered = _make_clustered_grid()
    random = _make_random_grid()

    moran_c, _ = label_onehot_autocorr(
        clustered, "celltype", min_spots=1, max_labels=10, n_perms=n_perms, n_jobs=1
    )
    moran_r, _ = label_onehot_autocorr(
        random, "celltype", min_spots=1, max_labels=10, n_perms=n_perms, n_jobs=1
    )

    i_clustered = moran_c.set_index("label_category")["I"].mean()
    i_random = moran_r.set_index("label_category")["I"].mean()
    assert i_clustered > i_random


def test_filter_drops_rare_categories(n_perms: int) -> None:
    """Rare label should collapse to Other_rare and not appear as its own Moran row."""
    xs = np.array([0, 0, 1, 1, 2, 2, 3, 3])
    ys = np.array([0, 1, 0, 1, 0, 1, 0, 1])
    labels = ["Common"] * 6 + ["Rare", "Rare"]
    obs = pd.DataFrame({"celltype": labels})
    adata = build_spatial_adata(obs, xs, ys, "celltype", n_neighs=4)

    moran, _ = label_onehot_autocorr(
        adata,
        label_col="celltype",
        min_spots=3,  # Rare (2 spots) -> Other_rare
        max_labels=10,
        n_perms=n_perms,
        n_jobs=1,
    )
    categories = set(moran["label_category"])
    assert "Rare" not in categories
    assert "Common" in categories
    assert "Other_rare" in categories


def test_filter_label_categories_collapses_rare() -> None:
    labels = pd.Series(["A", "A", "A", "B", "B"])
    filtered, kept = filter_label_categories(labels, min_spots=3, max_labels=10)
    assert kept == ["A"]
    assert set(filtered) == {"A", "Other_rare"}
    assert list(filtered).count("Other_rare") == 2
