"""Helpers for spatial_rabies_cursor Snakemake workflow."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def _as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"1", "true", "yes", "y"})


def ref_tag(config: dict) -> str:
    """Suffix for ref-dependent outputs (cell2location, transfer, loupe)."""
    tag = str(config.get("ref_tag", "v2")).strip()
    if not tag:
        raise ValueError("config ref_tag must be a non-empty string (e.g. v2, compv1_mmv2)")
    return tag


def ref_path_suffix(config: dict) -> str:
    """Filename suffix for ref-dependent outputs; legacy v2 omits loupe/transfer tag."""
    tag = ref_tag(config)
    return "" if tag == "v2" else f"_{tag}"


def c2l_tag(config: dict) -> str:
    """Directory tag for cell2location outputs (always includes ref_tag)."""
    return ref_tag(config)


def load_samples(config: dict) -> pd.DataFrame:
    project_dir = Path(config["project_dir"])
    samples_path = project_dir / config["samples_tsv"]
    df = pd.read_csv(samples_path, sep="\t", comment="#")
    df.columns = df.columns.str.strip()
    if "label" not in df.columns:
        raise ValueError(
            f"{samples_path}: missing 'label' column. "
            "Do not prefix the header row with '#'; only use # for comment lines."
        )
    df = df.set_index("label", drop=False)
    for col in (
        "run_cell2location",
        "run_seurat_transfer",
        "run_cellestimate",
        "run_morans_i",
    ):
        if col not in df.columns:
            df[col] = False
        else:
            df[col] = _as_bool(df[col])
    return df


def sr_out_dir(config: dict, row: pd.Series) -> str:
    return f"{config['plumm_processed_root']}/{row['sr_sample_dir']}/outs"


def sr_bin_dir(config: dict, row: pd.Series) -> str:
    return f"{sr_out_dir(config, row)}/binned_outputs/square_{row['bin_um']}"


def processed_dir(config: dict, row: pd.Series) -> Path:
    return Path(config["project_dir"]) / "data" / "processed" / row["cohort"]


def h5ad_path(config: dict, label: str) -> str:
    row = load_samples(config).loc[label]
    return str(processed_dir(config, row) / f"spatial_{label}.h5ad")


def rds_path(config: dict, label: str) -> str:
    row = load_samples(config).loc[label]
    return str(processed_dir(config, row) / f"spatial_{label}.rds")


def cellestimate_path(config: dict, label: str) -> str:
    row = load_samples(config).loc[label]
    return str(
        processed_dir(config, row)
        / f"spatial_{label}{config['cellestimate_cf']}_cellType.txt"
    )


def transfer_path(config: dict, label: str) -> str:
    row = load_samples(config).loc[label]
    suffix = ref_path_suffix(config)
    return str(
        processed_dir(config, row)
        / f"spatial_{label}{suffix}_predicted.id_transfer.txt"
    )


def cell2location_path(config: dict, label: str) -> str:
    tag = c2l_tag(config)
    return str(
        Path(config["project_dir"])
        / "results_cell2location"
        / f"c2l_spatial_{label}_{tag}"
        / "sp.h5ad"
    )


def loupe_path(config: dict, label: str) -> str:
    suffix = ref_path_suffix(config)
    return str(
        Path(config["project_dir"]) / "for_loupe" / f"spatial_{label}{suffix}_for_loupe.csv"
    )
