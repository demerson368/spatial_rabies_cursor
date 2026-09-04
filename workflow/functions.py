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
    """Filename suffix for ref-tagged annotation outputs (transfer, loupe, cellestimate, morans). Legacy v2 omits tag."""
    tag = ref_tag(config)
    return "" if tag == "v2" else f"_{tag}"


def cellestimate_filename_suffix(config: dict) -> str:
    """Fragment between spatial_{label} and _cellType.txt (includes cf; ref tag when set)."""
    cf = config["cellestimate_cf"]
    tag_suffix = ref_path_suffix(config)
    if tag_suffix:
        return f"{tag_suffix}_{cf}"
    return cf


def cellestimate_cf_arg(config: dict) -> str:
    """--cf value for cte4so.R (untagged .rds input; ref tag encoded in cf when set)."""
    cf = config["cellestimate_cf"]
    tag_suffix = ref_path_suffix(config)
    if tag_suffix:
        return f"{tag_suffix}_{cf}"
    return cf


def c2l_tag(config: dict) -> str:
    """Directory tag for cell2location outputs (always includes ref_tag)."""
    return ref_tag(config)


def _run_nontumor_filter(samples: pd.DataFrame, label: str) -> bool:
    """Return whether non-tumor filtering is enabled for a sample label."""
    if "run_nontumor_filter" not in samples.columns:
        return False
    val = samples.loc[label, "run_nontumor_filter"]
    return str(val).lower() in {"1", "true", "yes", "y"}


def rds_stem_suffix(config: dict, samples: pd.DataFrame, label: str) -> str:
    """Return '_nontumor' when filtering enabled, else ''."""
    if _run_nontumor_filter(samples, label):
        return "_nontumor"
    return ""


def annotation_path_suffix(config: dict, samples: pd.DataFrame, label: str) -> str:
    """REF_SUFFIX fragment; appends '_nontumor' when filtering enabled for this sample."""
    base = ref_path_suffix(config)
    if _run_nontumor_filter(samples, label):
        return f"{base}_nontumor"
    return base


def c2l_dir_suffix(config: dict, samples: pd.DataFrame, label: str) -> str:
    """Return '_nontumor' or '' for cell2location output directory name."""
    if _run_nontumor_filter(samples, label):
        return "_nontumor"
    return ""


def c2l_n_cells_suffix(config: dict) -> str:
    """Directory suffix when ``cell2location_n_cells_per_location`` is not 2.

    Args:
        config: Snakemake config dict.

    Returns:
        str: Empty string for default 2; otherwise ``_n{n}`` (e.g. ``_n1``).
    """
    n_cells = int(config.get("cell2location_n_cells_per_location", 2))
    if n_cells == 2:
        return ""
    return f"_n{n_cells}"


def annotation_h5ad_rel(config: dict, samples: pd.DataFrame, label: str) -> str:
    """Relative path under data/processed/ for R/CLI (cohort/filename)."""
    row = samples.loc[label]
    cohort = row["cohort"]
    stem = "_nontumor" if _run_nontumor_filter(samples, label) else ""
    return f"{cohort}/spatial_{label}{stem}.h5ad"


def annotation_h5ad_input(config: dict, samples: pd.DataFrame, label: str) -> str:
    """Repo-relative h5ad path for Snakemake inputs/outputs."""
    return f"data/processed/{annotation_h5ad_rel(config, samples, label)}"


def annotation_h5ad(config: dict, samples: pd.DataFrame, label: str) -> str:
    """Absolute preprocess or filtered h5ad path for CLI params."""
    return str(Path(config["project_dir"]) / annotation_h5ad_input(config, samples, label))


def nontumor_filter_labels(samples: pd.DataFrame) -> list[str]:
    """Labels with run_nontumor_filter == true."""
    if "run_nontumor_filter" not in samples.columns:
        return []
    mask = _as_bool(samples["run_nontumor_filter"])
    return samples.loc[mask].index.tolist()


def zr_fusion_max_score(samples: pd.DataFrame, label: str, config: dict) -> float:
    """Per-sample ZR fusion cutoff; falls back to config default."""
    default = float(config.get("zr_fusion_max_score", 0.0))
    if "zr_fusion_max_score" not in samples.columns:
        return default
    val = samples.loc[label, "zr_fusion_max_score"]
    if pd.isna(val):
        return default
    return float(val)


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
        "run_nontumor_filter",
    ):
        if col not in df.columns:
            df[col] = False
        else:
            df[col] = _as_bool(df[col])
    if "zr_fusion_max_score" not in df.columns:
        df["zr_fusion_max_score"] = float(config.get("zr_fusion_max_score", 0.0))
    return df


def sr_out_dir(config: dict, row: pd.Series) -> str:
    return f"{config['plumm_processed_root']}/{row['sr_sample_dir']}/outs"


def sr_bin_dir(config: dict, row: pd.Series) -> str:
    return f"{sr_out_dir(config, row)}/binned_outputs/square_{row['bin_um']}"


def processed_dir(config: dict, row: pd.Series) -> Path:
    return Path(config["project_dir"]) / "data" / "processed" / row["cohort"]


def h5ad_path(config: dict, label: str, samples: pd.DataFrame | None = None) -> str:
    """Base preprocess h5ad path (always unfiltered, repo-relative)."""
    if samples is None:
        samples = load_samples(config)
    row = samples.loc[label]
    return f"data/processed/{row['cohort']}/spatial_{label}.h5ad"


def nontumor_h5ad_path(config: dict, label: str, samples: pd.DataFrame | None = None) -> str:
    """Filtered non-tumor h5ad path (repo-relative)."""
    if samples is None:
        samples = load_samples(config)
    row = samples.loc[label]
    return f"data/processed/{row['cohort']}/spatial_{label}_nontumor.h5ad"


def rds_path(config: dict, label: str, samples: pd.DataFrame | None = None) -> str:
    """Seurat RDS path; collision-safe between filtered and unfiltered modes."""
    if samples is None:
        samples = load_samples(config)
    row = samples.loc[label]
    stem = rds_stem_suffix(config, samples, label)
    return f"data/processed/{row['cohort']}/spatial_{label}{stem}.rds"


def cellestimate_path(
    config: dict,
    label: str,
    samples: pd.DataFrame | None = None,
) -> str:
    """Cellestimate cellType.txt path with optional _nontumor suffix."""
    if samples is None:
        samples = load_samples(config)
    row = samples.loc[label]
    stem_suffix = cellestimate_filename_suffix(config)
    if _run_nontumor_filter(samples, label):
        return f"data/processed/{row['cohort']}/spatial_{label}_nontumor{stem_suffix}_cellType.txt"
    return f"data/processed/{row['cohort']}/spatial_{label}{stem_suffix}_cellType.txt"


def transfer_path(config: dict, label: str, samples: pd.DataFrame | None = None) -> str:
    """Label transfer output path with per-sample annotation suffix."""
    if samples is None:
        samples = load_samples(config)
    row = samples.loc[label]
    suffix = annotation_path_suffix(config, samples, label)
    return f"data/processed/{row['cohort']}/spatial_{label}{suffix}_predicted.id_transfer.txt"


def cell2location_path(
    config: dict,
    label: str,
    samples: pd.DataFrame | None = None,
) -> str:
    """Cell2location sp.h5ad path; collision-safe between filter modes."""
    if samples is None:
        samples = load_samples(config)
    tag = c2l_tag(config)
    n_cells_suffix = c2l_n_cells_suffix(config)
    dir_suffix = c2l_dir_suffix(config, samples, label)
    return f"results_cell2location/c2l_spatial_{label}_{tag}{n_cells_suffix}{dir_suffix}/sp.h5ad"


def cell2location_dir(
    config: dict,
    label: str,
    samples: pd.DataFrame | None = None,
) -> str:
    """Cell2location output directory (absolute path for CLI)."""
    rel = Path(cell2location_path(config, label, samples)).parent
    return str(Path(config["project_dir"]) / rel)


def h5ad_ext_path(config: dict, label: str, samples: pd.DataFrame | None = None) -> str:
    """Extended h5ad path from assemble_loupe."""
    if samples is None:
        samples = load_samples(config)
    row = samples.loc[label]
    suffix = annotation_path_suffix(config, samples, label)
    return f"data/processed/{row['cohort']}/spatial_{label}{suffix}_ext.h5ad"


def moran_csv_path(config: dict, label: str, samples: pd.DataFrame | None = None) -> str:
    """Moran's I CSV path with per-sample annotation suffix."""
    if samples is None:
        samples = load_samples(config)
    row = samples.loc[label]
    suffix = annotation_path_suffix(config, samples, label)
    return f"data/processed/{row['cohort']}/spatial_{label}{suffix}_moranI.csv"


def moran_svg_path(config: dict, label: str, samples: pd.DataFrame | None = None) -> str:
    """Moran's I SVG plot path with per-sample annotation suffix."""
    if samples is None:
        samples = load_samples(config)
    row = samples.loc[label]
    suffix = annotation_path_suffix(config, samples, label)
    return f"data/processed/{row['cohort']}/spatial_svgs_{label}{suffix}.png"


def loupe_path(config: dict, label: str, samples: pd.DataFrame | None = None) -> str:
    """Loupe CSV path with per-sample annotation suffix."""
    if samples is None:
        samples = load_samples(config)
    suffix = annotation_path_suffix(config, samples, label)
    return f"for_loupe/spatial_{label}{suffix}_for_loupe.csv"


def log_suffix(config: dict, samples: pd.DataFrame, label: str) -> str:
    """'_nontumor' log fragment when filtering enabled, else ''."""
    return "_nontumor" if _run_nontumor_filter(samples, label) else ""
