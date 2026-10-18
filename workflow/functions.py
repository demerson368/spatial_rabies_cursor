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
        "run_rctd",
        "run_seurat_transfer",
        "run_cellestimate",
        "run_morans_i",
        "run_nontumor_filter",
    ):
        if col not in df.columns:
            df[col] = False
        else:
            df[col] = _as_bool(df[col])
    if "rctd_mode" not in df.columns:
        df["rctd_mode"] = ""
    if "zr_fusion_max_score" not in df.columns:
        df["zr_fusion_max_score"] = float(config.get("zr_fusion_max_score", 0.0))
    return df


def plumm_processed_base(config: dict) -> str:
    """Parent of cohort Processed dirs (.../common)."""
    if "plumm_processed_base" in config:
        return str(config["plumm_processed_base"]).rstrip("/")
    root = str(config["plumm_processed_root"]).rstrip("/")
    if root.endswith("/Processed"):
        return root[: -len("/Processed")].rsplit("/", 1)[0]
    return root


def cohort_for_label(samples: pd.DataFrame, label: str) -> str:
    """Return cohort string for a sample label."""
    return str(samples.loc[label, "cohort"])


def processed_dir_rel(config: dict, label: str, samples: pd.DataFrame | None = None) -> str:
    """Repo-relative processed directory for a sample (data/processed/{cohort})."""
    if samples is None:
        samples = load_samples(config)
    return f"data/processed/{cohort_for_label(samples, label)}"


def qc_plot_path(config: dict, label: str, samples: pd.DataFrame | None = None) -> str:
    """Preprocess QC plot path (repo-relative)."""
    if samples is None:
        samples = load_samples(config)
    cohort = cohort_for_label(samples, label)
    return f"data/processed/{cohort}/spatial_qc_{label}.png"


def filter_report_path(config: dict, label: str, samples: pd.DataFrame | None = None) -> str:
    """Non-tumor filter report CSV path (repo-relative)."""
    if samples is None:
        samples = load_samples(config)
    cohort = cohort_for_label(samples, label)
    return f"data/processed/{cohort}/spatial_{label}_nontumor_filter_report.csv"


def filter_qc_path(config: dict, label: str, samples: pd.DataFrame | None = None) -> str:
    """Non-tumor filter QC plot path under results/ (repo-relative)."""
    if samples is None:
        samples = load_samples(config)
    cohort = cohort_for_label(samples, label)
    return f"results/{cohort}/spatial_{label}_zr_fusion_nontumor_qc.png"


def sr_out_dir(config: dict, row: pd.Series) -> str:
    base = plumm_processed_base(config)
    return f"{base}/{row['cohort']}/Processed/{row['sr_sample_dir']}/outs"


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
    row = samples.loc[label]
    tag = c2l_tag(config)
    n_cells_suffix = c2l_n_cells_suffix(config)
    dir_suffix = c2l_dir_suffix(config, samples, label)
    cohort = row["cohort"]
    return (
        f"results_cell2location/{cohort}/c2l_spatial_{label}_{tag}"
        f"{n_cells_suffix}{dir_suffix}/sp.h5ad"
    )


def cell2location_ref_rel(config: dict, cohort: str) -> str:
    """Reference RegressionModel cache directory (repo-relative, per cohort).

    Args:
        config: Snakemake config dict.
        cohort: Cohort name from samples.tsv.

    Returns:
        str: Path like ``results_cell2location/{cohort}/ref_{ref_label}_{ref_tag}``.
    """
    ref_label = config["ref_label"]
    tag = c2l_tag(config)
    return f"results_cell2location/{cohort}/ref_{ref_label}_{tag}"


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
    row = samples.loc[label]
    suffix = annotation_path_suffix(config, samples, label)
    return f"for_loupe/{row['cohort']}/spatial_{label}{suffix}_for_loupe.csv"


def log_suffix(config: dict, samples: pd.DataFrame, label: str) -> str:
    """'_nontumor' log fragment when filtering enabled, else ''."""
    return "_nontumor" if _run_nontumor_filter(samples, label) else ""


def rctd_mode_suffix(mode: str) -> str:
    """Directory suffix for RCTD full mode outputs."""
    return "" if mode == "doublet" else "_full"


def resolve_rctd_mode(config: dict, samples: pd.DataFrame, label: str) -> str:
    """Resolve RCTD doublet/full mode for a sample label."""
    if "rctd_mode" in samples.columns:
        val = str(samples.loc[label, "rctd_mode"]).strip().lower()
        if val in {"doublet", "full"}:
            return val
    if config.get("rctd_auto_mode_by_bin", True):
        bin_um = str(samples.loc[label, "bin_um"]).strip()
        if bin_um == "016um":
            return "full"
        if bin_um == "008um":
            return "doublet"
    return str(config.get("rctd_doublet_mode", "doublet")).strip().lower()


def rctd_reference_path(config: dict, cohort: str) -> str:
    """Shared RCTD reference RDS path (repo-relative, per cohort).

    Args:
        config: Snakemake config dict.
        cohort: Cohort name from samples.tsv.

    Returns:
        str: Path to ``reference.rds`` under ``results_rctd/{cohort}/ref_...``.
    """
    ref_label = config["ref_label"]
    tag = ref_tag(config)
    return f"results_rctd/{cohort}/ref_{ref_label}_{tag}/reference.rds"


def rctd_path(
    config: dict,
    label: str,
    samples: pd.DataFrame | None = None,
) -> str:
    """RCTD sp.h5ad path; collision-safe between filter modes and RCTD modes."""
    if samples is None:
        samples = load_samples(config)
    row = samples.loc[label]
    tag = ref_tag(config)
    mode_suffix = rctd_mode_suffix(resolve_rctd_mode(config, samples, label))
    dir_suffix = c2l_dir_suffix(config, samples, label)
    cohort = row["cohort"]
    return (
        f"results_rctd/{cohort}/rctd_spatial_{label}_{tag}"
        f"{mode_suffix}{dir_suffix}/sp.h5ad"
    )


def rctd_dir(
    config: dict,
    label: str,
    samples: pd.DataFrame | None = None,
) -> str:
    """RCTD output directory (absolute path for CLI)."""
    rel = Path(rctd_path(config, label, samples)).parent
    return str(Path(config["project_dir"]) / rel)
