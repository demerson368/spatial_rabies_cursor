"""Orchestrate RCTD deconvolution: prep inputs, run spacexr, export sp.h5ad."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent


def ensure_dir(path: str | Path) -> None:
    os.makedirs(path, exist_ok=True)


def rctd_mode_suffix(doublet_mode: str) -> str:
    return "" if doublet_mode == "doublet" else "_full"


def resolve_doublet_mode(cli_mode: str | None, bin_um: str, config_mode: str, auto_by_bin: bool) -> str:
    if cli_mode:
        return cli_mode
    if auto_by_bin:
        if bin_um == "016um":
            return "full"
        if bin_um == "008um":
            return "doublet"
    return config_mode or "doublet"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref_out", type=str, required=True)
    parser.add_argument("--ref_label", type=str, required=True)
    parser.add_argument("--ref_file", type=str, required=True)
    parser.add_argument("--ref_tag", type=str, default="v2")
    parser.add_argument("--sr_out", type=str, required=True)
    parser.add_argument("--label", type=str, required=True)
    parser.add_argument("--cohort", type=str, default=None, help="Cohort subdir under data/processed/")
    parser.add_argument("--bin_um", type=str, default="008um")
    parser.add_argument("--min_counts_per_bin", type=int, default=0)
    parser.add_argument("--meta_h5ad", type=str, default=None)
    parser.add_argument("--spatial_dir", type=str, default=None)
    parser.add_argument("--doublet_mode", type=str, default=None, choices=["doublet", "full"])
    parser.add_argument("--max_cores", type=int, default=8)
    parser.add_argument("--umi_min", type=int, default=0)
    parser.add_argument("--umi_max", type=int, default=20000)
    parser.add_argument("--rscript", type=str, default="Rscript")
    parser.add_argument("--python", type=str, default=sys.executable)
    parser.add_argument(
        "--export_only",
        action="store_true",
        help="Skip prep and run.RCTD; re-export CSV/JSON from rctd_object.rds and run h5ad export",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    doublet_mode = resolve_doublet_mode(args.doublet_mode, args.bin_um, "doublet", True)
    mode_suffix = rctd_mode_suffix(doublet_mode)

    cohort = args.cohort
    meta_h5ad = args.meta_h5ad or str(
        PROJECT_ROOT / f"data/processed/{cohort or 'Tuyu_Apr2026'}/spatial_{args.label}.h5ad"
    )
    if cohort is None and "/processed/" in meta_h5ad:
        cohort = meta_h5ad.split("/processed/")[1].split("/")[0]
    if cohort is None:
        raise ValueError(
            "Could not infer cohort; pass --cohort or use --meta_h5ad under data/processed/{cohort}/"
        )

    ref_dir = PROJECT_ROOT / f"results_rctd/{cohort}/ref_{args.ref_label}_{args.ref_tag}"
    if args.spatial_dir:
        spatial_dir = Path(args.spatial_dir)
    else:
        spatial_dir = (
            PROJECT_ROOT
            / f"results_rctd/{cohort}/rctd_spatial_{args.label}_{args.ref_tag}{mode_suffix}"
        )

    ensure_dir(ref_dir)
    ensure_dir(spatial_dir)

    reference_rds = ref_dir / "reference.rds"

    if not args.export_only:
        prep = SCRIPT_DIR / "prepare_rctd_inputs.py"
        subprocess.run(
            [
                args.python,
                str(prep),
                "--mode",
                "reference",
                "--ref_out",
                args.ref_out,
                "--ref_label",
                args.ref_label,
                "--ref_file",
                args.ref_file,
                "--ref_tag",
                args.ref_tag,
                "--ref_dir",
                str(ref_dir),
                "--rscript",
                args.rscript,
            ],
            check=True,
        )

        subprocess.run(
            [
                args.python,
                str(prep),
                "--mode",
                "spatial",
                "--sr_out",
                args.sr_out,
                "--label",
                args.label,
                "--bin_um",
                args.bin_um,
                "--meta_h5ad",
                meta_h5ad,
                "--min_counts_per_bin",
                str(args.min_counts_per_bin),
                "--spatial_dir",
                str(spatial_dir),
                "--reference_rds",
                str(reference_rds),
                "--rscript",
                args.rscript,
            ],
            check=True,
        )

    rscript_cmd = [
        args.rscript,
        str(SCRIPT_DIR / "run_rctd.R"),
        "--reference_rds",
        str(reference_rds),
        "--spatial_input_rds",
        str(spatial_dir / "spatial_input.rds"),
        "--out_dir",
        str(spatial_dir),
        "--doublet_mode",
        doublet_mode,
        "--max_cores",
        str(args.max_cores),
        "--umi_min",
        str(args.umi_min),
        "--umi_max",
        str(args.umi_max),
    ]
    if args.export_only:
        rscript_cmd.append("--export_only")
    subprocess.run(rscript_cmd, check=True)

    subprocess.run(
        [
            args.python,
            str(SCRIPT_DIR / "export_rctd_h5ad.py"),
            "--rctd_dir",
            str(spatial_dir),
            "--meta_h5ad",
            meta_h5ad,
            "--out_h5ad",
            str(spatial_dir / "sp.h5ad"),
            "--label",
            args.label,
            "--ref_label",
            args.ref_label,
            "--ref_tag",
            args.ref_tag,
            "--doublet_mode",
            doublet_mode,
            "--min_counts_per_bin",
            str(args.min_counts_per_bin),
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
