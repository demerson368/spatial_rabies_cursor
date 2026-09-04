# Step 2 — cell2location (GPU; current script hardcodes 008um bin)

rule run_cell2location:
    input:
        h5ad=H5AD_OUT,
        ref_h5ad=REF_H5AD,
    output:
        sp_h5ad=C2L_OUT,
    params:
        ref_out=config["ref_out"],
        ref_label=config["ref_label"],
        ref_file=config["ref_file"],
        ref_tag=REF_TAG,
        meta_h5ad=lambda wildcards: str(
            Path(config["project_dir"]) / H5AD_OUT.format(label=wildcards.label)
        ),
        spatial_dir=lambda wildcards: str(
            Path(config["project_dir"])
            / f"results_cell2location/c2l_spatial_{wildcards.label}_{REF_TAG}{C2L_N_CELLS_SUFFIX}"
        ),
        sr_out=lambda wildcards: sr_out_dir(config, SAMPLES.loc[wildcards.label]),
        min_counts=lambda wildcards: int(SAMPLES.loc[wildcards.label, "min_counts_c2l"]),
        bin_um=lambda wildcards: SAMPLES.loc[wildcards.label, "bin_um"],
        n_cells=config.get("cell2location_n_cells_per_location", 2),
        detection_alpha=config.get("cell2location_detection_alpha", 200),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
    log:
        f"logs/cell2location_{{label}}_{REF_TAG}.log",
    wildcard_constraints:
        label="|".join(unfiltered_cell2location_labels()) if unfiltered_cell2location_labels() else "NEVER_MATCH",
    resources:
        mem_mb=200000,
        gpu=1,
        queue="gpu",
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.python} run_cell2location_v2.py \
            --ref_out {params.ref_out}/ \
            --ref_label {params.ref_label} \
            --ref_file {params.ref_file} \
            --ref_tag {params.ref_tag} \
            --sr_out {params.sr_out} \
            --label {wildcards.label} \
            --min_counts_per_bin {params.min_counts} \
            --bin_um {params.bin_um} \
            --meta_h5ad {params.meta_h5ad} \
            --spatial_dir {params.spatial_dir} \
            --n_cells_per_location {params.n_cells} \
            --detection_alpha {params.detection_alpha} \
            > {params.project_dir}/{log} 2>&1
        """


rule run_cell2location_nontumor:
    input:
        h5ad=H5AD_NONTUMOR,
        ref_h5ad=REF_H5AD,
    output:
        sp_h5ad=C2L_NONTUMOR_OUT,
    params:
        ref_out=config["ref_out"],
        ref_label=config["ref_label"],
        ref_file=config["ref_file"],
        ref_tag=REF_TAG,
        meta_h5ad=lambda wildcards: str(
            Path(config["project_dir"]) / H5AD_NONTUMOR.format(label=wildcards.label)
        ),
        spatial_dir=lambda wildcards: str(
            Path(config["project_dir"])
            / f"results_cell2location/c2l_spatial_{wildcards.label}_{REF_TAG}{C2L_N_CELLS_SUFFIX}_nontumor"
        ),
        sr_out=lambda wildcards: sr_out_dir(config, SAMPLES.loc[wildcards.label]),
        min_counts=lambda wildcards: int(SAMPLES.loc[wildcards.label, "min_counts_c2l"]),
        bin_um=lambda wildcards: SAMPLES.loc[wildcards.label, "bin_um"],
        n_cells=config.get("cell2location_n_cells_per_location", 2),
        detection_alpha=config.get("cell2location_detection_alpha", 200),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
    log:
        f"logs/cell2location_{{label}}_{REF_TAG}_nontumor.log",
    wildcard_constraints:
        label="|".join(filtered_cell2location_labels()) if filtered_cell2location_labels() else "NEVER_MATCH",
    resources:
        mem_mb=200000,
        gpu=1,
        queue="gpu",
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.python} run_cell2location_v2.py \
            --ref_out {params.ref_out}/ \
            --ref_label {params.ref_label} \
            --ref_file {params.ref_file} \
            --ref_tag {params.ref_tag} \
            --sr_out {params.sr_out} \
            --label {wildcards.label} \
            --min_counts_per_bin {params.min_counts} \
            --bin_um {params.bin_um} \
            --meta_h5ad {params.meta_h5ad} \
            --spatial_dir {params.spatial_dir} \
            --n_cells_per_location {params.n_cells} \
            --detection_alpha {params.detection_alpha} \
            > {params.project_dir}/{log} 2>&1
        """
