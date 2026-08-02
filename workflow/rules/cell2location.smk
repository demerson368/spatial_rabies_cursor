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
        sr_out=lambda wildcards: sr_out_dir(config, SAMPLES.loc[wildcards.label]),
        min_counts=lambda wildcards: int(SAMPLES.loc[wildcards.label, "min_counts_c2l"]),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
    log:
        f"logs/cell2location_{{label}}_{REF_TAG}.log",
    wildcard_constraints:
        label="|".join(cell2location_labels()) if cell2location_labels() else "NEVER_MATCH",
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
            > {params.project_dir}/{log} 2>&1
        """
