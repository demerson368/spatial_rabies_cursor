# Step 1 — preprocess Space Ranger binned output → normalized h5ad + QC plots

rule preprocess_spatial:
    input:
        sr_bin=lambda wildcards: sr_bin_dir(config, SAMPLES.loc[wildcards.label]),
    output:
        h5ad=H5AD_OUT,
        qc_plot=QC_PLOT,
    params:
        label=lambda wildcards: wildcards.label,
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
        runtime = "5h",
    log:
        "logs/preprocess_{label}.log",
    threads: 4
    resources:
        mem_mb=100000,
        runtime = 300    # minutes — LSF will kill after this
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.python} spatial_pipeline.py {input.sr_bin} {params.label} \
            > {params.project_dir}/{log} 2>&1
        """
