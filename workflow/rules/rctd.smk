# Step 2b — RCTD (spacexr, CPU)

RCTD_REF_RDS = f"results_rctd/{{cohort}}/ref_{config['ref_label']}_{REF_TAG}/reference.rds"

RCTD_SHELL = """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.python} run_rctd.py \
            --ref_out {params.ref_out}/ \
            --ref_label {params.ref_label} \
            --ref_file {params.ref_file} \
            --ref_tag {params.ref_tag} \
            --sr_out {params.sr_out} \
            --label {wildcards.label} \
            --cohort {params.cohort} \
            --min_counts_per_bin {params.min_counts} \
            --bin_um {params.bin_um} \
            --meta_h5ad {params.meta_h5ad} \
            --spatial_dir {params.spatial_dir} \
            --doublet_mode {params.doublet_mode} \
            --max_cores {params.max_cores} \
            --umi_min {params.umi_min} \
            --umi_max {params.umi_max} \
            --rscript {params.rscript} \
            --python {params.python} \
            > {params.project_dir}/{log} 2>&1
"""


rule build_rctd_reference:
    input:
        ref_h5ad=REF_H5AD,
    output:
        reference_rds=RCTD_REF_RDS,
    params:
        ref_out=config["ref_out"],
        ref_label=config["ref_label"],
        ref_file=config["ref_file"],
        ref_tag=REF_TAG,
        ref_dir=lambda wildcards: str(
            Path(config["project_dir"])
            / f"results_rctd/{wildcards.cohort}/ref_{config['ref_label']}_{REF_TAG}"
        ),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
        rscript=RSCRIPT,
    log:
        f"logs/rctd_reference_{{cohort}}_{REF_TAG}.log",
    wildcard_constraints:
        cohort=COHORTS,
    resources:
        mem_mb=100000,
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.python} prepare_rctd_inputs.py \
            --mode reference \
            --ref_out {params.ref_out}/ \
            --ref_label {params.ref_label} \
            --ref_file {params.ref_file} \
            --ref_tag {params.ref_tag} \
            --ref_dir {params.ref_dir} \
            --rscript {params.rscript} \
            > {params.project_dir}/{log} 2>&1
        """


rule run_rctd:
    input:
        h5ad=H5AD_OUT,
        reference_rds=RCTD_REF_RDS,
    output:
        sp_h5ad=RCTD_OUT,
    params:
        ref_out=config["ref_out"],
        ref_label=config["ref_label"],
        ref_file=config["ref_file"],
        ref_tag=REF_TAG,
        cohort=lambda wildcards: wildcards.cohort,
        meta_h5ad=lambda wildcards: str(
            Path(config["project_dir"]) / H5AD_OUT.format(cohort=wildcards.cohort, label=wildcards.label)
        ),
        spatial_dir=lambda wildcards: str(
            Path(config["project_dir"])
            / Path(RCTD_OUT.format(cohort=wildcards.cohort, label=wildcards.label)).parent
        ),
        sr_out=lambda wildcards: sr_out_dir(config, SAMPLES.loc[wildcards.label]),
        min_counts=lambda wildcards: int(SAMPLES.loc[wildcards.label, "min_counts_c2l"]),
        bin_um=lambda wildcards: SAMPLES.loc[wildcards.label, "bin_um"],
        doublet_mode="doublet",
        max_cores=config.get("rctd_max_cores", 8),
        umi_min=config.get("rctd_umi_min", 0),
        umi_max=config.get("rctd_umi_max", 20000),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
        rscript=RSCRIPT,
    log:
        f"logs/rctd_{{cohort}}_{{label}}_{REF_TAG}.log",
    wildcard_constraints:
        cohort=COHORTS,
        label="|".join(unfiltered_rctd_doublet_labels()) if unfiltered_rctd_doublet_labels() else "NEVER_MATCH",
    threads: config.get("rctd_max_cores", 8),
    resources:
        mem_mb=200000,
        gpu=0,
        queue="none",
    shell:
        RCTD_SHELL


rule run_rctd_full:
    input:
        h5ad=H5AD_OUT,
        reference_rds=RCTD_REF_RDS,
    output:
        sp_h5ad=RCTD_OUT_FULL,
    params:
        ref_out=config["ref_out"],
        ref_label=config["ref_label"],
        ref_file=config["ref_file"],
        ref_tag=REF_TAG,
        cohort=lambda wildcards: wildcards.cohort,
        meta_h5ad=lambda wildcards: str(
            Path(config["project_dir"]) / H5AD_OUT.format(cohort=wildcards.cohort, label=wildcards.label)
        ),
        spatial_dir=lambda wildcards: str(
            Path(config["project_dir"])
            / Path(RCTD_OUT_FULL.format(cohort=wildcards.cohort, label=wildcards.label)).parent
        ),
        sr_out=lambda wildcards: sr_out_dir(config, SAMPLES.loc[wildcards.label]),
        min_counts=lambda wildcards: int(SAMPLES.loc[wildcards.label, "min_counts_c2l"]),
        bin_um=lambda wildcards: SAMPLES.loc[wildcards.label, "bin_um"],
        doublet_mode="full",
        max_cores=config.get("rctd_max_cores", 8),
        umi_min=config.get("rctd_umi_min", 0),
        umi_max=config.get("rctd_umi_max", 20000),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
        rscript=RSCRIPT,
    log:
        f"logs/rctd_{{cohort}}_{{label}}_{REF_TAG}_full.log",
    wildcard_constraints:
        cohort=COHORTS,
        label="|".join(unfiltered_rctd_full_labels()) if unfiltered_rctd_full_labels() else "NEVER_MATCH",
    threads: config.get("rctd_max_cores", 8),
    resources:
        mem_mb=200000,
        gpu=0,
        queue="none",
    shell:
        RCTD_SHELL


rule run_rctd_nontumor:
    input:
        h5ad=H5AD_NONTUMOR,
        reference_rds=RCTD_REF_RDS,
    output:
        sp_h5ad=RCTD_NONTUMOR_OUT,
    params:
        ref_out=config["ref_out"],
        ref_label=config["ref_label"],
        ref_file=config["ref_file"],
        ref_tag=REF_TAG,
        cohort=lambda wildcards: wildcards.cohort,
        meta_h5ad=lambda wildcards: str(
            Path(config["project_dir"])
            / H5AD_NONTUMOR.format(cohort=wildcards.cohort, label=wildcards.label)
        ),
        spatial_dir=lambda wildcards: str(
            Path(config["project_dir"])
            / Path(RCTD_NONTUMOR_OUT.format(cohort=wildcards.cohort, label=wildcards.label)).parent
        ),
        sr_out=lambda wildcards: sr_out_dir(config, SAMPLES.loc[wildcards.label]),
        min_counts=lambda wildcards: int(SAMPLES.loc[wildcards.label, "min_counts_c2l"]),
        bin_um=lambda wildcards: SAMPLES.loc[wildcards.label, "bin_um"],
        doublet_mode="doublet",
        max_cores=config.get("rctd_max_cores", 8),
        umi_min=config.get("rctd_umi_min", 0),
        umi_max=config.get("rctd_umi_max", 20000),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
        rscript=RSCRIPT,
    log:
        f"logs/rctd_{{cohort}}_{{label}}_{REF_TAG}_nontumor.log",
    wildcard_constraints:
        cohort=COHORTS,
        label="|".join(filtered_rctd_doublet_labels()) if filtered_rctd_doublet_labels() else "NEVER_MATCH",
    threads: config.get("rctd_max_cores", 8),
    resources:
        mem_mb=200000,
        gpu=0,
        queue="none",
    shell:
        RCTD_SHELL


rule run_rctd_nontumor_full:
    input:
        h5ad=H5AD_NONTUMOR,
        reference_rds=RCTD_REF_RDS,
    output:
        sp_h5ad=RCTD_NONTUMOR_OUT_FULL,
    params:
        ref_out=config["ref_out"],
        ref_label=config["ref_label"],
        ref_file=config["ref_file"],
        ref_tag=REF_TAG,
        cohort=lambda wildcards: wildcards.cohort,
        meta_h5ad=lambda wildcards: str(
            Path(config["project_dir"])
            / H5AD_NONTUMOR.format(cohort=wildcards.cohort, label=wildcards.label)
        ),
        spatial_dir=lambda wildcards: str(
            Path(config["project_dir"])
            / Path(RCTD_NONTUMOR_OUT_FULL.format(cohort=wildcards.cohort, label=wildcards.label)).parent
        ),
        sr_out=lambda wildcards: sr_out_dir(config, SAMPLES.loc[wildcards.label]),
        min_counts=lambda wildcards: int(SAMPLES.loc[wildcards.label, "min_counts_c2l"]),
        bin_um=lambda wildcards: SAMPLES.loc[wildcards.label, "bin_um"],
        doublet_mode="full",
        max_cores=config.get("rctd_max_cores", 8),
        umi_min=config.get("rctd_umi_min", 0),
        umi_max=config.get("rctd_umi_max", 20000),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
        rscript=RSCRIPT,
    log:
        f"logs/rctd_{{cohort}}_{{label}}_{REF_TAG}_full_nontumor.log",
    wildcard_constraints:
        cohort=COHORTS,
        label="|".join(filtered_rctd_full_labels()) if filtered_rctd_full_labels() else "NEVER_MATCH",
    threads: config.get("rctd_max_cores", 8),
    resources:
        mem_mb=200000,
        gpu=0,
        queue="none",
    shell:
        RCTD_SHELL
