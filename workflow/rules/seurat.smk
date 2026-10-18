# Step 3a — h5ad → Seurat rds (required for label transfer + Cellestimate)

rule create_seurat:
    input:
        h5ad=H5AD_OUT,
    output:
        rds=RDS_OUT,
    params:
        h5ad_rel=lambda wildcards: H5AD_OUT.format(cohort=wildcards.cohort, label=wildcards.label).replace("data/processed/", ""),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        rscript=RSCRIPT,
    log:
        "logs/create_seurat_{cohort}_{label}.log",
    wildcard_constraints:
        cohort=COHORTS,
        label="|".join(unfiltered_seurat_labels()) if unfiltered_seurat_labels() else "NEVER_MATCH",
    resources:
        mem_mb=50000,
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.rscript} create_seurat_from_h5ad.R --h5ad ../data/processed/{params.h5ad_rel} \
            > {params.project_dir}/{log} 2>&1
        """


rule create_seurat_nontumor:
    input:
        h5ad=H5AD_NONTUMOR,
    output:
        rds=RDS_NONTUMOR,
    params:
        h5ad_rel=lambda wildcards: H5AD_NONTUMOR.format(cohort=wildcards.cohort, label=wildcards.label).replace("data/processed/", ""),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        rscript=RSCRIPT,
    log:
        "logs/create_seurat_{cohort}_{label}_nontumor.log",
    wildcard_constraints:
        cohort=COHORTS,
        label="|".join(filtered_seurat_labels()) if filtered_seurat_labels() else "NEVER_MATCH",
    resources:
        mem_mb=50000,
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.rscript} create_seurat_from_h5ad.R --h5ad ../data/processed/{params.h5ad_rel} \
            > {params.project_dir}/{log} 2>&1
        """


# Step 3b — Seurat label transfer from merged scRNA reference

rule label_transfer:
    input:
        rds=RDS_OUT,
        ref_rds=REF_RDS,
    output:
        transfer=TRANSFER_OUT,
    params:
        rds_abs=lambda wildcards: str(
            Path(config["project_dir"]) / RDS_OUT.format(cohort=wildcards.cohort, label=wildcards.label)
        ),
        out_abs=lambda wildcards: str(
            Path(config["project_dir"]) / TRANSFER_OUT.format(cohort=wildcards.cohort, label=wildcards.label)
        ),
        ref_rds=str(Path(config["ref_out"]) / config["ref_rds"]),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        rscript=RSCRIPT,
    log:
        f"logs/label_transfer_{{cohort}}_{{label}}_{REF_TAG}.log",
    wildcard_constraints:
        cohort=COHORTS,
        label="|".join(unfiltered_transfer_labels()) if unfiltered_transfer_labels() else "NEVER_MATCH",
    resources:
        mem_mb=50000,
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.rscript} label_transfer.R \
            --project_dir {params.project_dir} \
            --vis_file {params.rds_abs} \
            --ref_rds {params.ref_rds} \
            --out_file {params.out_abs} \
            > {params.project_dir}/{log} 2>&1
        """


rule label_transfer_nontumor:
    input:
        rds=RDS_NONTUMOR,
        ref_rds=REF_RDS,
    output:
        transfer=TRANSFER_NONTUMOR_OUT,
    params:
        rds_abs=lambda wildcards: str(
            Path(config["project_dir"]) / RDS_NONTUMOR.format(cohort=wildcards.cohort, label=wildcards.label)
        ),
        out_abs=lambda wildcards: str(
            Path(config["project_dir"]) / TRANSFER_NONTUMOR_OUT.format(cohort=wildcards.cohort, label=wildcards.label)
        ),
        ref_rds=str(Path(config["ref_out"]) / config["ref_rds"]),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        rscript=RSCRIPT,
    log:
        f"logs/label_transfer_{{cohort}}_{{label}}_{REF_TAG}_nontumor.log",
    wildcard_constraints:
        cohort=COHORTS,
        label="|".join(filtered_transfer_labels()) if filtered_transfer_labels() else "NEVER_MATCH",
    resources:
        mem_mb=50000,
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.rscript} label_transfer.R \
            --project_dir {params.project_dir} \
            --vis_file {params.rds_abs} \
            --ref_rds {params.ref_rds} \
            --out_file {params.out_abs} \
            > {params.project_dir}/{log} 2>&1
        """
