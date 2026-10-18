# Step 4 — Cellestimate on Spatial assay (requires Seurat rds)

rule run_cellestimate:
    input:
        rds=RDS_OUT,
    output:
        celltype=CELLESTIMATE_OUT,
    params:
        script=config["cellestimate_script"],
        db=config["cellestimate_db"],
        cf=CELLESTIMATE_CF_ARG,
        sheet=config.get("cellestimate_sheet", "subsclass_immune"),
        one_class_split=config.get("cellestimate_one_class_split", "TRUE"),
        rds_abs=lambda wildcards: str(
            Path(config["project_dir"]) / RDS_OUT.format(cohort=wildcards.cohort, label=wildcards.label)
        ),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        rscript=RSCRIPT,
    log:
        f"logs/cellestimate_{{cohort}}_{{label}}_{REF_TAG}.log",
    wildcard_constraints:
        cohort=COHORTS,
        label="|".join(unfiltered_cellestimate_labels()) if unfiltered_cellestimate_labels() else "NEVER_MATCH",
    resources:
        mem_mb=50000,
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        {params.rscript} {params.script} \
            --rds {params.rds_abs} \
            --db {params.db} \
            --assay Spatial \
            --cf {params.cf} \
            --sheet {params.sheet} \
            --one_class_split {params.one_class_split} \
            > {params.project_dir}/{log} 2>&1
        """


rule run_cellestimate_nontumor:
    input:
        rds=RDS_NONTUMOR,
    output:
        celltype=CELLESTIMATE_NONTUMOR_OUT,
    params:
        script=config["cellestimate_script"],
        db=config["cellestimate_db"],
        cf=CELLESTIMATE_CF_ARG,
        sheet=config.get("cellestimate_sheet", "subsclass_immune"),
        one_class_split=config.get("cellestimate_one_class_split", "TRUE"),
        rds_abs=lambda wildcards: str(
            Path(config["project_dir"]) / RDS_NONTUMOR.format(cohort=wildcards.cohort, label=wildcards.label)
        ),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        rscript=RSCRIPT,
    log:
        f"logs/cellestimate_{{cohort}}_{{label}}_{REF_TAG}_nontumor.log",
    wildcard_constraints:
        cohort=COHORTS,
        label="|".join(filtered_cellestimate_labels()) if filtered_cellestimate_labels() else "NEVER_MATCH",
    resources:
        mem_mb=50000,
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        {params.rscript} {params.script} \
            --rds {params.rds_abs} \
            --db {params.db} \
            --assay Spatial \
            --cf {params.cf} \
            --sheet {params.sheet} \
            --one_class_split {params.one_class_split} \
            > {params.project_dir}/{log} 2>&1
        """
