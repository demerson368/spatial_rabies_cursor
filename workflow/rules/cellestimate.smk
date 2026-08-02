# Step 4 — Cellestimate on Spatial assay (requires Seurat rds)

rule run_cellestimate:
    input:
        rds=RDS_OUT,
    output:
        celltype=CELLESTIMATE_OUT,
    params:
        script=config["cellestimate_script"],
        db=config["cellestimate_db"],
        cf=config["cellestimate_cf"],
        rds_abs=lambda wildcards: str(
            Path(config["project_dir"]) / RDS_OUT.format(label=wildcards.label)
        ),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        rscript=RSCRIPT,
    log:
        "logs/cellestimate_{label}.log",
    wildcard_constraints:
        label="|".join(cellestimate_labels()) if cellestimate_labels() else "NEVER_MATCH",
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
            > {params.project_dir}/{log} 2>&1
        """
