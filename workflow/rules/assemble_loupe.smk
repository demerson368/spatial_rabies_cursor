# Step 5 — tlscore / NNLS / integrate all methods → Loupe CSV + extended h5ad

def assemble_inputs(wildcards):
    label = wildcards.label
    row = SAMPLES.loc[label]
    inputs = {
        "h5ad": H5AD_OUT.format(label=label),
    }
    if row.run_cellestimate:
        inputs["cellestimate"] = CELLESTIMATE_OUT.format(label=label)
    if row.run_seurat_transfer:
        inputs["transfer"] = TRANSFER_OUT.format(label=label)
    if row.run_cell2location:
        inputs["cell2location"] = C2L_OUT.format(label=label)
    return inputs


rule assemble_loupe:
    input:
        unpack(assemble_inputs),
        marker_excel=config["marker_excel"],
    output:
        loupe=LOUPE_OUT,
        h5ad_ext=H5AD_EXT,
    params:
        label=lambda wildcards: wildcards.label,
        ref_tag=REF_TAG,
        cohort=COHORT,
        cellestimate_cf=config["cellestimate_cf"],
        marker_excel=config["marker_excel"],
        marker_sheet=config["marker_sheet"],
        marker_cluster_col=config["marker_cluster_col"],
        marker_class_col=config["marker_class_col"],
        h5ad_abs=lambda wildcards: str(
            Path(config["project_dir"]) / H5AD_OUT.format(label=wildcards.label)
        ),
        h5ad_ext_abs=lambda wildcards: str(
            Path(config["project_dir"]) / H5AD_EXT.format(label=wildcards.label)
        ),
        loupe_abs=lambda wildcards: str(
            Path(config["project_dir"]) / LOUPE_OUT.format(label=wildcards.label)
        ),
        c2l_abs=lambda wildcards: str(
            Path(config["project_dir"]) / C2L_OUT.format(label=wildcards.label)
        ),
        transfer_abs=lambda wildcards: str(
            Path(config["project_dir"]) / TRANSFER_OUT.format(label=wildcards.label)
        ),
        cellestimate_abs=lambda wildcards: str(
            Path(config["project_dir"]) / CELLESTIMATE_OUT.format(label=wildcards.label)
        ),
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
    log:
        f"logs/assemble_loupe_{{label}}_{REF_TAG}.log",
    wildcard_constraints:
        label="|".join(assemble_labels()) if assemble_labels() else "NEVER_MATCH",
    resources:
        mem_mb=50000,
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.python} assemble_loupe_file.py \
            --h5ad {params.h5ad_abs} \
            --label {params.label} \
            --ref_tag {params.ref_tag} \
            --cohort {params.cohort} \
            --cellestimate_cf {params.cellestimate_cf} \
            --h5ad_ext {params.h5ad_ext_abs} \
            --loupe_out {params.loupe_abs} \
            --c2l_h5ad {params.c2l_abs} \
            --transfer_txt {params.transfer_abs} \
            --cellestimate_txt {params.cellestimate_abs} \
            --marker_excel {params.marker_excel} \
            --marker_sheet {params.marker_sheet} \
            --marker_cluster_col {params.marker_cluster_col} \
            --marker_class_col {params.marker_class_col} \
            > {params.project_dir}/{log} 2>&1
        """
