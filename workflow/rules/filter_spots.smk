# Optional — score ZR fusion and subset to non-tumor spots after preprocess

rule filter_nontumor_spots:
    input:
        h5ad=H5AD_OUT,
    output:
        h5ad_nontumor=H5AD_NONTUMOR,
        filter_report=FILTER_REPORT,
        filter_qc=FILTER_QC,
    params:
        max_score=lambda wc: zr_fusion_max_score(SAMPLES, wc.label, config),
        h5ad_in=lambda wc: str(Path(config["project_dir"]) / H5AD_OUT.format(cohort=wc.cohort, label=wc.label)),
        h5ad_out=lambda wc: str(Path(config["project_dir"]) / H5AD_NONTUMOR.format(cohort=wc.cohort, label=wc.label)),
        report_out=lambda wc: str(Path(config["project_dir"]) / FILTER_REPORT.format(cohort=wc.cohort, label=wc.label)),
        qc_out=lambda wc: str(Path(config["project_dir"]) / FILTER_QC.format(cohort=wc.cohort, label=wc.label)),
        cohort=lambda wc: wc.cohort,
        label=lambda wildcards: wildcards.label,
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
    log:
        "logs/filter_nontumor_{cohort}_{label}.log",
    wildcard_constraints:
        cohort=COHORTS,
        label="|".join(nontumor_filter_labels(SAMPLES)) if nontumor_filter_labels(SAMPLES) else "NEVER_MATCH",
    resources:
        mem_mb=50000,
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.python} filter_nontumor_spots.py \
            --h5ad {params.h5ad_in} \
            --label {params.label} \
            --cohort {params.cohort} \
            --max_score {params.max_score} \
            --out_h5ad {params.h5ad_out} \
            --report_csv {params.report_out} \
            --qc_plot {params.qc_out} \
            > {params.project_dir}/{log} 2>&1
        """
