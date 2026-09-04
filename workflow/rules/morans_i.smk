# Optional — Moran's I / spatially variable genes (slow on Visium HD)
# Does NOT block cell2location, Seurat, or assemble_loupe.

rule run_morans_i:
    input:
        h5ad=H5AD_OUT,
    output:
        moran_csv=MORAN_CSV,
        svg_plot=SVG_PLOT,
    params:
        out_dir=PROCESSED,
        n_perms=config.get("morans_n_perms", 100),
        n_jobs=config.get("morans_n_jobs", 4),
        n_top=config.get("morans_n_top", 6),
        n_genes=config.get("morans_n_genes", 3000),
        ref_suffix=REF_SUFFIX,
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
    log:
        f"logs/morans_i_{{label}}_{REF_TAG}.log",
    wildcard_constraints:
        label="|".join(unfiltered_morans_labels()) if unfiltered_morans_labels() else "NEVER_MATCH",
    threads: 4
    resources:
        mem_mb=100000,
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.python} run_morans_i.py \
            --h5ad {params.project_dir}/{input.h5ad} \
            --label {wildcards.label} \
            --out_dir {params.project_dir}/{params.out_dir} \
            --n_perms {params.n_perms} \
            --n_jobs {params.n_jobs} \
            --n_top {params.n_top} \
            --n_genes {params.n_genes} \
            --ref_suffix {params.ref_suffix} \
            > {params.project_dir}/{log} 2>&1
        """


rule run_morans_i_nontumor:
    input:
        h5ad=H5AD_NONTUMOR,
    output:
        moran_csv=MORAN_CSV_NONTUMOR,
        svg_plot=SVG_PLOT_NONTUMOR,
    params:
        out_dir=PROCESSED,
        n_perms=config.get("morans_n_perms", 100),
        n_jobs=config.get("morans_n_jobs", 4),
        n_top=config.get("morans_n_top", 6),
        n_genes=config.get("morans_n_genes", 3000),
        ref_suffix=f"{REF_SUFFIX}_nontumor",
        scripts_dir=str(Path(config["project_dir"]) / "scripts"),
        project_dir=config["project_dir"],
        conda_prefix=CONDA_PREFIX,
        python=PYTHON,
    log:
        f"logs/morans_i_{{label}}_{REF_TAG}_nontumor.log",
    wildcard_constraints:
        label="|".join(filtered_morans_labels()) if filtered_morans_labels() else "NEVER_MATCH",
    threads: 4
    resources:
        mem_mb=100000,
    shell:
        """
        set -eu
        export LD_LIBRARY_PATH="{params.conda_prefix}/lib:${{LD_LIBRARY_PATH:-}}"
        mkdir -p {params.project_dir}/logs
        cd {params.scripts_dir}
        {params.python} run_morans_i.py \
            --h5ad {params.project_dir}/{input.h5ad} \
            --label {wildcards.label} \
            --out_dir {params.project_dir}/{params.out_dir} \
            --n_perms {params.n_perms} \
            --n_jobs {params.n_jobs} \
            --n_top {params.n_top} \
            --n_genes {params.n_genes} \
            --ref_suffix {params.ref_suffix} \
            > {params.project_dir}/{log} 2>&1
        """


rule all_morans:
    input:
        expand(MORAN_CSV, label=unfiltered_morans_labels()),
        expand(SVG_PLOT, label=unfiltered_morans_labels()),
        expand(MORAN_CSV_NONTUMOR, label=filtered_morans_labels()),
        expand(SVG_PLOT_NONTUMOR, label=filtered_morans_labels()),
