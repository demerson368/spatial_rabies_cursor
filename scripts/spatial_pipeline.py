#conda activate spatial

import os
import argparse
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Load directly from SpaceRanger output folder
def main():

    output_dir = '../data/processed/'

    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str, help="folder to space ranger output")
    parser.add_argument('label', type=str, help="name of spatial output")
    parser.add_argument('--cohort', type=str, required=True, help="cohort subdir under data/processed/")
    args = parser.parse_args()

    path = args.path
    label = args.label
    full_output = os.path.join(output_dir, args.cohort)
    os.makedirs(full_output, exist_ok=True)

    adata_st = sc.read_visium(
        path=path,          # folder containing spatial/ and matrix/
        count_file="filtered_feature_bc_matrix.h5",
        load_images=True,
    )

    #squidy recommended?
    #adata_st = sq.read.visium(path,counts_file="filtered_feature_bc_matrix.h5")

    adata_st.var_names_make_unique()

    # Confirm spatial coordinates loaded
    print(adata_st.obsm["spatial"])          # (n_spots, 2) pixel coordinates
    print(adata_st.uns["spatial"])           # image + scalefactors

    # Quick tissue plot
    sc.pp.calculate_qc_metrics(adata_st, inplace=True)
    sq.pl.spatial_scatter(adata_st, color="total_counts", size=1.5)
    plt.savefig(
        os.path.join(full_output, f"spatial_total_counts_{label}.png"),
        dpi=300, bbox_inches="tight"
    )
    plt.close()

    #Step 2 — QC and normalization
    sc.pp.calculate_qc_metrics(adata_st, inplace=True)

    # Spatial QC — plot metrics on tissue
    sq.pl.spatial_scatter(adata_st, color=["total_counts", "n_genes_by_counts"], size=1.5)
    plt.savefig(
        os.path.join(full_output, f"spatial_qc_{label}.png"),
        dpi=300, bbox_inches="tight"
    )
    plt.close()


    # Filter low-quality spots
    sc.pp.filter_cells(adata_st, min_counts=500)
    sc.pp.filter_genes(adata_st, min_cells=10)

    # Stash raw counts
    adata_st.layers["counts"] = adata_st.X.copy()

    sc.pp.normalize_total(adata_st, target_sum=1e4)
    sc.pp.log1p(adata_st)


    # Spatial neighbors (used by downstream spatial stats if needed)
    sq.gr.spatial_neighbors(adata_st, coord_type="grid", n_neighs=6)

    # Moran's I / SVG plots are a separate optional Snakemake rule:
    #   scripts/run_morans_i.py  →  rule run_morans_i / all_morans
    # Keep them out of this critical-path preprocess step.

    adata_st.write_h5ad(os.path.join(full_output, f"spatial_{label}.h5ad"))
    print(f"Saved: spatial_{label}.h5ad")


if __name__ == '__main__':
    main()
    print("Script complete — exiting")
    import sys
    sys.exit(0)   # force clean exit
