#conda activate spatial

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Load directly from SpaceRanger output folder

path = "/research/groups/mackgrp/home/common/Collaborations/RABV/3356629_D49_D1_SteveM_Tuyu_GFP_v2/3356629_D49_D1_SteveM_Tuyu_GFP_v2/outs/binned_outputs/square_016um"

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
plt.savefig("spatial_total_counts_test.png",dpi=300,bbox_inches="tight")
