#Note using cell_type from integrated*_mito_doublet_removal_SCT_update (corresponds to cell_type_post from the raw where each batch is kept separate (integrated_ref_rawcounts_with_labels.rds))
#assignment to cell_type after integration of samples and renormalization and then cell_type biomarker calling (cell_type_post))



#create rabies info from spatial data
python scripts/count_rabies_gfp_cells.py \
    --h5ad "/research/groups/mackgrp/home/demerson/spatial_rabies_cursor/data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext.h5ad" \
    --no-batch --cell-type-col cellestimate --annotate-h5ad


python scripts/count_rabies_gfp_cells.py \
    --h5ad "/research/groups/mackgrp/home/demerson/spatial_rabies_cursor/data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext.h5ad" \
     --no-batch --cell-type-col cellestimate --annotate-h5ad
