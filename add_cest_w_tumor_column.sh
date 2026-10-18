python scripts/add_cest_w_tumor_column.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --tumor_threshold 0.1

python scripts/add_cest_w_tumor_column.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --tumor_threshold 0.1
