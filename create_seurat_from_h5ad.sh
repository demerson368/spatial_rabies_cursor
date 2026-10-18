#bsub -o create_seurat_from_h5ad.o -e create_seurat_from_h5ad.e -R "rusage[mem=50000]" Rscript scripts/create_seurat_from_h5ad.R --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um.h5ad

#bsub -o create_seurat_from_h5ad.o -e create_seurat_from_h5ad.e -R "rusage[mem=50000]" Rscript scripts/create_seurat_from_h5ad.R --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um.h5ad

#bsub -o create_seurat_from_h5ad.o -e create_seurat_from_h5ad.e -R "rusage[mem=50000]" Rscript scripts/create_seurat_from_h5ad.R --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J74_3_A1_GFP_R1_008um.h5ad

#bsub -o create_seurat_from_h5ad.o -e create_seurat_from_h5ad.e -R "rusage[mem=50000]" Rscript scripts/create_seurat_from_h5ad.R --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J75_3_D1_GFP_R1_008um.h5ad


#bsub -o create_seurat_from_h5ad.o -e create_seurat_from_h5ad.e -R "rusage[mem=50000]" Rscript scripts/create_seurat_from_h5ad.R --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad

#bsub -o create_seurat_from_h5ad.o -e create_seurat_from_h5ad.e -R "rusage[mem=50000]" Rscript scripts/create_seurat_from_h5ad.R --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad

bsub -o create_seurat_from_h5ad.o -e create_seurat_from_h5ad.e -R "rusage[mem=50000]" Rscript scripts/create_seurat_from_h5ad.R --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad

bsub -o create_seurat_from_h5ad.o -e create_seurat_from_h5ad.e -R "rusage[mem=50000]" Rscript scripts/create_seurat_from_h5ad.R --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad
