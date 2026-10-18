#bsub -o plot_gene_score_celltype_distribution.o -e plot_gene_score_celltype_distribution.e -R "rusage[mem=50000]" python scripts/plot_gene_score_celltype_distribution.py \
#  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
#  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
#  --gene Grin2d \
#  --celltype "Microglia NN" \
#  --compare_all \
#  --label pooled_D90_5_A1_D90_7_D1


#bsub -o plot_gene_score_celltype_distribution.o -e plot_gene_score_celltype_distribution.e -R "rusage[mem=50000]" python scripts/plot_gene_score_celltype_distribution.py \
#  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
#  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
#  --gene Grin2d \
#  --all_celltypes \
#  --label pooled_D90_5_A1_D90_7_D1




bsub -o plot_gene_score_celltype_distribution.o -e plot_gene_score_celltype_distribution.e -R "rusage[mem=50000]" python scripts/plot_gene_score_celltype_distribution.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --gene Grin2d \
  --celltype_col cest_w_tumor \
  --all_celltypes \
  --include_celltype GABA \
  --include_celltype Glut \
  --include_celltype "Microglia NN" \
  --include_celltype tumor \
  --label pooled_D90_5_A1_D90_7_D1 \
  --score_log log1p


bsub -o plot_gene_score_celltype_distribution.o -e plot_gene_score_celltype_distribution.e -R "rusage[mem=50000]" python scripts/plot_gene_score_celltype_distribution.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --gene Chrm3 \
  --celltype_col cest_w_tumor \
  --all_celltypes \
  --include_celltype GABA \
  --include_celltype Glut \
  --include_celltype "Microglia NN" \
  --include_celltype tumor \
  --label pooled_D90_5_A1_D90_7_D1 \
  --score_log log1p

bsub -o plot_gene_score_celltype_distribution.o -e plot_gene_score_celltype_distribution.e -R "rusage[mem=50000]" python scripts/plot_gene_score_celltype_distribution.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --gene Grm4 \
  --celltype_col cest_w_tumor \
  --all_celltypes \
  --include_celltype GABA \
  --include_celltype Glut \
  --include_celltype "Microglia NN" \
  --include_celltype tumor \
  --label pooled_D90_5_A1_D90_7_D1 \
  --score_log log1p

bsub -o plot_gene_score_celltype_distribution.o -e plot_gene_score_celltype_distribution.e -R "rusage[mem=50000]" python scripts/plot_gene_score_celltype_distribution.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --gene Nrxn1 \
  --celltype_col cest_w_tumor \
  --all_celltypes \
  --include_celltype GABA \
  --include_celltype Glut \
  --include_celltype "Microglia NN" \
  --include_celltype tumor \
  --label pooled_D90_5_A1_D90_7_D1 \
  --score_log log1p

bsub -o plot_gene_score_celltype_distribution.o -e plot_gene_score_celltype_distribution.e -R "rusage[mem=50000]" python scripts/plot_gene_score_celltype_distribution.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --gene Gabra3 \
  --celltype_col cest_w_tumor \
  --all_celltypes \
  --include_celltype GABA \
  --include_celltype Glut \
  --include_celltype "Microglia NN" \
  --include_celltype tumor \
  --label pooled_D90_5_A1_D90_7_D1 \
  --score_log log1p


bsub -o plot_gene_score_celltype_distribution.o -e plot_gene_score_celltype_distribution.e -R "rusage[mem=50000]" python scripts/plot_gene_score_celltype_distribution.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --gene nrxn2 \
  --celltype_col cest_w_tumor \
  --all_celltypes \
  --include_celltype GABA \
  --include_celltype Glut \
  --include_celltype "Microglia NN" \
  --include_celltype tumor \
  --label pooled_D90_5_A1_D90_7_D1 \


bsub -o plot_gene_score_celltype_distribution.o -e plot_gene_score_celltype_distribution.e -R "rusage[mem=50000]" python scripts/plot_gene_score_celltype_distribution.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --gene Chrnb2 \
  --celltype_col cest_w_tumor \
  --all_celltypes \
  --include_celltype GABA \
  --include_celltype Glut \
  --include_celltype "Microglia NN" \
  --include_celltype tumor \
  --label pooled_D90_5_A1_D90_7_D1 \

bsub -o plot_gene_score_celltype_distribution.o -e plot_gene_score_celltype_distribution.e -R "rusage[mem=50000]" python scripts/plot_gene_score_celltype_distribution.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment_cest_w_tumor.h5ad \
  --gene Chrna4 \
  --celltype_col cest_w_tumor \
  --all_celltypes \
  --include_celltype GABA \
  --include_celltype Glut \
  --include_celltype "Microglia NN" \
  --include_celltype tumor \
  --label pooled_D90_5_A1_D90_7_D1 \
