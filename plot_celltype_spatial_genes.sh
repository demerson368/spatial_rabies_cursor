#bsub -o plot_celltype_spatial_genes.o -e plot_celltype_spatial_genes.e  -R "rusage[mem=40000]" python scripts/plot_celltype_spatial_genes.py --loupe_csv for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_for_loupe.csv --h5ad  data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_ext.h5ad



# Default signature file: reference/signature_genes_by_cellName_wide.csv
# For merged broad objects use --celltype_value "GABA" (matches cellName column).
# Legacy P GABA Yao Allen list: --gene_set_file reference/p_gaba_biomarkers.txt --celltype_value "P GABA"
#
#bsub -o plot_celltype_spatial_genes.o -e plot_celltype_spatial_genes.e  -R "rusage[mem=40000]" python scripts/plot_celltype_spatial_genes.py --loupe_csv for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_for_loupe.csv --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_ext.h5ad --gene_set_file reference/p_gaba_biomarkers.txt --celltype_value "P GABA" --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um --top_n 12


#bsub -o plot_celltype_spatial_genes.o -e plot_celltype_spatial_genes.e  -R "rusage[mem=40000]" python scripts/plot_celltype_spatial_genes.py --loupe_csv for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_for_loupe.csv --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext.h5ad --gene_set_file reference/signature_genes_by_cellName_wide.csv --celltype_value GABA --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad --top_n 20

#bsub -o plot_celltype_spatial_genes.o -e plot_celltype_spatial_genes.e  -R "rusage[mem=40000]" python scripts/plot_celltype_spatial_genes.py --loupe_csv for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_for_loupe.csv --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext.h5ad --gene_set_file reference/signature_genes_by_cellName_wide.csv --celltype_col cellestimate --celltype_value GABA --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad --top_n 20

#bsub -o plot_celltype_spatial_genes.o -e plot_celltype_spatial_genes.e  -R "rusage[mem=40000]" python scripts/plot_celltype_spatial_genes.py --loupe_csv for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_for_loupe.csv --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext.h5ad --gene_set_file reference/signature_genes_by_cellName_wide.csv --celltype_value Glut --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad --top_n 20

bsub -o plot_celltype_spatial_genes.o -e plot_celltype_spatial_genes.e  -R "rusage[mem=40000]" python scripts/plot_celltype_spatial_genes.py --loupe_csv for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_for_loupe.csv --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext.h5ad --gene_set_file reference/signature_genes_by_cellName_wide.csv --celltype_col cellestimate --celltype_value Glut --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad --top_n 20

