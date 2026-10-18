#Note: use scripts/celltype_biomarker_enrichment.py for any cell2location label.
# P GABA wrapper (same defaults): scripts/p_gaba_biomarker_enrichment.py
#
# Example — MB Dopa:
# python scripts/celltype_biomarker_enrichment.py \
#   --loupe_csv for_loupe/spatial_{label}_for_loupe.csv \
#   --h5ad data/processed/Tuyu_Apr2026/spatial_{label}_ext.h5ad \
#   --celltype_value "MB Dopa" \
#   --gene_set_file reference/mb_dopa_genes_allen_yao_et_al.txt \
#   --label {label}
#
# Output stem: spatial_{label}_{celltype_slug}_*  (e.g. ..._mb_dopa_gsea.csv)
#
# P GABA signature-gene spatial ranking (default: reference/p_gaba_biomarkers.txt):
# python scripts/plot_celltype_spatial_genes.py \
#   --loupe_csv for_loupe/spatial_{label}_for_loupe.csv \
#   --h5ad data/processed/Tuyu_Apr2026/spatial_{label}_ext.h5ad \
#   --gene_set_file reference/p_gaba_biomarkers.txt \
#   --celltype_value "P GABA" --label {label}

#spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_ext.h5ad same as spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um.h5ad except added spot celltype annotations
#bsub -o p_gaba_biomarker_enrichment.o -e p_gaba_biomarker_enrichment.e -R "rusage[mem=20000]" -q gpu -gpu "num=1/host" python scripts/p_gaba_biomarker_enrichment.py --loupe_csv for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_for_loupe.csv  --gene_set_file reference/p_gaba_genes_allen_yao_et_al.txt --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_ext.h5ad --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um

#spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_ext.h5ad same as spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um.h5ad except added spot celltype annotations
bsub -o p_gaba_biomarker_enrichment.o -e p_gaba_biomarker_enrichment.e -R "rusage[mem=40000]" -q gpu -gpu "num=1/host" python scripts/p_gaba_biomarker_enrichment.py --loupe_csv for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_for_loupe.csv  --gene_set_file reference/p_gaba_genes_allen_yao_et_al.txt --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_ext.h5ad --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um


#bsub -o p_gaba_biomarker_enrichment.o -e p_gaba_biomarker_enrichment.e -R "rusage[mem=30000]" -q gpu -gpu "num=1/host" python scripts/p_gaba_biomarker_enrichment.py --loupe_csv for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J74_3_A1_GFP_R1_008um_for_loupe.csv  --gene_set_file reference/p_gaba_genes_allen_yao_et_al.txt --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J74_3_A1_GFP_R1_008um_ext.h5ad --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J74_3_A1_GFP_R1_008um


#bsub -o p_gaba_biomarker_enrichment.o -e p_gaba_biomarker_enrichment.e -R "rusage[mem=30000]" -q gpu -gpu "num=1/host" python scripts/p_gaba_biomarker_enrichment.py --loupe_csv for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J75_3_D1_GFP_R1_008um_for_loupe.csv  --gene_set_file reference/p_gaba_genes_allen_yao_et_al.txt --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J75_3_D1_GFP_R1_008um_ext.h5ad --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J75_3_D1_GFP_R1_008um




