#Note: use scripts/celltype_biomarker_enrichment.py for any cell2location label.


my_array=("Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um" "Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um" "Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J74_3_A1_GFP_R1_008um" "Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J75_3_D1_GFP_R1_008um")

for label in "${my_array[@]}"; do
    echo "$label"
    bsub -o p_gaba_biomarker_enrichment.o -e p_gaba_biomarker_enrichment.e -R "rusage[mem=40000]" -q gpu -gpu "num=1/host" python scripts/celltype_biomarker_enrichment.py \
    --loupe_csv for_loupe/spatial_${label}_for_loupe.csv \
    --h5ad data/processed/Tuyu_Apr2026/spatial_${label}_ext.h5ad \
    --celltype_value "MB Dopa" \
    --gene_set_file reference/mb_dopa_genes_allen_yao_et_al.txt \
    --label "$label"
done





