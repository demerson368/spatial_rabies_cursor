bsub -o assemble_loupe_file.o -e assemble_loupe_file.e -R "rusage[mem=50000]" python assemble_loupe_file.py --h5ad ../data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um.h5ad --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um

#bsub -o assemble_loupe_file.o -e assemble_loupe_file.e -R "rusage[mem=50000]" python assemble_loupe_file.py --h5ad ../data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J74_3_A1_GFP_R1_008um.h5ad --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J74_3_A1_GFP_R1_008um

#bsub -o assemble_loupe_file.o -e assemble_loupe_file.e -R "rusage[mem=50000]" python assemble_loupe_file.py --h5ad ../data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J75_3_D1_GFP_R1_008um.h5ad --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J75_3_D1_GFP_R1_008um 
