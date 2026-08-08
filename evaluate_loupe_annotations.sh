#bsub -o evaluate_loupe_annotations.o -e evaluate_loupe_annotations.e -R "rusage[mem=30000]" python scripts/evaluate_loupe_annotations.py --input_dir for_loupe --samples Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um --label_cols cell2location --n_perms 10  --skip_autocorr --skip_enrichment

#bsub -o evaluate_loupe_annotations.o -e evaluate_loupe_annotations.e -R "rusage[mem=30000]" python scripts/evaluate_loupe_annotations.py --input_dir for_loupe --samples Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um --n_perms 10  --skip_autocorr --skip_enrichment

#bsub -o evaluate_loupe_annotations.o -e evaluate_loupe_annotations.e -R "rusage[mem=30000]" python scripts/evaluate_loupe_annotations.py --input_dir for_loupe --samples Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um  --n_perms 5 --skip_autocorr

#bsub -o evaluate_loupe_annotations.o -e evaluate_loupe_annotations.e -R "rusage[mem=30000]" python scripts/evaluate_loupe_annotations.py --input_dir for_loupe --samples Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um  --n_perms 5 --skip_enrichment

#bsub -o evaluate_loupe_annotations.o -e evaluate_loupe_annotations.e -R "rusage[mem=10000]" python scripts/evaluate_loupe_annotations.py --input_dir for_loupe --samples Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um  --n_perms 5 --skip_enrichment

bsub -o evaluate_loupe_annotations.o -e evaluate_loupe_annotations.e -R "rusage[mem=10000]" python scripts/evaluate_loupe_annotations.py --input_dir for_loupe --samples Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J74_3_A1_GFP_R1_008um  --n_perms 5 --skip_enrichment

bsub -o evaluate_loupe_annotations.o -e evaluate_loupe_annotations.e -R "rusage[mem=10000]" python scripts/evaluate_loupe_annotations.py --input_dir for_loupe --samples Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J75_3_D1_GFP_R1_008um  --n_perms 5 --skip_enrichment

