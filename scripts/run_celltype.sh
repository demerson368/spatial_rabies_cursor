#/research/groups/mackgrp/home/demerson/CellTypeEstimate/cte4so.R --rds ../data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um.rds --db Yao_et_al_mouse_format_addendum --assay Spatial --cf Yao_et_al_mouse_format

bsub -o run_celltype.o -e run_celltype.e -R "rusage[mem=50000]" /research/groups/mackgrp/home/demerson/CellTypeEstimate/cte4so.R --rds ../data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um.rds --db Yao_et_al_mouse_format_addendum --assay Spatial --cf Yao_et_al_mouse_format

#bsub -o run_celltype.o -e run_celltype.e -R "rusage[mem=50000]" /research/groups/mackgrp/home/demerson/CellTypeEstimate/cte4so.R --rds ../data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J75_3_D1_GFP_R1_008um.rds  --db Yao_et_al_mouse_format_addendum --assay Spatial --cf Yao_et_al_mouse_format

#bsub -o run_celltype.o -e run_celltype.e -R "rusage[mem=50000]" /research/groups/mackgrp/home/demerson/CellTypeEstimate/cte4so.R --rds ../data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J74_3_A1_GFP_R1_008um.rds  --db Yao_et_al_mouse_format_addendum --assay Spatial --cf Yao_et_al_mouse_format
