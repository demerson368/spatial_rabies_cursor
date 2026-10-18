bsub -o plot_score_celltype_html.o -e plot_score_celltype_html.e -R "rusage[mem=30000]" python scripts/plot_score_celltype_html.py \
    --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um \
    --gene EGFP \
    --loupe_csv for_loupe/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_for_loupe.csv \
    --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
    --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
    --highlight_celltype Glut  \
    --if_downsample 0.25  \
    --score_vmin 0.1  \
    --highlight_opacity 0.3  \
    --if_channel rgb

bsub -o plot_score_celltype_html.o -e plot_score_celltype_html.e -R "rusage[mem=30000]" python scripts/plot_score_celltype_html.py \
    --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um \
    --gene EGFP \
    --loupe_csv for_loupe/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_for_loupe.csv \
    --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
    --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
    --highlight_celltype GABA  \
    --if_downsample 0.25  \
    --score_vmin 0.1  \
    --highlight_opacity 0.3  \
    --if_channel rgb

bsub -o plot_score_celltype_html.o -e plot_score_celltype_html.e -R "rusage[mem=30000]"  python scripts/plot_score_celltype_html.py \
    --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um \
    --gene EGFP \
    --loupe_csv for_loupe/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_for_loupe.csv \
    --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
    --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
    --highlight_celltype "Microglia NN"  \
    --if_downsample 0.25  \
    --score_vmin 0.1  \
    --highlight_opacity 0.3  \
    --if_channel rgb





#bsub -o plot_score_celltype_html.o -e plot_score_celltype_html.e -R "rusage[mem=30000]" python scripts/plot_score_celltype_html.py \
#    --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um \
#    --signature zr_fusion \
#  --loupe_csv for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_for_loupe.csv \
#  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
#  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
#  --highlight_celltype Glut  \
#  --if_downsample 0.25  \
#  --score_vmin 0.1  \
#  --highlight_opacity 0.3  \
#  --if_channel rgb
#
#    --celltype_col 
