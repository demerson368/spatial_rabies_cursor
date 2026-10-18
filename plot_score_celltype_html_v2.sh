#bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=30000]" python scripts/plot_score_celltype_html_v2.py \
#  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
#  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
#  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
#  --gene Ptprc \
#  --celltype_col rabies_assignment \
#  --celltype_col cellestimate \
#  --highlight_celltype positive \
#  --highlight_celltype Glut


#python scripts/plot_score_celltype_html_v2.py \
#  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
#  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
#  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
#  --no_score \
#  --celltype_col cellestimate \
#  --celltype_col rabies_assignment \
#  --highlight_celltype Glut \
#  --highlight_celltype positive \
#  --highlight_opacity 0.2 \
#  --highlight_opacity 0.3




#python scripts/plot_score_celltype_html_v2.py \
#  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
#  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
#  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-7-D1-GFP_R1_IF_transform_reoriented.tif  \
#  --no_score \
#  --celltype_col cellestimate \
#  --celltype_col rabies_assignment \
#  --highlight_celltype Glut \
#  --highlight_celltype positive \
#  --highlight_opacity 0.2 \
#  --highlight_opacity 0.3



bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
  --signature nrxn_signature \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5  \
  --spot_size 3.0
#  --score_log log1p
#  --score_vmin 2  \


bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
  --signature nlgn_signature \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_log log1p \
  --spot_size 3.0
##  --score_vmin 0.1  \


bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
  --signature signatures_set2 \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5  \
  --spot_size 3.0



bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
  --signature signatures_synaptic \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5 \
  --spot_size 3.0 

bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Grin2d \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5 \
  --spot_size 3.0

bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Nrxn1 \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5 \
  --spot_size 3.0


bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Grin2d \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --extra_score_col tltumscore_C11orf95_relafus1 \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5  \
  --extra_score_vmin 0.1 \
  --spot_size 3.0



bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Grm4 \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5 \
  --spot_size 3.0

bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Gabra3 \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5 \
  --spot_size 3.0



bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Chrnb2 \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5 \
  --spot_size 3.0 

bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Chrna4 \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5 \
  --spot_size 3.0

#bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
#  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
#  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
#  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-7-D1-GFP_R1_IF_transform_reoriented.tif  \
#  --gene Grin2d \
#  --celltype_col cellestimate \
#  --celltype_col rabies_assignment \
#  --highlight_celltype Glut \
#  --highlight_celltype positive \
#  --highlight_opacity 0.3 \
#  --highlight_opacity 1.0 \
#  --score_vmin 1.5  \
#  --spot_size 3.0

bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-7-D1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Grm4 \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5  \
  --spot_size 3.0

bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]"  python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-7-D1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Gabra3 \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5  \
  --spot_size 3.0

bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]"  python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-7-D1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Nrxn1 \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5  \
  --spot_size 3.0



bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-7-D1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Grin2d \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --extra_score_col tltumscore_C11orf95_relafus1 \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5  \
  --extra_score_vmin 0.1 \
  --spot_size 3.0



bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-7-D1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Chrnb2 \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5  \
  --spot_size 3.0

bsub -o plot_score_celltype_html_v2.o -e plot_score_celltype_html_v2.e -R "rusage[mem=50000]" python scripts/plot_score_celltype_html_v2.py \
  --h5ad data/processed/Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_merged_nontumor_v1_broad_ext_rabies_assignment.h5ad \
  --if_dir /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images/  \
  --if_image 846531_HD_BA_WM_C_DNB_SM_D90-7-D1-GFP_R1_IF_transform_reoriented.tif  \
  --gene Chrna4 \
  --celltype_col cellestimate \
  --celltype_col rabies_assignment \
  --highlight_celltype Glut \
  --highlight_celltype positive \
  --highlight_opacity 0.3 \
  --highlight_opacity 1.0 \
  --score_vmin 1.5  \
  --spot_size 3.0



