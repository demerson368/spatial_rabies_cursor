#run_rctd.py includes run_rctd.R
#bsub -o run_rctd.o -e run_rctd.e -R "rusage[mem=100000]" python scripts/run_rctd.py \
#  --ref_out /research/groups/mackgrp/home/common/Analysis/seqData.Tuyu/RABV_scRNAseq_202507/Dan_Emerson_analysis/merged_dJ33_J51_D87_J58_J59_J61_D85_J52_broad/ \
#  --ref_label merged_dJ33_J51_D87_J58_J59_J61_D85_J52_broad \
#  --ref_file integrated_ref_rawcounts_with_labels.h5ad \
#  --ref_tag merged_nontumor_v1_broad \
#  --sr_out /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/Processed/846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1/outs \
#  --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um \
#  --cohort Tuyu_Apr2026 \
#  --bin_um 008um


#bsub -o run_rctd.o -e run_rctd.e -R "rusage[mem=100000]"  Rscript scripts/run_rctd.R --reference_rds results_rctd/ref_merged_dJ33_J51_D87_J58_J59_J61_D85_J52_broad_merged_nontumor_v1_broad/reference.rds --spatial_input_rds results_rctd/rctd_spatial_Tuyu_Apr2026_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad/spatial_input.rds --out_dir results_rctd/rctd_spatial_Tuyu_Apr2026_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad --doublet_mode doublet --max_cores 8 --umi_min 0 --umi_max 20000 --export_only TRUE



#bsub -o run_rctd.o -e run_rctd.e -R "rusage[mem=100000]" python scripts/run_rctd.py \
#  --ref_out /research/groups/mackgrp/home/common/Analysis/seqData.Tuyu/RABV_scRNAseq_202507/Dan_Emerson_analysis/merged_dJ33_J51_D87_J58_J59_J61_D85_J52_tuyusig_8_22/ \
#  --ref_label merged_dJ33_J51_D87_J58_J59_J61_D85_J52_tuyusig_8_22 \
#  --ref_file integrated_ref_rawcounts_with_labels.h5ad \
#  --ref_tag merged_dJ33_J51_D87_J58_J59_J61_D85_J52_tuyusig_8_22 \
#  --sr_out /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/Processed/846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1/outs \
#  --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um \
#  --cohort Tuyu_Apr2026 \
#  --bin_um 008um



#bsub -o run_rctd.o -e run_rctd.e -R "rusage[mem=100000]" python scripts/run_rctd.py \
#  --ref_out /research/groups/mackgrp/home/common/Analysis/seqData.Tuyu/RABV_scRNAseq_202507/Dan_Emerson_analysis/merged_tumor_v1_broad/ \
#  --ref_label merged_tumor_v1_broad \
#  --ref_file integrated_tumor_v1_broad_ref_rawcounts_with_labels.h5ad \
#  --ref_tag merged_tumor_v1_broad \
#  --sr_out /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/Processed/846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1/outs \
#  --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um \
#  --cohort Tuyu_Apr2026 \
#  --bin_um 008um



bsub -o run_rctd.o -e run_rctd.e -R "rusage[mem=100000]" python scripts/run_rctd.py \
  --ref_out /research/groups/mackgrp/home/common/Analysis/seqData.Tuyu/RABV_scRNAseq_202507/Dan_Emerson_analysis/merged_dJ33_J51_D87_J58_J59_J61_D85_J52_broad/ \
  --ref_label merged_dJ33_J51_D87_J58_J59_J61_D85_J52_broad \
  --ref_file integrated_ref_rawcounts_with_labels_neuron_ifn_drop.h5ad \
  --ref_tag merged_dJ33_J51_D87_J58_J59_J61_D85_J52_broad_neuron_ifn_drop \
  --sr_out /research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/Processed/846531_HD_BA_WM_C_DNB_SM_D90-5-A1-GFP_R1/outs \
  --label Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um \
  --cohort Tuyu_Apr2026 \
  --bin_um 008um

