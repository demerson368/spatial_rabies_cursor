# Analyze GFP+ neighborhoods (default annotation: cellestimate) and optionally plot.
#python scripts/gfp_neighbor_cell2location.py \
#  --input_dir for_loupe/ \
#  --output_dir results/Tuyu_Apr2026/ \
#  --annotation_col cellestimate \
#  --write_per_seed \
#  --plot

# Explicit cell2location annotation:
# python scripts/gfp_neighbor_cell2location.py \
#   --annotation_col cell2location \
#   --output_dir results/Tuyu_Apr2026/

# Plot only (from existing *_gfp_neighbor_{slug}_all_radii.csv files):
python scripts/gfp_neighbor_cell2location.py \
   --plot_only \
   --annotation_col cellestimate \
   --samples Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um
   --top_n 10 \
   --output_dir results/Tuyu_Apr2026/ \
   --plot_metric fraction


#python scripts/gfp_neighbor_cell2location.py \
#   --plot_only \
#   --annotation_col cell2location \
#   --output_dir results/Tuyu_Apr2026/ \
#   --top_n 10 \
#   --plot_metric fraction



# Permutation test for neighborhood enrichment (default: cellestimate):
# python scripts/gfp_neighbor_cell2location_permutation.py \
#   --input_dir for_loupe/ \
#   --output_dir results/Tuyu_Apr2026/ \
#   --annotation_col cellestimate \
#   --n_perm 999 \
#   --seed 42
