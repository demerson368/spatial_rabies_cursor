# Analyze GFP+ neighborhoods and optionally plot top cell2location bar charts.
#python scripts/gfp_neighbor_cell2location.py \
#  --input_dir for_loupe/ \
#  --output_dir results/Tuyu_Apr2026/ \
#  --write_per_seed \
#  --plot

# Plot only (from existing *_gfp_neighbor_c2l_all_radii.csv files):
python scripts/gfp_neighbor_cell2location.py \
   --plot_only \
   --output_dir results/Tuyu_Apr2026/ \
   --top_n 10 \
   --plot_metric fraction
