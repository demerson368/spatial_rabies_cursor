# Permutation test for neighborhood cell2location enrichment:
bsub -o gfp_neighbor_cell2location_permutation.o -e gfp_neighbor_cell2location_permutation.e -R "rusage[mem=100000]" python scripts/gfp_neighbor_cell2location_permutation.py --input_dir for_loupe/ --output_dir results/Tuyu_Apr2026/ --n_perm 999 --seed 42
