bsub -o filtering_postprocess.o -e filtering_postprocess.e -R "rusage[mem=50000]" Rscript filtering_postprocess.R --outs_path GSE196771_RAW/GSM5900316/ --filtered_file filtered_feature_bc_matrix --raw_file raw_feature_bc_matrix --proj_name GSM5900316_proj

