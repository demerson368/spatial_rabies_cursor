#snakemake -s workflow/Snakefile --configfile config/config.yaml \
#  --profile profiles/lsf \
#  results_cell2location/c2l_spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_v2/sp.h5ad \
#  --forcerun run_cell2location


#snakemake -s workflow/Snakefile --configfile config/config.yaml --profile profiles/lsf results_cell2location/c2l_spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_v2/sp.h5ad --forcerun run_cell2location --cluster "$(pwd)/profiles/lsf/cluster-submit.sh {dependencies}" --cluster-status "$(pwd)/profiles/lsf/status_lsf.sh {jobid}" --immediate-submit --notemp --cores 10 --jobs 50


#snakemake -s workflow/Snakefile --configfile config/config.yaml --profile profiles/lsf results_cell2location/c2l_spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_v2/sp.h5ad --forcerun run_cell2location --cluster "$(pwd)/profiles/lsf/cluster-submit.sh {dependencies}" --cluster-status "$(pwd)/profiles/lsf/status_lsf.sh {jobid}" --immediate-submit --notemp --cores 10 --jobs 50

snakemake -s workflow/Snakefile --configfile config/config.yaml --profile profiles/lsf results_cell2location/c2l_spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_v2/sp.h5ad --forcerun run_cell2location --cluster "$(pwd)/profiles/lsf/cluster-submit.sh {dependencies}" --cluster-status "$(pwd)/profiles/lsf/status_lsf.sh {jobid}" --immediate-submit --notemp --cores 10 --jobs 50
