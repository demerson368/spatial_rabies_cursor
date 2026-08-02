# spatial_rabies

A spatial transcriptomics workflow for analyzing rabies virus infection patterns in mouse brain.

## Overview

This project provides:

- Integration of space ranger transcription and spatial output into h5 object
- Basic QC and normalization of integrated h5 object
- Spot annotation against class marker genes
   
### Set up environment
conda env create -f spatial_environment.yml  
conda activate spatial


### Step1: process and normalize spatial counts

```bash
python scripts/spatial_pipeline.py {spaceranger outs pathway}  {name_label}  
```

output:  
spatial_*.h5ad - contains log normalized low counts filtered spots, spatial_autocorrelation, and spatial data  


### Step2: cell2location

Use raw scRNA raw counts and raw spot counts (step 1 not necessary)

```bash

python run_cell2location_v2.py --ref_out {ref_dir} --ref_label {ref_label}  --ref_file {ref h5ad} --sr_out {spatial_outs} --label  {name_label}
sh scripts/run_cell2location_v2.sh
```


### Step3: celltype label transfer method

Note:  Not working well so far  
convert spatial .h5ad output in step 1 to seurat .rds and perform label transfer from ref scRNA to spatial.    
Using the reference scRNA merged scRNA nontumor:  
integrated_dJ33_J51_D87_J58_J59_J61_D85_J52_J53_ambient_corrected_v2_mito_doublet_removal_SCT_update.rds  


```bash
Rscript create_seurat_from_h5ad.R --h5ad {h5ad}
sh scripts/create_seurat_from_h5ad.sh

Rscript label_transfer.R --vis_file {converted_rds}
sh scripts/label_transfer.sh

```

create transfer label output as "*predicted.id_transfer.txt" located at data/processed/Tuyu_Apr2026/.   
Label spots from transfer prediction from scRNA 


### Step4: cellestimate method on spatial

use .rds object (step 3) conversion from .h5ad (step 1) with Spatial assay  


```bash
/research/groups/mackgrp/home/demerson/CellTypeEstimate/cte4so.R --rds {../data/processed/*/spatial_*.rds} --db {Yao_et_al_mouse_format_addendum} --assay Spatial --cf {Yao_et_al_mouse_format}

sh scripts/run_celltype.sh

```

cellestimate output:  
data/processed/*/spatial_*Yao_et_al_mouse_format_cellType.txt  


### Step5: complete celltype methods and integrate into h5ad and for_loupe

create or import celltype on spots either direct calling or transfer/modeling:

*tlscore* - score spot (current, step 5)  
*tlscore_norm* - score on spot. Normalized (current, step 5)  
*NNLS* - deconvolution on spot.  set to max fraction per spot (current, step 5)  
*cellestimate* - import cellestimate calls per spot  (step 4)  
*predicted.id* - import celltype transfer from scRNA to spot (step 3)  
*cell2location* - import cell2location calls (step 2)  


```bash
python assemble_loupe_file.py --h5ad {h5ad} --label {label}
sh scripts/assemble_loupe_file.sh
```

add features to .h5ad and export loupe ../for_loupe/*.csv   


## Snakemake workflow (orchestration)

Install Snakemake into the **`spatial`** env (once):

```bash
conda activate spatial
pip install -r config/snakemake_requirements.txt
```

The repo includes a Snakemake scaffold that wraps the steps above. Run from the **project root** on the cluster:

```bash
conda activate spatial
cd /research/groups/mackgrp/home/demerson/spatial_rabies_cursor

# Dry run (login node — use 1 core)
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --cores 1

# Submit to LSF (profile sets cores: 10, jobs: 50)
snakemake -s workflow/Snakefile --configfile config/config.yaml --profile profiles/lsf

# One sample end-to-end
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 1 \
  for_loupe/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um_for_loupe.csv
```

| Path | Purpose |
|------|---------|
| `config/samples.tsv` | Sample × bin matrix (edit to add samples) |
| `config/config.yaml` | Reference paths, conda env (`spatial`) |
| `config/snakemake_requirements.txt` | Snakemake + LSF plugin + compatible PuLP |
| `workflow/Snakefile` | Main DAG |
| `workflow/rules/*.smk` | One rule file per pipeline stage |
| `profiles/lsf/config.yaml` | LSF memory/GPU defaults |

**Troubleshooting:**
- `pulp` has no attribute `list_solvers` — pin PuLP: `pip install 'pulp<2.8'` (see `config/snakemake_requirements.txt` for Python 3.10).
- `you need to specify the maximum number of CPU cores` — add `--cores 1` (dry run) or use `--profile profiles/lsf` (sets `cores: 10`).

**Notes:**
- `assemble_loupe_file.py` currently requires cell2location output; samples with `run_cell2location=false` are excluded from `rule all`.
- `run_cell2location_v2.py` hardcodes `008um` bins — only enable cell2location for 008µm rows in `samples.tsv`.
