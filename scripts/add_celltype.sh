bsub -o add_celltype.o -e add_celltype.e -R "rusage[mem=50000]" Rscript add_celltype.R
