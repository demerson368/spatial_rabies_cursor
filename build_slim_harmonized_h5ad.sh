bsub -o build_slim_harmonized_h5ad.o -e build_slim_harmonized_h5ad.e -R "rusage[mem=150000]" python scripts/build_slim_harmonized_h5ad.py 
