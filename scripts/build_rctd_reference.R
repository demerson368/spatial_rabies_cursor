#!/usr/bin/env Rscript
# Build cached RCTD reference RDS from Python-exported Matrix Market bundle.

suppressPackageStartupMessages({
  library(Matrix)
  library(GetoptLong)
})

bundle_dir <- ""
out_rds <- ""

GetoptLong(
  "bundle_dir=s", "Directory with reference_*.mtx and sidecar files",
  "out_rds=s", "Output reference.rds path"
)

if (!nzchar(bundle_dir) || !nzchar(out_rds)) {
  stop("--bundle_dir and --out_rds are required", call. = FALSE)
}

counts <- Matrix::readMM(file.path(bundle_dir, "reference_counts.mtx"))
counts <- as(counts, "CsparseMatrix")
genes <- readLines(file.path(bundle_dir, "reference_genes.txt"))
barcodes <- readLines(file.path(bundle_dir, "reference_barcodes.txt"))
rownames(counts) <- genes
colnames(counts) <- barcodes

cell_type_df <- read.csv(file.path(bundle_dir, "reference_cell_types.csv"), stringsAsFactors = FALSE)
cell_types <- factor(cell_type_df$cell_type[match(barcodes, cell_type_df$barcode)])
names(cell_types) <- barcodes
nUMI <- Matrix::colSums(counts)
names(nUMI) <- barcodes

meta <- jsonlite::fromJSON(file.path(bundle_dir, "reference_meta.json"))

reference_inputs <- list(
  counts = counts,
  cell_types = cell_types,
  nUMI = nUMI,
  genes = genes,
  barcodes = barcodes,
  meta = meta
)

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(reference_inputs, out_rds)
message("Wrote reference RDS: ", out_rds)
