#!/usr/bin/env Rscript
# Build cached RCTD spatial input RDS from Python-exported Matrix Market bundle.

suppressPackageStartupMessages({
  library(Matrix)
  library(GetoptLong)
})

bundle_dir <- ""
out_rds <- ""

GetoptLong(
  "bundle_dir=s", "Directory with spatial_*.mtx and sidecar files",
  "out_rds=s", "Output spatial_input.rds path"
)

if (!nzchar(bundle_dir) || !nzchar(out_rds)) {
  stop("--bundle_dir and --out_rds are required", call. = FALSE)
}

counts <- Matrix::readMM(file.path(bundle_dir, "spatial_counts.mtx"))
counts <- as(counts, "CsparseMatrix")
genes <- readLines(file.path(bundle_dir, "spatial_genes.txt"))
barcodes <- readLines(file.path(bundle_dir, "spatial_barcodes.txt"))
rownames(counts) <- genes
colnames(counts) <- barcodes

coords <- read.csv(file.path(bundle_dir, "spatial_coords.csv"), row.names = 1, check.names = FALSE)
coords <- coords[barcodes, , drop = FALSE]
colnames(coords) <- c("x", "y")

nUMI <- Matrix::colSums(counts)
names(nUMI) <- barcodes

meta <- jsonlite::fromJSON(file.path(bundle_dir, "spatial_meta.json"))

spatial_inputs <- list(
  counts = counts,
  coords = coords,
  nUMI = nUMI,
  barcodes = barcodes,
  meta = meta
)

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(spatial_inputs, out_rds)
message("Wrote spatial input RDS: ", out_rds)
