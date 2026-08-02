#!/usr/bin/env Rscript
# Seurat label transfer: scRNA reference → Visium spatial query.
# Paths are passed on the CLI (no hardcoded setwd).

library(Seurat)
library(GetoptLong)

# Defaults (overridden by CLI / Snakemake)
project_dir <- "/research/groups/mackgrp/home/demerson/spatial_rabies_cursor"
vis_file <- ""
ref_rds <- ""
out_file <- ""

GetoptLong(
  "project_dir=s", "project root (spatial_rabies_cursor)",
  "vis_file=s", "absolute or project-relative path to spatial .rds",
  "ref_rds=s", "absolute path to reference Seurat .rds",
  "out_file=s", "absolute path for predicted.id transfer table"
)

if (!nzchar(vis_file)) {
  stop("--vis_file is required")
}
if (!nzchar(ref_rds)) {
  stop("--ref_rds is required")
}

# Resolve vis_file relative to project if needed
if (!startsWith(vis_file, "/")) {
  vis_file <- file.path(project_dir, vis_file)
}
if (!nzchar(out_file)) {
  out_file <- sub("\\.rds$", "_predicted.id_transfer.txt", vis_file)
}

message("vis_file: ", vis_file)
message("ref_rds:  ", ref_rds)
message("out_file: ", out_file)

if (!file.exists(vis_file)) {
  stop("Spatial RDS not found: ", vis_file)
}
if (!file.exists(ref_rds)) {
  stop("Reference RDS not found: ", ref_rds)
}

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

ref <- LoadSeuratRds(ref_rds)
vis <- LoadSeuratRds(vis_file)

DefaultAssay(vis) <- "Spatial"
vis <- NormalizeData(vis, verbose = FALSE)
vis <- FindVariableFeatures(vis, verbose = FALSE)
vis <- ScaleData(vis, verbose = FALSE)
vis <- RunPCA(vis, verbose = FALSE)
# Skip UMAP/clustering — not required for TransferData and slows HD jobs

DefaultAssay(ref) <- "RNA"
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, verbose = FALSE)
ref <- ScaleData(ref, verbose = FALSE)
ref <- RunPCA(ref, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = ref,
  query = vis,
  normalization.method = "LogNormalize",
  dims = 1:30
)
predictions <- TransferData(
  anchorset = anchors,
  refdata = ref$cell_type,
  dims = 1:30
)

vis <- AddMetaData(vis, metadata = predictions)

# Non-interactive: skip SpatialDimPlot (can fail / hang without a display)
# SpatialDimPlot(vis, group.by = "predicted.id")

write.table(
  vis@meta.data["predicted.id"],
  file = out_file,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
message("Wrote: ", out_file)
