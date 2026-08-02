#conda activate scrna
library(SeuratDisk)
library(Seurat)
library(anndata)

ad <- read_h5ad("spatial_3356629_D49_D1_16um_test.h5ad")
obj <- CreateSeuratObject(
	     counts = t(ad$X), # AnnData is (cells x genes), Seurat needs (genes x cells)
	     meta.data = as.data.frame(ad$obs)
       )

# Spatial data in H5AD is usually stored in ad$obsm['spatial']
if ("spatial" %in% names(ad$obsm)) {
  message("Extracting spatial coordinates...")
  spatial_coords <- as.data.frame(ad$obsm[["spatial"]])
  colnames(spatial_coords) <- c("x", "y")
  rownames(spatial_coords) <- rownames(obj@meta.data)

  # Add coordinates to the Seurat object metadata
  obj$spatial_x <- spatial_coords$x
  obj$spatial_y <- spatial_coords$y

  # Optional: Create a formal Seurat spatial image/centroid object if needed
  # For Seurat v5, simple metadata coordinates often suffice for custom plotting,
  # but you can also build a Centroids object:
  # obj[["spatial"]] <- CreateCentroids(spatial_coords)
}
obj <- NormalizeData(obj, assay = "RNA")
obj <- FindVariableFeatures(obj, assay = "RNA")
obj <- ScaleData(obj, assay = "RNA")


saveRDS(obj, file = "spatial_3356629_D49_D1_16um_test.rds")
#seu <- LoadSeuratRds(rds_file) #LoadSeuratRds("GSM5900316_proj_raw.rds")
#counts  <- GetAssayData(seu, assay = "RNA", layer = "counts")

#genes <- rownames(counts)
#barcodes <- colnames(counts)
#meta     <- seu@meta.data



#rscript_path <- file.path(R.home("bin"), "Rscript")
#cell_type_script <- "/research/groups/mackgrp/home/demerson/CellTypeEstimate/cte4so.R"
#output <-  "spatial_3356629_D49_D1_16um_test.rds"
#print(paste(rscript_path, cell_type_script, ' --rds ',output,' --assay SCT --db Yao_et_al_mouse_format_addendum --cf Yao_et_al_mouse_format'))
#system(paste(rscript_path, cell_type_script, ' --rds ',output,' --assay SCT --db Yao_et_al_mouse_format_addendum --cf Yao_et_al_mouse_format'))
