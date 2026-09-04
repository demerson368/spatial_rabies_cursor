#conda activate scrna45
library(Seurat)
library(SeuratObject)
library(Matrix)
library(reticulate)
library(GetoptLong)

reticulate::use_condaenv("scrna45", required = TRUE)
ad <- import("anndata")
sc <- import("scipy.sparse")

#default
h5ad <- "Tuyu_Apr2026/spatial_Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um.h5ad"
out <- ""

GetoptLong(
    "h5ad=s",       "input spatial h5ad for conversion to seurat",
    "out=s",        "output Seurat rds path (default: h5ad with .rds extension)"
)

output <- if (nzchar(out)) out else sub("\\.h5ad$", ".rds", h5ad)
adata <- ad$read_h5ad(h5ad)
obs_names <- unlist(adata$obs_names$to_list(), use.names = FALSE)
var_names <- unlist(adata$var_names$to_list(), use.names = FALSE)

# Pull X into R once (may already be dgCMatrix/dgRMatrix/matrix)
X <- py_to_r(adata$X)
cat("X class:", paste(class(X), collapse = ", "), " dims:", paste(dim(X), collapse = " x "), "\n")

print("counts creation")

# AnnData X = cells x genes; Seurat = genes x cells
Xt <- Matrix::t(X)
print("dgCMatrix")
counts <- as(as(Xt, "TsparseMatrix"), "CsparseMatrix")  # -> dgCMatrix

colnames(counts) <- obs_names
rownames(counts) <- var_names
stopifnot(inherits(counts, "dgCMatrix"))
# Optional: drop empty bins early
counts <- counts[, Matrix::colSums(counts) > 0, drop = FALSE]

print("create seu")
seu <- CreateSeuratObject(
  counts = counts,
  assay = "Spatial",
  min.cells = 0,
  min.features = 0
)


print("add meta")
# Add obs metadata safely (row order matches cells)
meta <- py_to_r(adata$obs)
rownames(meta) <- obs_names
seu <- AddMetaData(seu, metadata = meta[colnames(seu), , drop = FALSE])

print("match coords")
# Spatial coordinates — rownames MUST match cell barcodes
coords <- py_to_r(adata$obsm["spatial"])
coords <- as.data.frame(coords)
colnames(coords) <- c("x", "y")
rownames(coords) <- obs_names

coords <- coords[colnames(seu), , drop = FALSE]  # same cells, same order
stopifnot(all(rownames(coords) == colnames(seu)))

print("fov creation")
fov <- CreateFOV(
  coords = coords,
  type = "centroids",
  assay = "Spatial"
)

print("assign fov")
# Use a named FOV slot (not "spatial")
seu[["slice1"]] <- fov

#table(seu$nCount_Spatial > 0)
#TRUE 
#78284

threshold <- 0
subset(seu, nCount_Spatial >= threshold)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.8)


SaveSeuratRds(seu, file = output)



