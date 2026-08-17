#!/usr/bin/env Rscript
# One-time install for R packages used by CellTypeEstimate and Snakemake R rules.
# Run after creating/recreating the spatial conda env:
#   /home/demerson/.conda/envs/spatial/bin/Rscript scripts/install_r_pipeline_packages.R

cran <- "https://cloud.r-project.org"

need_cran <- c("HGNChelper", "GetoptLong", "readxl", "openxlsx")
missing <- need_cran[!vapply(need_cran, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  message("Installing from CRAN: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = cran)
} else {
  message("All CRAN packages already installed.")
}

for (pkg in need_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package still missing after install: ", pkg, call. = FALSE)
  }
  message("OK: ", pkg)
}
