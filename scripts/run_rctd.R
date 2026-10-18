#!/usr/bin/env Rscript
# Run RCTD (spacexr) on prepared reference and spatial RDS inputs.

suppressPackageStartupMessages({
  library(Matrix)
  library(GetoptLong)
})

rctd_weights_to_df <- function(w, barcodes, cell_type_names) {
  if (is.data.frame(w)) {
    dense <- as.matrix(w)
  } else if (inherits(w, "sparseMatrix") || is.matrix(w)) {
    dense <- as.matrix(w)
  } else {
    dense <- as.matrix(w)
  }
  n_spots <- length(barcodes)
  if (nrow(dense) == n_spots) {
    weights <- as.data.frame(dense, stringsAsFactors = FALSE)
  } else if (ncol(dense) == n_spots) {
    weights <- as.data.frame(t(dense), stringsAsFactors = FALSE)
  } else {
    stop(
      "RCTD weights dimensions (", paste(dim(dense), collapse = " x "), ") ",
      "do not match n_spots (", n_spots, ")", call. = FALSE
    )
  }
  if (is.null(rownames(weights))) {
    rownames(weights) <- barcodes
  }
  if (ncol(weights) != length(cell_type_names)) {
    stop(
      "RCTD weights has ", ncol(weights), " columns but expected ",
      length(cell_type_names), " cell types", call. = FALSE
    )
  }
  colnames(weights) <- cell_type_names
  weights
}

write_rctd_outputs <- function(
  my_rctd,
  spat_inputs,
  ref_inputs,
  out_dir,
  doublet_mode,
  max_cores,
  umi_min,
  umi_max
) {
  results <- my_rctd@results
  cell_type_names <- my_rctd@cell_type_info$info[[2]]
  # Spot x cell-type matrix (dgCMatrix in recent spacexr). weights_doublet is n_spots x 2.
  weights_src <- results$weights
  if (is.null(weights_src)) {
    stop("RCTD results missing weights matrix", call. = FALSE)
  }
  weights <- rctd_weights_to_df(weights_src, spat_inputs$barcodes, cell_type_names)
  write.csv(weights, file.path(out_dir, "rctd_weights.csv"), row.names = TRUE)

  results_df <- results$results_df
  if (!is.null(results_df)) {
    write.csv(results_df, file.path(out_dir, "rctd_results_df.csv"), row.names = TRUE)
  }

  run_meta <- list(
    doublet_mode = doublet_mode,
    max_cores = max_cores,
    umi_min = umi_min,
    umi_max = umi_max,
    n_pixels = length(spat_inputs$barcodes),
    n_genes = nrow(spat_inputs$counts),
    cell_type_names = cell_type_names,
    spatial_meta = spat_inputs$meta,
    reference_meta = ref_inputs$meta
  )
  jsonlite::write_json(
    run_meta,
    file.path(out_dir, "rctd_run_meta.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )
}

reference_rds <- ""
spatial_input_rds <- ""
out_dir <- ""
doublet_mode <- "doublet"
max_cores <- 8L
umi_min <- 0L
umi_max <- 20000L
export_only <- FALSE

GetoptLong(
  "reference_rds=s", "Cached reference RDS from build_rctd_reference.R",
  "spatial_input_rds=s", "Cached spatial input RDS",
  "out_dir=s", "Output directory for RCTD results",
  "doublet_mode=s", "RCTD mode: doublet or full (default: doublet)",
  "max_cores=i", "Number of cores for RCTD (default: 8)",
  "umi_min=i", "Minimum UMI per pixel (default: 0)",
  "umi_max=i", "Maximum UMI per pixel (default: 20000)",
  "export_only!", "Only export CSV/JSON from existing rctd_object.rds (skip run.RCTD)"
)

if (any(grepl("^--export_only($|=)", commandArgs(trailingOnly = TRUE)))) {
  export_only <- TRUE
}

if (!nzchar(reference_rds) || !nzchar(spatial_input_rds) || !nzchar(out_dir)) {
  stop("--reference_rds, --spatial_input_rds, and --out_dir are required", call. = FALSE)
}

doublet_mode <- tolower(doublet_mode)
if (!doublet_mode %in% c("doublet", "full")) {
  stop("doublet_mode must be 'doublet' or 'full'", call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ref_inputs <- readRDS(reference_rds)
spat_inputs <- readRDS(spatial_input_rds)

if (export_only) {
  rctd_path <- file.path(out_dir, "rctd_object.rds")
  if (!file.exists(rctd_path)) {
    stop("export_only requires existing file: ", rctd_path, call. = FALSE)
  }
  message("Export only: loading ", rctd_path)
  my_rctd <- readRDS(rctd_path)
  write_rctd_outputs(
    my_rctd,
    spat_inputs,
    ref_inputs,
    out_dir,
    doublet_mode,
    max_cores,
    umi_min,
    umi_max
  )
  message("RCTD export complete (export_only). Outputs in: ", out_dir)
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages(library(spacexr))

reference <- Reference(ref_inputs$counts, ref_inputs$cell_types, ref_inputs$nUMI)
puck <- SpatialRNA(spat_inputs$coords, spat_inputs$counts, spat_inputs$nUMI)

message(
  "Running RCTD: pixels=", ncol(spat_inputs$counts),
  ", genes=", nrow(spat_inputs$counts),
  ", cell_types=", length(levels(ref_inputs$cell_types)),
  ", mode=", doublet_mode
)

my_rctd <- create.RCTD(
  puck,
  reference,
  max_cores = as.integer(max_cores),
  UMI_min = as.integer(umi_min),
  UMI_max = as.integer(umi_max)
)
my_rctd <- run.RCTD(my_rctd, doublet_mode = doublet_mode)

saveRDS(my_rctd, file.path(out_dir, "rctd_object.rds"))

write_rctd_outputs(
  my_rctd,
  spat_inputs,
  ref_inputs,
  out_dir,
  doublet_mode,
  max_cores,
  umi_min,
  umi_max
)

message("RCTD complete. Outputs in: ", out_dir)
