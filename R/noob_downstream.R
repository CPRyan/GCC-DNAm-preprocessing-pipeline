estimate_cell_counts_noob <- function(config, chunk) {
  require_package("minfi")
  grset <- readRDS(noob_grset_file(config, chunk))
  counts <- minfi::estimateCellCounts2(grset)
  file <- study_file(config, "Noob_cell_counts", chunk = chunk, ext = "csv")
  utils::write.csv(as.data.frame(counts), file, row.names = TRUE)
  invisible(file)
}

run_pc_clocks_noob <- function(config, chunk) {
  stop("PC clock wrapper is not implemented yet. Configure external PC clock functions before enabling this step.", call. = FALSE)
}

run_dunedinpace_noob <- function(config, chunk) {
  stop("DunedinPACE wrapper is not implemented yet. Configure DunedinPACE before enabling this step.", call. = FALSE)
}

run_noob_downstream <- function(config, chunk) {
  opts <- config$noob_downstream %||% list()
  outputs <- list()
  if (isTRUE(opts$cell_counts)) {
    outputs$cell_counts <- estimate_cell_counts_noob(config, chunk)
  }
  if (isTRUE(opts$pc_clocks)) {
    outputs$pc_clocks <- run_pc_clocks_noob(config, chunk)
  }
  if (isTRUE(opts$dunedinpace)) {
    outputs$dunedinpace <- run_dunedinpace_noob(config, chunk)
  }
  invisible(outputs)
}
