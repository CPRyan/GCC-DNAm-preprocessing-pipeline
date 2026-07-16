source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

manifest <- load_and_validate_manifest(config)
if (any(!manifest$idat_pair_exists)) {
  stop("Cannot create targets until all IDAT pairs exist. Run scripts/00_validate_inputs.R for details.", call. = FALSE)
}

n_chunks <- if (isTRUE(config$chunking$enabled)) as.integer(config$chunking$n_chunks %||% 1) else 1L
manifest <- create_chunk_assignments(
  manifest,
  n_chunks = n_chunks,
  grouping_cols = config$chunking$grouping_cols,
  seed = as.integer(config$chunking$seed %||% 123)
)

utils::write.csv(manifest, targets_file(config), row.names = FALSE)
utils::write.csv(manifest, chunks_file(config), row.names = FALSE)

message("Created ", n_chunks, " chunk(s).")
message("Wrote: ", chunks_file(config))
