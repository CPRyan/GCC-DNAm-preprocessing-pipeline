merge_chunk_csvs <- function(config, suffix, output_suffix = suffix) {
  files <- vapply(seq_len(as.integer(config$chunking$n_chunks %||% 1)), function(chunk) {
    study_file(config, suffix, chunk = chunk, ext = "csv")
  }, character(1))
  files <- files[file.exists(files)]
  if (!length(files)) {
    warning("No files found for suffix: ", suffix, call. = FALSE)
    return(invisible(NULL))
  }
  merged <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE, check.names = FALSE))
  out <- study_file(config, paste0(output_suffix, "_all_chunks"), ext = "csv")
  utils::write.csv(merged, out, row.names = FALSE)
  invisible(out)
}

merge_standard_outputs <- function(config) {
  list(
    qc_metrics = merge_chunk_csvs(config, "sample_qc_metrics"),
    qc_pass_fail = merge_chunk_csvs(config, "sample_qc_pass_fail"),
    qc_filtering_decisions = merge_chunk_csvs(config, "sample_qc_filtering_decisions")
  )
}
