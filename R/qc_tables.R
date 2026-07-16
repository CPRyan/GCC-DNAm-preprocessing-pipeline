merge_by_sample <- function(...) {
  tables <- list(...)
  tables <- Filter(function(x) !is.null(x) && nrow(x) > 0, tables)
  if (!length(tables)) {
    return(data.frame())
  }
  Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE, sort = FALSE), tables)
}

build_qc_tables <- function(config, chunk, manifest, raw_qc, ewas_qc) {
  ewas_pass <- as_ewas_table(ewas_qc$passed_QC)
  ewas_values <- as_ewas_table(ewas_qc$QC_result_values)

  manifest_small <- manifest[, intersect(c("sample_id", "Basename", "chunk"), colnames(manifest)), drop = FALSE]
  metrics <- merge_by_sample(manifest_small, raw_qc$signal, raw_qc$beads, raw_qc$detection, ewas_values)

  standard_pass <- merge_by_sample(
    raw_qc$signal[, c("sample_id", "passed_signal_quality"), drop = FALSE],
    raw_qc$beads[, c("sample_id", "passed_beadcount"), drop = FALSE],
    raw_qc$detection[, c("sample_id", "passed_detection_cutoff"), drop = FALSE]
  )
  pass_fail <- merge_by_sample(manifest_small, standard_pass, ewas_pass)

  utils::write.csv(metrics, qc_metrics_file(config, chunk), row.names = FALSE)
  utils::write.csv(pass_fail, qc_pass_fail_file(config, chunk), row.names = FALSE)

  list(metrics = metrics, pass_fail = pass_fail)
}
