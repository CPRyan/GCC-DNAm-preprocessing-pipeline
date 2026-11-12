qc_value_failed <- function(x) {
  if (is.logical(x)) {
    return(is.na(x) | !x)
  }
  lx <- tolower(as.character(x))
  is.na(x) | lx %in% c("fail", "failed", "false", "f", "0", "no", "n")
}

selected_filter_columns <- function(config, pass_fail) {
  standard_map <- c(
    signal_quality = "passed_signal_quality",
    beadcount = "passed_beadcount",
    bisulfite_conversion = "passed_bisulf_conv",
    detection_p = "passed_detection_cutoff",
    snp_outlier = "passed_snp_outlier_bl4",
    sex_mismatch = "passed_sex_mismatch"
  )

  selected <- character()
  standard <- config$sample_filtering$standard_metrics %||% list()
  for (metric in names(standard)) {
    col <- standard_map[[metric]]
    if (isTRUE(standard[[metric]]) && !is.null(col) && col %in% colnames(pass_fail)) {
      selected <- c(selected, col)
    }
  }

  ewas <- config$sample_filtering$ewas_tools_metrics %||% list()
  for (metric in names(ewas)) {
    if (isTRUE(ewas[[metric]]) && metric %in% colnames(pass_fail)) {
      selected <- c(selected, metric)
    }
  }

  unique(selected)
}

# One list entry per row: the names of the columns that failed (character(0) if none).
# Tests each column in its own type, and always returns one entry per sample.
failed_names <- function(df) {
  failed <- lapply(df, qc_value_failed)
  lapply(seq_len(nrow(df)), function(i) {
    colnames(df)[vapply(failed, function(f) f[[i]], logical(1))]
  })
}

apply_sample_filtering <- function(config, chunk) {
  pass_fail <- utils::read.csv(qc_pass_fail_file(config, chunk), stringsAsFactors = FALSE, check.names = FALSE)
  selected_cols <- selected_filter_columns(config, pass_fail)
  if (!"sample_id" %in% colnames(pass_fail)) {
    stop("Pass/fail table is missing sample_id.", call. = FALSE)
  }
  
  detected_cols <- setdiff(colnames(pass_fail), c("sample_id", "Basename", "chunk"))
  failed_detected <- failed_names(pass_fail[, detected_cols, drop = FALSE])
  failed_used     <- failed_names(pass_fail[, selected_cols, drop = FALSE])
  
  decisions <- data.frame(
    sample_id = pass_fail$sample_id,
    failed_metrics_detected = vapply(failed_detected, paste, character(1), collapse = ";"),
    failed_metrics_used_for_filtering = vapply(failed_used, paste, character(1), collapse = ";"),
    removed_from_pipeline = lengths(failed_used) > 0,
    stringsAsFactors = FALSE
  )
  utils::write.csv(decisions, qc_filtering_decisions_file(config, chunk), row.names = FALSE)
  decisions
}

save_good_sample_rgset <- function(config, chunk) {
  require_package("minfi")
  decisions <- apply_sample_filtering(config, chunk)
  rgset <- readRDS(raw_rgset_file(config, chunk))
  keep <- decisions$sample_id[!decisions$removed_from_pipeline]
  saveRDS(rgset[, colnames(rgset) %in% keep], good_rgset_file(config, chunk))
  invisible(good_rgset_file(config, chunk))
}
