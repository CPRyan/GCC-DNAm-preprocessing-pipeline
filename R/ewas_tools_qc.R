run_ewas_tools_qc <- function(config) {
  if (!isTRUE(config$ewas_tools_qc$run)) {
    stop("EWAS tools QC is disabled in config, but this pipeline requires it.", call. = FALSE)
  }

  manifest <- if (file.exists(chunks_file(config))) {
    read_chunk_assignments(config)
  } else {
    load_and_validate_manifest(config)
  }

  qc <- control_metric_qc_from_manifest(manifest)
  saveRDS(qc, ewas_qc_file(config))
  qc
}

control_metric_qc_from_manifest <- function(manifest) {
  require_package("ewastools")

  if (!"Basename" %in% colnames(manifest)) {
    stop("Manifest must contain Basename for EWAS tools QC.", call. = FALSE)
  }
  if (!"sample_id" %in% colnames(manifest)) {
    stop("Manifest must contain sample_id for EWAS tools QC.", call. = FALSE)
  }

  basenames <- manifest$Basename
  sample_ids <- manifest$sample_id

  data_for_control_metrics <- ewastools::read_idats(basenames, quiet = FALSE)
  control_test_values <- ewastools::control_metrics(data_for_control_metrics)

  passed_qc <- data.frame(matrix(ncol = length(control_test_values), nrow = length(sample_ids)))
  qc_result_values <- data.frame(matrix(ncol = length(control_test_values), nrow = length(sample_ids)))
  rownames(passed_qc) <- sample_ids
  rownames(qc_result_values) <- sample_ids
  colnames(passed_qc) <- names(control_test_values)
  colnames(qc_result_values) <- names(control_test_values)

  for (metric_name in names(control_test_values)) {
    metric_value <- as.numeric(control_test_values[[metric_name]])
    threshold <- attr(control_test_values[[metric_name]], "threshold")
    if (length(metric_value) != length(sample_ids)) {
      stop(
        "EWAS tools metric length did not match manifest sample count for metric: ",
        metric_name,
        call. = FALSE
      )
    }
    passed_qc[[metric_name]] <- metric_value >= threshold
    qc_result_values[[metric_name]] <- metric_value
    attr(qc_result_values[[metric_name]], "threshold") <- threshold
  }

  list(
    passed_QC = passed_qc,
    QC_result_values = qc_result_values
  )
}

load_or_run_ewas_tools_qc <- function(config) {
  file <- ewas_qc_file(config)
  if (file.exists(file)) {
    return(readRDS(file))
  }
  run_ewas_tools_qc(config)
}

as_ewas_table <- function(x, sample_id_col = "sample_id") {
  if (is.null(x)) {
    return(NULL)
  }
  out <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  out[[sample_id_col]] <- rownames(out)
  rownames(out) <- NULL
  out
}
