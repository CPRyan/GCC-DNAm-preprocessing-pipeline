safe_study_name <- function(config) {
  study_name <- config$study_name %||% "DNAmStudy"
  gsub("[^A-Za-z0-9_]+", "_", study_name)
}

output_root <- function(config) {
  config_path(config, config$paths$output_dir %||% "Output")
}

output_data_dir <- function(config) {
  file.path(output_root(config), "data")
}

output_figure_dir <- function(config) {
  file.path(output_root(config), "figures")
}

output_log_dir <- function(config) {
  file.path(output_root(config), "logs")
}

ensure_output_dirs <- function(config) {
  dirs <- c(output_root(config), output_data_dir(config), output_figure_dir(config), output_log_dir(config))
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  invisible(dirs)
}

study_file <- function(config, suffix, chunk = NULL, ext = NULL, dir = output_data_dir(config)) {
  stem <- safe_study_name(config)
  if (!is.null(chunk)) {
    stem <- paste0(stem, "_", suffix, "_chunk", chunk)
  } else {
    stem <- paste0(stem, "_", suffix)
  }
  if (!is.null(ext)) {
    stem <- paste0(stem, ".", ext)
  }
  file.path(dir, stem)
}

targets_file <- function(config) study_file(config, "targets", ext = "csv")
chunks_file <- function(config) study_file(config, "targets_chunks", ext = "csv")
validated_manifest_file <- function(config) study_file(config, "sample_manifest_validated", ext = "csv")
ewas_qc_file <- function(config) study_file(config, "ewas_tools_qc", ext = "rds")
raw_rgset_file <- function(config, chunk) study_file(config, "RGset_raw", chunk = chunk, ext = "rds")
good_rgset_file <- function(config, chunk) study_file(config, "RGset_good_samples", chunk = chunk, ext = "rds")
noob_grset_file <- function(config, chunk) study_file(config, "Noob_GRset", chunk = chunk, ext = "rds")
qc_metrics_file <- function(config, chunk) study_file(config, "sample_qc_metrics", chunk = chunk, ext = "csv")
qc_pass_fail_file <- function(config, chunk) study_file(config, "sample_qc_pass_fail", chunk = chunk, ext = "csv")
qc_filtering_decisions_file <- function(config, chunk) study_file(config, "sample_qc_filtering_decisions", chunk = chunk, ext = "csv")
probe_filtered_grset_file <- function(config, chunk) study_file(config, "probe_filtered_GRset", chunk = chunk, ext = "rds")
bmiq_betas_file <- function(config, chunk) study_file(config, "BMIQ_betas", chunk = chunk, ext = "rds")
