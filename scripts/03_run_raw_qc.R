source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

message("Loading or running mandatory EWAS tools QC.")
ewas_qc <- load_or_run_ewas_tools_qc(config)
ctrl_dir <- NULL

for (chunk in get_chunks_to_run(config)) {
  message("Running raw QC for chunk ", chunk)
  rgset <- readRDS(raw_rgset_file(config, chunk))
  manifest <- manifest_for_chunk(config, chunk)
  raw_qc <- run_raw_qc_for_rgset(config, rgset)
  build_qc_tables(config, chunk, manifest, raw_qc, ewas_qc)

  # save this chunk's control-probe betas; PCA runs on all chunks after the loop
  beta_file <- study_file(config, "control_probe_beta", chunk = chunk, ext = "rds")
  saveRDS(raw_qc$control_beta, beta_file)
  ctrl_dir <- dirname(beta_file)

  try(plot_qc_upset(config, chunk), silent = TRUE)
  message("Wrote QC tables for chunk ", chunk)
  rm(rgset, raw_qc)
  gc()
}

# ---- control-probe PCA across all chunks ----
if (is.null(ctrl_dir)) {
  message("No chunks run this session; skipping merged control-probe PCA.")
} else {
  beta_files <- list.files(ctrl_dir, pattern = "control_probe_beta.*\\.rds$", full.names = TRUE)
  message("Merging control-probe betas from ", length(beta_files), " chunk file(s).")
  betas  <- lapply(beta_files, readRDS)
  common <- Reduce(intersect, lapply(betas, rownames))
  control_beta_all <- do.call(cbind, lapply(betas, function(b) b[common, , drop = FALSE]))
  if (anyDuplicated(colnames(control_beta_all))) stop("Duplicate sample_ids across chunks", call. = FALSE)
  
  control_pca <- calculate_control_probe_pca_from_beta(
    control_beta_all,
    config$qc_thresholds$control_pca_variance_threshold %||% 0.90
  )
  message("  ", ncol(control_beta_all), " samples, ", control_pca$n_control_probes_used,
          " control probes, ", control_pca$n_pcs_to_threshold, " PCs to reach threshold")
  
  utils::write.csv(control_pca$scores,   file.path(ctrl_dir, "control_probe_pca_scores_all_chunks.csv"),   row.names = FALSE)
  utils::write.csv(control_pca$variance, file.path(ctrl_dir, "control_probe_pca_variance_all_chunks.csv"), row.names = FALSE)
  saveRDS(control_pca$pca, file.path(ctrl_dir, "control_probe_pca_all_chunks.rds"))
}