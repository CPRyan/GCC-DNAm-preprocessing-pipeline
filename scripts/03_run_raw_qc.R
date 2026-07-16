source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

message("Loading or running mandatory EWAS tools QC.")
ewas_qc <- load_or_run_ewas_tools_qc(config)

for (chunk in get_chunks_to_run(config)) {
  message("Running raw QC for chunk ", chunk)
  rgset <- readRDS(raw_rgset_file(config, chunk))
  manifest <- manifest_for_chunk(config, chunk)
  raw_qc <- run_raw_qc_for_rgset(config, rgset)
  build_qc_tables(config, chunk, manifest, raw_qc, ewas_qc)

  utils::write.csv(raw_qc$control_pca$scores, study_file(config, "control_probe_pca_scores", chunk = chunk, ext = "csv"), row.names = FALSE)
  utils::write.csv(raw_qc$control_pca$variance, study_file(config, "control_probe_pca_variance", chunk = chunk, ext = "csv"), row.names = FALSE)
  saveRDS(raw_qc$control_pca$pca, study_file(config, "control_probe_pca", chunk = chunk, ext = "rds"))

  try(plot_qc_upset(config, chunk), silent = TRUE)
  message("Wrote QC tables for chunk ", chunk)
  rm(rgset, raw_qc)
  gc()
}
