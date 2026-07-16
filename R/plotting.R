plot_qc_upset <- function(config, chunk) {
  require_package("UpSetR")
  pass_fail <- utils::read.csv(qc_pass_fail_file(config, chunk), stringsAsFactors = FALSE, check.names = FALSE)
  metric_cols <- setdiff(colnames(pass_fail), c("sample_id", "Basename", "chunk"))
  upset_list <- list()
  for (col in metric_cols) {
    failed <- qc_value_failed(pass_fail[[col]])
    if (any(failed, na.rm = TRUE)) {
      upset_list[[col]] <- pass_fail$sample_id[failed]
    }
  }
  if (!length(upset_list)) {
    return(invisible(NULL))
  }
  file <- study_file(config, "sample_failure_upset", chunk = chunk, ext = "pdf", dir = output_figure_dir(config))
  grDevices::pdf(file, height = 7, width = 13)
  print(UpSetR::upset(UpSetR::fromList(upset_list), set_size.show = TRUE))
  grDevices::dev.off()
  invisible(file)
}
