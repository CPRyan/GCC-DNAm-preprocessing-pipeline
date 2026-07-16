source_pipeline_R <- function(pipeline_root) {
  r_dir <- file.path(pipeline_root, "R")
  files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  files <- files[basename(files) != "bootstrap.R"]
  files <- files[order(basename(files))]

  for (file in files) {
    source(file)
  }

  invisible(files)
}
