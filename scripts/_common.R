file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(file_arg)) {
  stop("Could not determine script path. Run pipeline steps with Rscript.", call. = FALSE)
}
pipeline_root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
source(file.path(pipeline_root, "R", "bootstrap.R"))
source_pipeline_R(pipeline_root)

config <- load_config_from_cli(pipeline_root)
ensure_output_dirs(config)

message("Study: ", safe_study_name(config))
message("Config: ", config$.config_file)
