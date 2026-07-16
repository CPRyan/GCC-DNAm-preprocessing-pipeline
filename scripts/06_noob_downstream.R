source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

for (chunk in get_chunks_to_run(config)) {
  message("Running Noob downstream modules for chunk ", chunk)
  run_noob_downstream(config, chunk)
  gc()
}
