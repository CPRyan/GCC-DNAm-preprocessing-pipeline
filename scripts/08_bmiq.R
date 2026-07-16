source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

message("BMIQ is a standard but explicit step. Running because this script was called.")
for (chunk in get_chunks_to_run(config)) {
  message("Running BMIQ for chunk ", chunk)
  out <- bmiq_chunk(config, chunk)
  message("Wrote: ", out)
  gc()
}
