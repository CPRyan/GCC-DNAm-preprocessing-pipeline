source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

for (chunk in get_chunks_to_run(config)) {
  message("Probe filtering chunk ", chunk)
  out <- probe_filter_chunk(config, chunk)
  message("Wrote: ", out)
  gc()
}
