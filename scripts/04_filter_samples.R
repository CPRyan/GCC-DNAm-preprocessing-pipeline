source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

for (chunk in get_chunks_to_run(config)) {
  message("Filtering samples for chunk ", chunk)
  out <- save_good_sample_rgset(config, chunk)
  message("Wrote: ", out)
  gc()
}
