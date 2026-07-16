source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

for (chunk in get_chunks_to_run(config)) {
  message("Loading RGset for chunk ", chunk)
  out <- load_and_save_rgset_chunk(config, chunk)
  message("Wrote: ", out)
  gc()
}
