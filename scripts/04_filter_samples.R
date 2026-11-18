source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

chunks <- get_chunks_to_run(config)
for (chunk in chunks) {
  message("Filtering samples for chunk ", chunk)
  out <- save_good_sample_rgset(config, chunk)
  message("Wrote: ", out)
  gc()
}

# which samples were removed and why, over all chunks
sample_filtering_summary(config, chunks)