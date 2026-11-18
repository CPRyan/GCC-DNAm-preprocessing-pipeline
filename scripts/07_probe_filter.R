source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

# The probes are filtered chunk by chunk; --chunk=N runs a single chunk.
run_probe_filtering(config)