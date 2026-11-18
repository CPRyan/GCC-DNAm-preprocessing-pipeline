source(file.path(
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])),
  "_common.R"
))

# BMIQ downstream: M values and FraminghamPACE per chunk, then merge and plots
# (R/bmiq_downstream.R, loaded by _common.R).

run_bmiq_downstream(config)