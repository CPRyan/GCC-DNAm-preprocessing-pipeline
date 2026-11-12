source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

# ---- Nothing to merge with a single chunk 
n_chunks <- as.integer(config$chunking$n_chunks %||% 1)
if (n_chunks <= 1) {
  message("Only ", n_chunks, " chunk configured: nothing to merge. ",
          "The chunk-level outputs (chunk1) are the final outputs.")
  quit(save = "no", status = 0)
}

# ---- Merge (2 or more chunks) 
outs <- merge_standard_outputs(config)
message("Merged standard chunk-level CSV outputs.")
print(outs)

clocks_outs <- merge_clock_outputs(config)
message("Merged clocks chunk-level CSV outputs.")
print(clocks_outs)

cell_counts_outs <- merge_cell_type_outputs(config)
message("Merged cell count chunk-level CSV outputs.")
print(cell_counts_outs)

bmiq_out <- merge_bmiq_betas(config)
message("Merged BMIQ beta chunk-level RDS outputs.")
print(bmiq_out)