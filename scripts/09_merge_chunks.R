source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

outs <- merge_standard_outputs(config)
message("Merged standard chunk-level CSV outputs.")
print(outs)
