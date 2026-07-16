source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

step <- config$.cli$step %||% stop("Provide --step=validate|targets|load_rgset|raw_qc|filter_samples|noob|noob_downstream|probe_filter|bmiq|merge", call. = FALSE)
script_map <- c(
  validate = "00_validate_inputs.R",
  targets = "01_create_targets_and_chunks.R",
  load_rgset = "02_load_rgset_chunks.R",
  raw_qc = "03_run_raw_qc.R",
  filter_samples = "04_filter_samples.R",
  noob = "05_noob_normalize.R",
  noob_downstream = "06_noob_downstream.R",
  probe_filter = "07_probe_filter.R",
  bmiq = "08_bmiq.R",
  merge = "09_merge_chunks.R"
)

if (!step %in% names(script_map)) {
  stop("Unknown step: ", step, call. = FALSE)
}

script <- file.path(pipeline_root, "scripts", script_map[[step]])
args <- commandArgs(trailingOnly = TRUE)
args <- args[!grepl("^--step=", args)]
status <- system2("Rscript", c(script, args))
quit(status = status)
