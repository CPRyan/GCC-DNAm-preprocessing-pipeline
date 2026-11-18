source(file.path(
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])),
  "_common.R"
))

# Runs the pipeline steps (scripts/NN_*.R) one after the other with one command:
#
#   Rscript scripts/run_all_steps.R --config=config/my_study.yml
#
# Options:
#   --from=05 --to=07    run only the steps numbered 05 to 07
#
# Every step runs in its own R process, exactly as if you had started it yourself, so
# each one writes its own log (<output_dir>/logs/<step name>.log) and frees its memory
# before the next one starts. This script has its own log, run_all_steps.log, with the
# start and end of each step and a summary. It stops at the first step that fails.
# All chunks are always run (single chunks are for running one step by hand).

scripts_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]]))

n_chunks <- as.integer(config$chunking$n_chunks %||% 1)
if (n_chunks > 1 && !identical(get_chunk_argument(config), "all")) {
  stop("run_all_steps.R always runs all chunks; remove --chunk.", call. = FALSE)
}

# the steps are the scripts that start with a two-digit number, in that order
steps   <- sort(list.files(scripts_dir, pattern = "^[0-9]{2}_.*[.]R$"))
step_no <- as.integer(substr(steps, 1, 2))
from    <- as.integer(config$.cli$from %||% min(step_no))
to      <- as.integer(config$.cli$to   %||% max(step_no))
steps   <- steps[step_no >= from & step_no <= to]
if (!length(steps)) stop("No steps found between ", from, " and ", to, call. = FALSE)

message("Steps to run (", length(steps), "): ", paste(sub("[.]R$", "", steps), collapse = ", "))

# the start of the run and its steps, for the report (step 09 reads them from the environment)
Sys.setenv(GCC_RUN_START = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
           GCC_RUN_STEPS = paste(substr(steps[1], 1, 2), "to", substr(steps[length(steps)], 1, 2)))

rscript <- file.path(R.home("bin"), "Rscript")
results <- data.frame(step = sub("[.]R$", "", steps), status = NA_integer_, minutes = NA_real_,
                      warnings = NA_integer_, failed_modules = NA_integer_,
                      stringsAsFactors = FALSE)
failed_modules <- character()

for (i in seq_along(steps)) {
  step <- results$step[i]
  message("\n======= Step ", step, " =======")
  t0 <- Sys.time()
  status <- system2(rscript, c(shQuote(file.path(scripts_dir, steps[i])),
                               shQuote(paste0("--config=", config$.config_file))))
  results$status[i]  <- status
  results$minutes[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  
  # modules that fail inside a step only produce a warning, so look for them in the step's log
  log_file <- file.path(config$paths$output_logs, paste0(step, ".log"))
  lines    <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character()
  fails    <- grep("^======= .* FAILED", lines, value = TRUE)
  results$warnings[i]       <- sum(grepl("^Warning: ", lines))
  results$failed_modules[i] <- length(fails)
  failed_modules <- c(failed_modules, fails)
  
  if (status != 0) {
    message("======= ", step, " FAILED (exit status ", status, "); see ", log_file, " =======")
    break
  }
  message("======= ", step, " completed (", results$minutes[i], " min)",
          if (length(fails)) paste0(", but ", length(fails), " module(s) FAILED") else "", " =======")
}

message("\nSummary:")
print(results[!is.na(results$status), ], row.names = FALSE)

if (any(results$status != 0, na.rm = TRUE)) {
  stop("The pipeline stopped at step ", results$step[which(results$status != 0)[1]], ".", call. = FALSE)
}
if (length(failed_modules)) {
  warning("The pipeline finished, but some modules failed:\n", paste(failed_modules, collapse = "\n"),
          call. = FALSE)
} else {
  message("The pipeline finished without errors.")
}