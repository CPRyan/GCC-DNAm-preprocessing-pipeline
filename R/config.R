`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

parse_named_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  if (!length(args)) {
    return(out)
  }
  
  for (arg in args) {
    if (!grepl("^--", arg)) {
      next
    }
    arg <- sub("^--", "", arg)
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    key <- parts[[1]]
    value <- if (length(parts) > 1) paste(parts[-1], collapse = "=") else TRUE
    out[[key]] <- value
  }
  
  out
}

as_logical_arg <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required package is not installed: ", package, call. = FALSE)
  }
}

# ----- Fixed locations inside the repo (not user-configurable) ----------------
FIXED_PATHS <- list(
  probe_list_dir  = "data/probe_lists",
  PCClock_data    = "data/PCClocks_data.qs2",
  SystemsAge_data = "data/SystemsAge_data.qs2"
)

# Output subfolders created inside paths$output_dir (not user-configurable)
OUTPUT_SUBDIRS <- c(output_data = "data", output_figures = "figures", output_logs = "logs")
# ------------------------------------------------------------------------------

load_config <- function(config_file) {
  require_package("yaml")
  if (!file.exists(config_file)) {
    stop("Config file does not exist: ", config_file, call. = FALSE)
  }
  # readLines.warn = FALSE: a missing newline at the end of the file is harmless,
  # so it should not give a warning
  config <- yaml::read_yaml(config_file, readLines.warn = FALSE)
  
  # Fixed paths
  ignored <- intersect(names(FIXED_PATHS), names(config$paths))
  if (length(ignored)) {
    warning("These paths are fixed and cannot be changed in the config: ",
            paste(ignored, collapse = ", "), call. = FALSE)
  }
  config$paths[names(FIXED_PATHS)] <- FIXED_PATHS
  
  # Output subfolders, derived from output_dir
  if (is.null(config$paths$output_dir)) {
    stop("paths$output_dir must be set in the config", call. = FALSE)
  }
  ignored <- intersect(names(OUTPUT_SUBDIRS), names(config$paths))
  if (length(ignored)) {
    warning("These paths are derived from output_dir and are ignored: ",
            paste(ignored, collapse = ", "), call. = FALSE)
  }
  config$paths[names(OUTPUT_SUBDIRS)] <- as.list(file.path(config$paths$output_dir, OUTPUT_SUBDIRS))
  
  config
}

load_config_from_cli <- function(pipeline_root) {
  cli <- parse_named_args()
  config_file <- cli$config %||% file.path(pipeline_root, "config", "example_config.yml")
  # warnings raised while reading the config happen before the log exists: keep them
  # (the console still shows them) and write them into the log header
  early_warnings <- character()
  config <- withCallingHandlers(
    load_config(config_file),
    warning = function(w) early_warnings <<- c(early_warnings, conditionMessage(w))
  )
  config$.config_file <- normalizePath(config_file, mustWork = TRUE)
  config$.pipeline_root <- pipeline_root
  config$.cli <- cli
  
  # Resolve the output folders to full paths and create them
  for (key in c("output_dir", names(OUTPUT_SUBDIRS))) {
    config$paths[[key]] <- config_path(config, path.expand(config$paths[[key]]))
    dir.create(config$paths[[key]], recursive = TRUE, showWarnings = FALSE)
  }
  
  # Start the log of this run: every script loads its config here, so every script
  # gets a log in <output_dir>/logs (see start_log() below)
  start_log(config, early_warnings)
  
  config
}

config_path <- function(config, ...) {
  path <- file.path(...)
  if (grepl("^/", path)) {
    return(path)
  }
  file.path(config$.pipeline_root, path)
}

get_chunk_argument <- function(config) {
  cli_chunk <- config$.cli$chunk %||% "all"
  if (identical(tolower(as.character(cli_chunk)), "all")) {
    return("all")
  }
  as.integer(cli_chunk)
}

get_chunks_to_run <- function(config) {
  chunk_arg <- get_chunk_argument(config)
  n_chunks <- as.integer(config$chunking$n_chunks %||% 1)
  if (identical(chunk_arg, "all")) {
    return(seq_len(n_chunks))
  }
  if (is.na(chunk_arg) || chunk_arg < 1 || chunk_arg > n_chunks) {
    stop("Invalid chunk. Expected 1-", n_chunks, " or all.", call. = FALSE)
  }
  chunk_arg
}

# Warnings that are known to be harmless. When one of them is written to a log, its note is
# added after it, so nobody wonders about it. Any other warning is logged as it is. To add
# one, add a part of the warning text (matched literally) and its note.
KNOWN_WARNINGS <- c(
  "underflow occurred in 'lgammacor'" =
    "a numerical note from the BMIQ model fit; it does not stop the run and can be ignored, unless one sample has far more of them than the others"
)

# ---------- LOG OF A SCRIPT RUN -----------------------------------------------
# start_log(config) is called at the end of load_config_from_cli(), so every script
# that loads its config through _common.R gets a log without being edited. From then
# on, everything shown on the console (output, messages, warnings, errors) is also
# written to <output_logs>/<script name>.log, e.g. 09_bmiq_downstream.log. Nothing is
# hidden from the console.
# - The name has no date or time. A new run overwrites the log of the previous run.
# - With --chunk=N the name ends in _chunkN, so chunk jobs do not overwrite each other.
# - The log starts with the script, study, config and start time (and any warning
#   raised while the config was read), and ends with the finish time and how long the
#   run took (and a note if it ended with an error).
# - Lines printed before the config is loaded are not in the log.
# - In an interactive session (no script name) nothing is logged.
# Needs R >= 4.0 (globalCallingHandlers); the pipeline already needs it.
start_log <- function(config, early_warnings = character()) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file_arg)) return(invisible(NULL))
  if (isTRUE(getOption("gcc.log_started"))) return(invisible(NULL))   # only once per run
  options(gcc.log_started = TRUE)
  
  script <- sub("[.][Rr]$", "", basename(sub("^--file=", "", file_arg[[1]])))
  chunk  <- config$.cli$chunk
  if (!is.null(chunk) && !identical(tolower(as.character(chunk)), "all")) {
    script <- paste0(script, "_chunk", chunk)
  }
  log_file <- file.path(config$paths$output_logs, paste0(script, ".log"))
  
  # Two connections to the same file, both appending. Standard output goes through the
  # sink (out_con), which also shows it on the console. Messages, warnings and errors are
  # written to the second one (msg_con): the console already shows them, and writing them
  # to the sink connection would print them a second time.
  file.create(log_file)   # start with an empty file
  out_con <- file(log_file, open = "at")
  msg_con <- file(log_file, open = "at")
  started <- Sys.time()
  failed  <- FALSE
  stamp   <- function(t) format(t, "%Y-%m-%d %H:%M:%S")
  write_out <- function(...) {
    cat(..., file = out_con, sep = "")
    flush(out_con)
  }
  to_log <- function(...) {
    flush(out_con)   # earlier standard output first, to keep the order of the lines
    cat(..., file = msg_con, sep = "")
    flush(msg_con)
  }
  
  # The same warning can be raised thousands of times. The first max_repeats are written to
  # the log; the total is written at the end. Known harmless warnings get their note.
  warned      <- new.env()
  max_repeats <- 3L
  note_for <- function(txt) {
    hit <- names(KNOWN_WARNINGS)[vapply(names(KNOWN_WARNINGS),
                                        function(p) grepl(p, txt, fixed = TRUE), logical(1))]
    if (length(hit)) paste0("  [", KNOWN_WARNINGS[[hit[1]]], "]") else ""
  }
  
  write_out("Script:  ", script, "\n",
            "Study:   ", config$study_name %||% "", "\n",
            "Config:  ", config$.config_file, "\n",
            "Started: ", stamp(started), "\n",
            "R:       ", R.version.string, "\n")
  for (w in early_warnings) write_out("Warning (while reading the config): ", w, "\n")
  write_out("\n")
  
  sink(out_con, split = TRUE)   # standard output (print, cat): console and log
  
  # messages, warnings and errors are copied to the log; the console still shows them
  globalCallingHandlers(
    message = function(m) {
      # package start-up text (attaching minfi, lists of masked functions, ...) and two
      # harmless notes stay on the console only
      txt <- conditionMessage(m)
      if (inherits(m, "packageStartupMessage") ||
          grepl("^(Setting options\\(|No methods found in package)", txt)) {
        return(invisible(NULL))
      }
      to_log(txt)
    },
    warning = function(w) {
      txt <- conditionMessage(w)
      n   <- (warned[[txt]] %||% 0L) + 1L
      assign(txt, n, envir = warned)
      if (n <= max_repeats) to_log("Warning: ", txt, note_for(txt), "\n")
    },
    error   = function(e) {
      failed <<- TRUE
      to_log("Error: ", conditionMessage(e), "\n")
    }
  )
  
  # when R quits (also after an error): finish and close the log
  reg.finalizer(globalenv(), function(e) {
    tryCatch({
      sink()
      flush(out_con)
      for (txt in sort(ls(warned, all.names = TRUE))) {
        n <- get(txt, envir = warned)
        if (n > max_repeats) {
          cat("Warning (repeated, ", format(n, big.mark = ","), " times in total; the first ", max_repeats,
              " are written above): ", txt, note_for(txt), "\n", file = msg_con, sep = "")
        }
      }
      cat("\nFinished: ", stamp(Sys.time()), " (",
          round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), " min)",
          if (failed) " - ended with an error", "\n", file = msg_con, sep = "")
      flush(msg_con)
      close(out_con)
      close(msg_con)
    }, error = function(err) NULL)
  }, onexit = TRUE)
  
  invisible(log_file)
}