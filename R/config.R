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

load_config <- function(config_file) {
  require_package("yaml")
  if (!file.exists(config_file)) {
    stop("Config file does not exist: ", config_file, call. = FALSE)
  }
  yaml::read_yaml(config_file)
}

get_script_pipeline_root <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

load_config_from_cli <- function(pipeline_root) {
  cli <- parse_named_args()
  config_file <- cli$config %||% file.path(pipeline_root, "config", "example_config.yml")
  config <- load_config(config_file)
  config$.config_file <- normalizePath(config_file, mustWork = TRUE)
  config$.pipeline_root <- pipeline_root
  config$.cli <- cli
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
