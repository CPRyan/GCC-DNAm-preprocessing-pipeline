read_illumina_sample_sheet <- function(sample_sheet) {
  lines <- readLines(sample_sheet, warn = FALSE)
  data_line <- grep("^\\[Data\\]", trimws(lines))

  if (length(data_line)) {
    csv_text <- paste(lines[(data_line[[1]] + 1):length(lines)], collapse = "\n")
    sample_sheet_data <- utils::read.csv(
      text = csv_text,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = "character",
      na.strings = c("", "NA")
    )
  } else {
    sample_sheet_data <- utils::read.csv(
      sample_sheet,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = "character",
      na.strings = c("", "NA")
    )
  }

  sample_sheet_data <- sample_sheet_data[rowSums(!is.na(sample_sheet_data)) > 0, , drop = FALSE]
  rownames(sample_sheet_data) <- NULL
  sample_sheet_data
}

list_idat_basenames <- function(idat_dir, recursive = TRUE) {
  red_files <- list.files(idat_dir, pattern = "_Red[.]idat$", full.names = TRUE, recursive = recursive)
  green_files <- list.files(idat_dir, pattern = "_Grn[.]idat$", full.names = TRUE, recursive = recursive)
  red_base <- sub("_Red[.]idat$", "", red_files)
  green_base <- sub("_Grn[.]idat$", "", green_files)
  paired <- intersect(red_base, green_base)
  stats::setNames(paired, basename(paired))
}

add_idat_basenames <- function(sample_sheet_data,
                               idat_dir,
                               basename_col = NULL,
                               sentrix_id_col = "Sentrix_ID",
                               sentrix_position_col = "Sentrix_Position",
                               recursive = TRUE) {
  idat_map <- list_idat_basenames(idat_dir, recursive = recursive)
  if (!length(idat_map)) {
    stop("No paired Red/Grn IDAT files were found in: ", idat_dir, call. = FALSE)
  }

  if (!is.null(basename_col) && nzchar(basename_col)) {
    if (!basename_col %in% colnames(sample_sheet_data)) {
      stop("Sample sheet is missing basename column: ", basename_col, call. = FALSE)
    }
    basename_key <- as.character(sample_sheet_data[[basename_col]])
  } else {
    required_cols <- c(sentrix_id_col, sentrix_position_col)
    missing_cols <- setdiff(required_cols, colnames(sample_sheet_data))
    if (length(missing_cols)) {
      stop("Sample sheet is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    basename_key <- paste0(sample_sheet_data[[sentrix_id_col]], "_", sample_sheet_data[[sentrix_position_col]])
  }

  basename_key <- sub("_Red[.]idat$", "", basename_key)
  basename_key <- sub("_Grn[.]idat$", "", basename_key)
  basename_key <- basename(basename_key)
  sample_sheet_data$Basename <- unname(idat_map[basename_key])
  sample_sheet_data
}

validate_idat_pairs <- function(sample_sheet_data) {
  if (!"Basename" %in% colnames(sample_sheet_data)) {
    stop("Sample sheet data must contain a Basename column.", call. = FALSE)
  }

  sample_sheet_data$Red_IDAT <- paste0(sample_sheet_data$Basename, "_Red.idat")
  sample_sheet_data$Green_IDAT <- paste0(sample_sheet_data$Basename, "_Grn.idat")
  sample_sheet_data$red_idat_exists <- file.exists(sample_sheet_data$Red_IDAT)
  sample_sheet_data$green_idat_exists <- file.exists(sample_sheet_data$Green_IDAT)
  sample_sheet_data$idat_pair_exists <- sample_sheet_data$red_idat_exists & sample_sheet_data$green_idat_exists
  sample_sheet_data
}

load_and_validate_manifest <- function(config) {
  sample_sheet_file <- config_path(config, config$paths$sample_sheet)
  idat_dir <- config_path(config, config$paths$idat_dir)
  ss <- config$sample_sheet

  manifest <- read_illumina_sample_sheet(sample_sheet_file)
  manifest <- add_idat_basenames(
    manifest,
    idat_dir = idat_dir,
    basename_col = ss$basename_col,
    sentrix_id_col = ss$sentrix_id_col %||% "Sentrix_ID",
    sentrix_position_col = ss$sentrix_position_col %||% "Sentrix_Position",
    recursive = isTRUE(ss$recursive_idat_search)
  )
  manifest <- validate_idat_pairs(manifest)

  if (!ss$sample_id_col %in% colnames(manifest)) {
    stop("Sample sheet is missing sample ID column: ", ss$sample_id_col, call. = FALSE)
  }

  manifest$sample_id <- make.unique(as.character(manifest[[ss$sample_id_col]]))
  manifest
}
