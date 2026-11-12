merge_chunk_csvs <- function(config, suffix, output_suffix = suffix) {
  files <- vapply(seq_len(as.integer(config$chunking$n_chunks %||% 1)), function(chunk) {
    study_file(config, suffix, chunk = chunk, ext = "csv")
  }, character(1))
  files <- files[file.exists(files)]
  if (!length(files)) {
    warning("No files found for suffix: ", suffix, call. = FALSE)
    return(invisible(NULL))
  }
  merged <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE, check.names = FALSE))
  if ("sample_id" %in% colnames(merged) && anyDuplicated(merged$sample_id)) {
    warning("Duplicate sample_id rows found when merging suffix: ", suffix, call. = FALSE)
  }
  out <- study_file(config, paste0(output_suffix, "_all_chunks"), ext = "csv")
  utils::write.csv(merged, out, row.names = FALSE)
  invisible(out)
}

merge_standard_outputs <- function(config) {
  list(
    qc_metrics = merge_chunk_csvs(config, "sample_qc_metrics"),
    qc_pass_fail = merge_chunk_csvs(config, "sample_qc_pass_fail"),
    qc_filtering_decisions = merge_chunk_csvs(config, "sample_qc_filtering_decisions")
  )
}

merge_clock_outputs <- function(config) {
  dunedinPACE_file <- merge_chunk_csvs(config, "DunedinPACE_noob")
  pc_clocks_file <- merge_chunk_csvs(config, "PCClocks_noob")
  
  dunedinPACE <- if (!is.null(dunedinPACE_file)) {
    df <- utils::read.csv(dunedinPACE_file, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
    df$sample_id <- rownames(df)
    rownames(df) <- NULL
    df
  } else NULL
  
  pc_clocks <- if (!is.null(pc_clocks_file)) {
    utils::read.csv(pc_clocks_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else NULL
  
  combined <- merge_by_sample(dunedinPACE, pc_clocks)
  combined$Var.3 <- NULL # remove the duplicated col
  if (!nrow(combined)) {
    warning("No clock data found to combine.", call. = FALSE)
    return(invisible(NULL))
  }
  if (anyDuplicated(combined$sample_id)) {
    warning("Duplicate sample_id rows found when combining clock outputs.", call. = FALSE)
  }
  
  out <- study_file(config, "clocks_combined_all_chunks", ext = "csv")
  utils::write.csv(combined, out, row.names = FALSE)
  invisible(out)
}
  
read_chunks_rownames_to_sample_id <- function(config, suffix) {
  files <- vapply(seq_len(as.integer(config$chunking$n_chunks %||% 1)), function(chunk) {
    study_file(config, suffix, chunk = chunk, ext = "csv")
  }, character(1))
  files <- files[file.exists(files)]
  if (!length(files)) {
    warning("No files found for suffix: ", suffix, call. = FALSE)
    return(NULL)
  }
  tables <- lapply(files, function(f) {
    df <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
    df$sample_id <- rownames(df)
    rownames(df) <- NULL
    df
  })
  do.call(rbind, tables)
}

add_suffix <- function(df, suffix) {
  if (is.null(df)) return(NULL)
  cols_to_suffix <- setdiff(colnames(df), "sample_id")
  colnames(df)[colnames(df) %in% cols_to_suffix] <- paste0(cols_to_suffix, "_", suffix)
  df
}

merge_cell_type_outputs <- function(config) {
  blood_extended <- read_chunks_rownames_to_sample_id(config, "BloodExtended_cell_counts")
  hepidish <- read_chunks_rownames_to_sample_id(config, "HEpiDISH_RPC_cell_counts")
  ewastools <- read_chunks_rownames_to_sample_id(config, "Ewastools_cell_counts")
  
  blood_extended <- add_suffix(blood_extended, "BE")
  hepidish <- add_suffix(hepidish, "HD")
  ewastools <- add_suffix(ewastools, "ET")
  
  combined <- merge_by_sample(blood_extended, hepidish, ewastools)
  if (!nrow(combined)) {
    warning("No cell count data found to combine.", call. = FALSE)
    return(invisible(NULL))
  }
  if (anyDuplicated(combined$sample_id)) {
    warning("Duplicate sample_id rows found when combining cell count outputs.", call. = FALSE)
  }
  
  out <- study_file(config, "cell_counts_combined_all_chunks", ext = "csv")
  utils::write.csv(combined, out, row.names = FALSE)
  invisible(out)
}

merge_bmiq_betas <- function(config, suffix = "BMIQ_betas", probes = c("intersect", "union")) {
  probes <- match.arg(probes)
  n_chunks <- as.integer(config$chunking$n_chunks %||% 1)
  files <- vapply(seq_len(n_chunks), function(chunk) {
    study_file(config, suffix, chunk = chunk, ext = "rds")
  }, character(1))
  
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) {
    stop("Missing BMIQ beta file(s): ", paste(basename(missing_files), collapse = ", "), call. = FALSE)
  }
  
  betas <- lapply(files, function(f) {
    b <- as.matrix(readRDS(f))
    if (is.null(rownames(b)) || is.null(colnames(b))) {
      stop("BMIQ beta matrix has no row/column names: ", basename(f), call. = FALSE)
    }
    b
  })
  
  sample_ids <- unlist(lapply(betas, colnames))
  if (anyDuplicated(sample_ids)) {
    stop("Duplicate sample IDs across BMIQ chunks: ",
         paste(head(unique(sample_ids[duplicated(sample_ids)]), 5), collapse = ", "), call. = FALSE)
  }
  
  probe_sets <- lapply(betas, rownames)
  probe_ids <- if (probes == "intersect") Reduce(intersect, probe_sets) else Reduce(union, probe_sets)
  for (i in seq_along(betas)) {
    message("  chunk ", i, ": ", nrow(betas[[i]]), " probes x ", ncol(betas[[i]]), " samples")
  }
  message("  ", length(probe_ids), " probes kept (", probes, " across chunks)")
  if (!length(probe_ids)) stop("No probes left after combining chunks.", call. = FALSE)
  
  # match() puts every chunk in the same probe order (and gives NA rows for union)
  merged <- do.call(cbind, lapply(betas, function(b) {
    m <- b[match(probe_ids, rownames(b)), , drop = FALSE]
    rownames(m) <- probe_ids
    m
  }))
  rm(betas); gc()
  
  out <- study_file(config, paste0(suffix, "_all_chunks"), ext = "rds")
  saveRDS(merged, out)
  message("  wrote ", nrow(merged), " probes x ", ncol(merged), " samples to ", out)
  invisible(out)
}