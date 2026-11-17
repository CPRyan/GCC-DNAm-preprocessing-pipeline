derived_variables_dir <- function(config) {
  dir <- file.path(config$paths$output_data, "Derived_Variables")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir
}

derived_file <- function(config, name, suffix = "", ext = "csv") {
  file.path(derived_variables_dir(config), paste0(study_prefix(config), name, suffix, ".", ext))
}

# Row-binds tables whose columns can differ (a missing column is filled with NA)
bind_chunk_tables <- function(tables) {
  cols <- unique(unlist(lapply(tables, colnames)))
  do.call(rbind, lapply(tables, function(t) {
    for (m in setdiff(cols, colnames(t))) t[[m]] <- NA
    t[, cols, drop = FALSE]
  }))
}

# The sample QC tables of all chunks, one file each
export_qc_tables <- function(config, chunks, suffix = "") {
  decisions <- bind_chunk_tables(lapply(chunks, function(k) {
    d <- utils::read.csv(
      qc_filtering_decisions_file(config, k), stringsAsFactors = FALSE, check.names = FALSE,
      colClasses = c(sample_id = "character", failed_metrics_detected = "character",
                     failed_metrics_used_for_filtering = "character")
    )
    d$chunk <- k
    d
  }))
  pass_fail <- bind_chunk_tables(lapply(chunks, function(k) {
    d <- utils::read.csv(qc_pass_fail_file(config, k), stringsAsFactors = FALSE, check.names = FALSE,
                         colClasses = c(sample_id = "character"))
    d$chunk <- k
    d
  }))
  sample_and_chunk_first <- function(d) {
    d[, c("sample_id", "chunk", setdiff(colnames(d), c("sample_id", "chunk"))), drop = FALSE]
  }
  file_decisions <- derived_file(config, "sample_qc_filtering_decisions", suffix)
  file_pass_fail <- derived_file(config, "sample_qc_pass_fail", suffix)
  utils::write.csv(add_name_column(config, sample_and_chunk_first(decisions)), file_decisions, row.names = FALSE)
  utils::write.csv(add_name_column(config, sample_and_chunk_first(pass_fail)), file_pass_fail, row.names = FALSE)
  message("Sample QC tables of ", length(chunks), " chunk(s) (", nrow(decisions), " samples) -> ",
          basename(file_decisions), ", ", basename(file_pass_fail))
  invisible(c(file_decisions, file_pass_fail))
}

# The BMIQ beta and M value matrices of all chunks, one RDS file each (probes in rows, samples
# in columns; the columns are named by the IDATid). The chunks are filtered separately and can
# have different probes, so the file has the probes that all chunks share. With one chunk, the
# chunk file is copied. If the merged matrix would need more than max_gb of memory, it is not
# made (the files of the chunks stay in the data folder).
export_bmiq_matrices <- function(config, chunks, suffix = "", max_gb = 6) {
  kinds <- c(betas = "BMIQ_betas", Mvalues = "BMIQ_Mvalues")
  written <- character()
  for (kind in names(kinds)) {
    files <- vapply(chunks, function(k) study_file(config, kinds[[kind]], chunk = k, ext = "rds"), character(1))
    if (!all(file.exists(files))) {
      warning("BMIQ ", kind, ": ", sum(!file.exists(files)), " of ", length(files),
              " chunk file(s) missing; the matrix was not made", call. = FALSE)
      next
    }
    dest <- derived_file(config, kinds[[kind]], suffix, "rds")
    
    if (length(files) == 1) {
      file.copy(files, dest, overwrite = TRUE)
      message("BMIQ ", kind, ": copied the file of chunk ", chunks, " -> ", basename(dest))
      written <- c(written, dest)
      next
    }
    
    merged <- readRDS(files[1])
    n_probes <- nrow(merged)
    est_gb <- nrow(merged) * ncol(merged) * length(files) * 8 / 1e9
    if (est_gb > max_gb) {
      message("BMIQ ", kind, ": not merged, it would need about ", round(est_gb, 1), " GB of memory. ",
              "The files of the chunks are in the data folder.")
      next
    }
    for (f in files[-1]) {
      m <- readRDS(f)
      if (anyDuplicated(c(colnames(merged), colnames(m)))) stop("Duplicate sample IDs across chunks", call. = FALSE)
      n_probes <- c(n_probes, nrow(m))
      shared <- intersect(rownames(merged), rownames(m))
      merged <- cbind(merged[shared, , drop = FALSE], m[shared, , drop = FALSE])
      rm(m)
      gc()
    }
    saveRDS(merged, dest)
    message("BMIQ ", kind, ": ", format(nrow(merged), big.mark = ","), " probes x ", ncol(merged), " samples",
            if (length(unique(n_probes)) > 1) paste0(" (shared by all chunks; the chunks have ",
                                                     paste(format(n_probes, big.mark = ","), collapse = ", "), " probes)"),
            " -> ", basename(dest))
    written <- c(written, dest)
    rm(merged)
    gc()
  }
  invisible(written)
}