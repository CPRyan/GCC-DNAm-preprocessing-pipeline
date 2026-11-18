qc_value_failed <- function(x) {
  if (is.logical(x)) {
    return(is.na(x) | !x)
  }
  lx <- tolower(as.character(x))
  is.na(x) | lx %in% c("fail", "failed", "false", "f", "0", "no", "n")
}

# The QC columns that are used for removing samples, as set in the config. A metric that
# is switched on in the config but has no column in the QC table cannot be used: it is
# reported with a warning instead of being skipped silently.
selected_filter_columns <- function(config, pass_fail) {
  # Bisulfite conversion has no entry here: it is checked by the three "Bisulfite Conversion"
  # metrics of EWAS tools (see sample_filtering.ewas_tools_metrics in the config).
  standard_map <- c(
    signal_quality = "passed_signal_quality",
    beadcount = "passed_beadcount",
    detection_p = "passed_detection_cutoff",
    snp_outlier = "passed_snp_outlier_bl4",
    sex_mismatch = "passed_sex_mismatch"
  )
  
  selected <- character()
  ignored  <- character()   # switched on in the config, but not in the QC table
  
  standard <- config$sample_filtering$standard_metrics %||% list()
  for (metric in names(standard)) {
    if (!isTRUE(standard[[metric]])) next
    col <- unname(standard_map[metric])   # NA for a name that is not known
    if (!is.na(col) && col %in% colnames(pass_fail)) {
      selected <- c(selected, col)
    } else {
      ignored <- c(ignored, metric)
    }
  }
  
  ewas <- config$sample_filtering$ewas_tools_metrics %||% list()
  for (metric in names(ewas)) {
    if (!isTRUE(ewas[[metric]])) next
    if (metric %in% colnames(pass_fail)) {
      selected <- c(selected, metric)
    } else {
      ignored <- c(ignored, metric)
    }
  }
  
  if (length(ignored)) {
    warning("Switched on in the config but not available in the QC table, so NOT used for ",
            "filtering samples: ", paste(ignored, collapse = ", "), call. = FALSE)
  }
  
  unique(selected)
}

# One list entry per row: the names of the columns that failed (character(0) if none).
# Tests each column in its own type, and always returns one entry per sample.
failed_names <- function(df) {
  failed <- lapply(df, qc_value_failed)
  lapply(seq_len(nrow(df)), function(i) {
    colnames(df)[vapply(failed, function(f) f[[i]], logical(1))]
  })
}

apply_sample_filtering <- function(config, chunk) {
  pass_fail <- utils::read.csv(qc_pass_fail_file(config, chunk), stringsAsFactors = FALSE, check.names = FALSE)
  selected_cols <- selected_filter_columns(config, pass_fail)
  if (!"sample_id" %in% colnames(pass_fail)) {
    stop("Pass/fail table is missing sample_id.", call. = FALSE)
  }
  
  detected_cols <- setdiff(colnames(pass_fail), c("sample_id", "Basename", "chunk"))
  failed_detected <- failed_names(pass_fail[, detected_cols, drop = FALSE])
  failed_used     <- failed_names(pass_fail[, selected_cols, drop = FALSE])
  
  decisions <- data.frame(
    sample_id = pass_fail$sample_id,
    failed_metrics_detected = vapply(failed_detected, paste, character(1), collapse = ";"),
    failed_metrics_used_for_filtering = vapply(failed_used, paste, character(1), collapse = ";"),
    removed_from_pipeline = lengths(failed_used) > 0,
    stringsAsFactors = FALSE
  )
  utils::write.csv(decisions, qc_filtering_decisions_file(config, chunk), row.names = FALSE)
  
  # what was used and what happened, for the log (a missing value counts as a failure)
  n_fail   <- vapply(pass_fail[, detected_cols, drop = FALSE],
                     function(x) sum(qc_value_failed(x)), numeric(1))
  failing  <- n_fail[n_fail > 0]
  message("Sample filtering (chunk ", chunk, "): ", nrow(decisions), " samples, ",
          sum(decisions$removed_from_pipeline), " removed, ",
          sum(!decisions$removed_from_pipeline), " kept")
  if (length(failing)) {
    message("  samples failing each metric: ",
            paste0(names(failing), " ", failing,
                   ifelse(names(failing) %in% selected_cols, " (used)", " (not used)"),
                   collapse = "; "))
  } else {
    message("  no sample fails any metric")
  }
  
  # ---- the identity and control checks, one line each ----
  removed <- decisions$removed_from_pipeline
  report_check <- function(label, pattern) {
    cols <- grep(pattern, detected_cols, value = TRUE)
    if (!length(cols)) {
      message("  ", label, ": not available (this check was not run)")
      return(invisible(NULL))
    }
    failed <- vapply(failed_detected, function(f) any(f %in% cols), logical(1))
    if (!any(failed)) {
      message("  ", label, ": 0 sample(s) failed")
      return(invisible(failed))
    }
    if (any(cols %in% selected_cols)) {
      message("  ", label, ": ", sum(failed), " sample(s) failed -> removed")
    } else {
      message("  ", label, ": ", sum(failed), " sample(s) failed -> NOT removed (not used for filtering); ",
              sum(failed & !removed), " of them are still in the sample")
    }
    invisible(failed)
  }
  message("Sample checks (chunk ", chunk, "):")
  report_check("bisulfite conversion (EWAS tools)", "^Bisulfite Conversion")
  report_check("SNP outliers", "^passed_snp_outlier")
  report_check("sex mismatch", "^passed_sex_mismatch")
  metrics_file <- qc_metrics_file(config, chunk)
  if (file.exists(metrics_file)) {
    m <- utils::read.csv(metrics_file, stringsAsFactors = FALSE, check.names = FALSE)
    if ("sex_check" %in% colnames(m)) {
      tab <- table(m$sex_check)
      message("  sex check results: ", paste(names(tab), tab, sep = " = ", collapse = ", "))
    }
  }
  
  decisions
}

save_good_sample_rgset <- function(config, chunk) {
  require_package("minfi")
  decisions <- apply_sample_filtering(config, chunk)
  rgset <- readRDS(raw_rgset_file(config, chunk))
  keep <- decisions$sample_id[!decisions$removed_from_pipeline]
  saveRDS(rgset[, colnames(rgset) %in% keep], good_rgset_file(config, chunk))
  invisible(good_rgset_file(config, chunk))
}

# Summary of step 04 for the log, over the chunks that were run: how many samples were
# removed and why, with their IDs. Samples that failed a metric that is NOT used for
# filtering are listed as "flagged but kept". At most max_listed samples are listed in
# each group; the filtering decisions tables have all of them.
# Returns the removed samples invisibly.
# The sample sheet column with a name or ID of the sample, shown next to the IDATid in the
# logs and the report (sample_sheet: name_col in the config; default Sample_Name)
name_column <- function(config) {
  config$sample_sheet$name_col %||% "Sample_Name"
}

# sample_id -> name, from the sample sheet (NULL if the sheet has no such column)
sample_name_lookup <- function(config) {
  manifest <- tryCatch(read_chunk_assignments(config), error = function(e) NULL)
  col <- name_column(config)
  if (is.null(manifest) || !col %in% colnames(manifest)) return(NULL)
  stats::setNames(as.character(manifest[[col]]), manifest$sample_id)
}

# For the output tables: the sample name column (Sample_Name) first, then the IDATid, which is
# written under the name of the sample ID column of the sample sheet instead of "sample_id"
add_name_column <- function(config, tbl) {
  if (!"sample_id" %in% colnames(tbl)) return(tbl)
  id_col <- config$sample_sheet$sample_id_col %||% "IDATid"
  col    <- name_column(config)
  lookup <- sample_name_lookup(config)
  use_name <- !is.null(lookup) && !identical(col, id_col)
  if (use_name) tbl[[col]] <- unname(lookup[as.character(tbl$sample_id)])
  colnames(tbl)[colnames(tbl) == "sample_id"] <- id_col
  first <- c(if (use_name) col, id_col)
  tbl[, c(first, setdiff(colnames(tbl), first)), drop = FALSE]
}

# Readable names for the QC columns, for the reasons of removed and flagged samples
reason_text <- function(x) {
  labels <- c(passed_signal_quality = "signal intensity", passed_beadcount = "bead count",
              passed_detection_cutoff = "detection p-value", passed_snp_outlier_bl4 = "SNP outlier",
              passed_sex_mismatch = "sex mismatch")
  vapply(strsplit(x, ";", fixed = TRUE), function(p) {
    p <- ifelse(p %in% names(labels), unname(labels[p]), p)
    paste(p, collapse = ", ")
  }, character(1))
}

# "reported female, predicted male" for every sample that has a predicted sex (named by sample ID)
sex_detail_from_metrics <- function(config, chunks) {
  out <- character()
  lab <- function(x) ifelse(x %in% "f", "female", ifelse(x %in% "m", "male", "unknown"))
  for (k in chunks) {
    file <- qc_metrics_file(config, k)
    if (!file.exists(file)) next
    m <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c("sample_id", "reported_sex", "predicted_sex") %in% colnames(m))) {
      out[as.character(m$sample_id)] <- paste0("reported ", lab(m$reported_sex), ", predicted ", lab(m$predicted_sex))
    }
  }
  out
}

# The reasons of samples (raw = the ";"-separated QC columns), with the reported and predicted
# sex added to "sex mismatch"
reason_with_sex <- function(raw, ids, sex_detail) {
  why <- reason_text(raw)
  for (j in which(grepl("passed_sex_mismatch", raw, fixed = TRUE) & ids %in% names(sex_detail))) {
    why[j] <- sub("sex mismatch", paste0("sex mismatch (", sex_detail[[ids[j]]], ")"), why[j], fixed = TRUE)
  }
  why
}

sample_filtering_summary <- function(config, chunks, max_listed = 100) {
  read_decisions <- function(k) {
    d <- utils::read.csv(
      qc_filtering_decisions_file(config, k), stringsAsFactors = FALSE, check.names = FALSE,
      colClasses = c(sample_id = "character", failed_metrics_detected = "character",
                     failed_metrics_used_for_filtering = "character")
    )
    d$chunk <- k
    d
  }
  decisions <- do.call(rbind, lapply(chunks, read_decisions))
  decisions$failed_metrics_detected[is.na(decisions$failed_metrics_detected)] <- ""
  decisions$failed_metrics_used_for_filtering[is.na(decisions$failed_metrics_used_for_filtering)] <- ""
  
  # the name or ID next to the IDAT ID, when the sample sheet has that column
  manifest <- tryCatch(read_chunk_assignments(config), error = function(e) NULL)
  name_col <- name_column(config)
  name_of <- function(ids) {
    if (!is.null(manifest) && name_col %in% colnames(manifest)) {
      paste0(ids, " (", name_col, " ", manifest[[name_col]][match(ids, manifest$sample_id)], ")")
    } else {
      ids
    }
  }
  
  removed <- decisions$removed_from_pipeline
  message("Sample filtering summary (chunk(s) ", paste(chunks, collapse = ", "), "): ",
          nrow(decisions), " samples, ", sum(removed), " removed, ", sum(!removed), " kept")
  
  # which metrics were used (the same in every chunk, so they are listed once, here)
  pass_fail <- utils::read.csv(qc_pass_fail_file(config, chunks[1]), stringsAsFactors = FALSE, check.names = FALSE)
  used      <- suppressWarnings(selected_filter_columns(config, pass_fail))
  not_used  <- setdiff(setdiff(colnames(pass_fail), c("sample_id", "Basename", "chunk")), used)
  message("  metrics used for filtering (", length(used), "): ",
          if (length(used)) paste(used, collapse = ", ") else "none")
  if (length(not_used)) {
    message("  metrics collected but not used for filtering: ", paste(not_used, collapse = ", "))
  }
  
  if (any(removed)) {
    why <- unlist(strsplit(decisions$failed_metrics_used_for_filtering[removed], ";", fixed = TRUE))
    tab <- sort(table(why), decreasing = TRUE)
    message("  removed, by failed metric (a sample can fail more than one): ",
            paste(reason_text(names(tab)), tab, sep = " = ", collapse = "; "))
  }
  
  sex_detail <- sex_detail_from_metrics(config, chunks)
  
  list_samples <- function(title, rows, metrics) {
    ids <- decisions$sample_id[rows]
    if (!length(ids)) return(invisible(NULL))
    message("  ", title, " (", length(ids), "):")
    why <- reason_with_sex(metrics[rows], ids, sex_detail)
    for (i in utils::head(seq_along(ids), max_listed)) {
      message("    ", name_of(ids[i]), " [chunk ", decisions$chunk[rows][i], "]: ", why[i])
    }
    if (length(ids) > max_listed) {
      message("    ... and ", length(ids) - max_listed, " more (see the filtering decisions tables)")
    }
  }
  list_samples("removed samples", removed, decisions$failed_metrics_used_for_filtering)
  list_samples("flagged but kept (failed a check that is not used for filtering, for example a sex mismatch)",
               !removed & nzchar(decisions$failed_metrics_detected),
               decisions$failed_metrics_detected)
  
  invisible(decisions[removed, , drop = FALSE])
}