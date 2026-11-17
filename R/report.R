# ---------- REPORT: THE NUMBERS --------------------------------------------------------
# The report collects the numbers of a run for a methods section: how many samples were
# removed and why, how many probes were removed by each filter, the settings and the
# software versions. Everything is read from the files that the earlier steps wrote (the
# sample sheet with its chunks, the QC pass/fail and filtering decisions tables, the probe
# filtering table) and from the config, so it needs no extra input.
#
# collect_report_numbers() returns the numbers as a list of tables;
# write_methods_summary() writes them as a plain text file:
#   <study_name>_methods_summary.txt   (<study_name>_methods_summary_chunkN.txt for one chunk)
# in the data folder, and the PDF report
#   <study_name>_report.pdf            (<study_name>_report_chunkN.pdf for one chunk)
# in the Derived_Variables folder (see R/derived_variables.R).

read_table_chr <- function(file, ...) {
  utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

# Readable name of a QC column
metric_label <- function(col) {
  labels <- c(
    passed_signal_quality   = "Signal intensity",
    passed_beadcount        = "Bead count",
    passed_detection_cutoff = "Detection p-value",
    passed_snp_outlier_bl4  = "SNP outliers",
    passed_sex_mismatch     = "Sex mismatch"
  )
  ifelse(col %in% names(labels), unname(labels[col]), paste0("EWAS tools: ", col))
}

# The rule behind a QC column, written with the thresholds of this run
metric_rule <- function(config, col) {
  q <- config$qc_thresholds %||% list()
  switch(
    col,
    passed_signal_quality = sprintf(
      "median log2 methylated or unmethylated intensity below %s",
      q$signal_intensity_threshold %||% 10.5),
    passed_beadcount = sprintf(
      "more than %s%% of probes with fewer than %s beads, or a mean bead count of %s or less",
      100 * (q$bead_fail_fraction %||% 0.05), q$bead_min_count %||% 5, q$mean_bead_threshold %||% 2),
    passed_detection_cutoff = sprintf(
      "more than %s%% of probes with a detection p-value above %s",
      100 * (q$detection_p_fail_fraction %||% 0.01), q$detection_p_threshold %||% 0.05),
    passed_snp_outlier_bl4 = sprintf(
      "SNP outlier score (ewastools) above %s", -(q$snp_outlier_threshold %||% 4)),
    passed_sex_mismatch =
      "sex predicted from the X and Y chromosome probes (ewastools) differs from the reported sex",
    "EWAS tools control metric (the threshold is set by the ewastools package)"
  )
}

# The numbers of the chunks that were run, as a list of tables
collect_report_numbers <- function(config, chunks) {
  targets <- read_chunk_assignments(config)
  
  # ---- samples ----
  read_decisions <- function(k) {
    d <- read_table_chr(
      qc_filtering_decisions_file(config, k),
      colClasses = c(sample_id = "character", failed_metrics_detected = "character",
                     failed_metrics_used_for_filtering = "character")
    )
    d$chunk <- k
    d
  }
  decisions <- do.call(rbind, lapply(chunks, read_decisions))
  decisions$failed_metrics_detected[is.na(decisions$failed_metrics_detected)] <- ""
  decisions$failed_metrics_used_for_filtering[is.na(decisions$failed_metrics_used_for_filtering)] <- ""
  pass_fail <- do.call(rbind, lapply(chunks, function(k) read_table_chr(qc_pass_fail_file(config, k))))
  removed   <- decisions$removed_from_pipeline
  
  per_chunk <- do.call(rbind, lapply(chunks, function(k) {
    r <- removed[decisions$chunk == k]
    data.frame(chunk = as.character(k), samples = length(r), removed = sum(r),
               kept = sum(!r), stringsAsFactors = FALSE)
  }))
  overview <- rbind(
    per_chunk,
    data.frame(chunk = "all", samples = sum(per_chunk$samples), removed = sum(per_chunk$removed),
               kept = sum(per_chunk$kept), stringsAsFactors = FALSE)
  )
  
  # one row per QC check: how many samples fail it, is it used, how many removed because of it
  metric_cols <- setdiff(colnames(pass_fail), c("sample_id", "Basename", "chunk"))
  used_cols   <- suppressWarnings(selected_filter_columns(config, pass_fail))
  why_used    <- strsplit(decisions$failed_metrics_used_for_filtering, ";", fixed = TRUE)
  metrics <- do.call(rbind, lapply(metric_cols, function(col) {
    data.frame(
      check = metric_label(col),
      column = col,
      rule = metric_rule(config, col),
      used = col %in% used_cols,
      failed = sum(qc_value_failed(pass_fail[[col]])),
      removed_because = sum(vapply(why_used[removed], function(w) col %in% w, logical(1))),
      stringsAsFactors = FALSE
    )
  }))
  
  # the EWAS tools control metrics as one group: samples failing at least one, and the detail
  ewas_cols <- metrics$column[grepl("^EWAS tools", metrics$check)]
  failed_ewas <- vapply(strsplit(decisions$failed_metrics_detected, ";", fixed = TRUE),
                        function(f) any(f %in% ewas_cols), logical(1))
  ew_rows <- metrics[metrics$column %in% ewas_cols, , drop = FALSE]
  ew_fail <- ew_rows[ew_rows$failed > 0, , drop = FALSE]
  ewas_summary <- list(
    n = length(ewas_cols), n_used = sum(ewas_cols %in% used_cols),
    failed = sum(failed_ewas), removed = sum(failed_ewas & removed),
    detail = if (nrow(ew_fail)) paste0(sub("^EWAS tools: ", "", ew_fail$check), " = ", ew_fail$failed, collapse = "; ")
    else "no sample fails any of them"
  )
  
  # the checks that matter most for sample identity and quality: how many samples failed each,
  # how many were removed, and how many failed but are still in the sample
  check_row <- function(label, pattern, used_pattern) {
    available <- any(grepl(pattern, metric_cols))
    failed_by <- grepl(pattern, decisions$failed_metrics_detected)
    data.frame(check = label, available = available, used = any(grepl(used_pattern, used_cols)),
               failed = sum(failed_by), removed = sum(failed_by & removed),
               still_in_sample = sum(failed_by & !removed), stringsAsFactors = FALSE)
  }
  checks <- rbind(
    check_row("Bisulfite conversion (EWAS tools, 3 metrics)", "^Bisulfite Conversion", "^Bisulfite Conversion"),
    check_row("SNP outliers (possible contamination or failed sample)", "^passed_snp_outlier", "^passed_snp_outlier"),
    check_row("Sex mismatch (predicted against reported sex)", "^passed_sex_mismatch", "^passed_sex_mismatch")
  )
  # result of the sex check for every sample (match, mismatch, inconclusive, reported sex missing)
  sex_check <- unlist(lapply(chunks, function(k) {
    t <- read_table_chr(qc_metrics_file(config, k))
    if ("sex_check" %in% colnames(t)) t$sex_check else NULL
  }))
  sex_results <- if (length(sex_check)) table(sex_check) else NULL
  
  # the removed samples, and the ones that failed a check that is not used for filtering
  name_col <- name_column(config)
  if (!name_col %in% colnames(targets)) {
    warning("The sample name column '", name_col, "' is not in the chunk table of step 01, so the names are empty ",
            "in the report. Rename the column in the sample sheet or set sample_sheet: name_col, then rerun from step 00.",
            call. = FALSE)
  }
  decisions$Sample_Name <- if (name_col %in% colnames(targets)) {
    targets[[name_col]][match(decisions$sample_id, targets$sample_id)]
  } else {
    rep(NA_character_, nrow(decisions))
  }
  sex_detail  <- sex_detail_from_metrics(config, chunks)
  removed_tbl <- decisions[removed, c("sample_id", "Sample_Name", "chunk", "failed_metrics_used_for_filtering"), drop = FALSE]
  colnames(removed_tbl)[4] <- "reasons"
  removed_tbl$reasons <- reason_with_sex(removed_tbl$reasons, removed_tbl$sample_id, sex_detail)
  flagged <- !removed & nzchar(decisions$failed_metrics_detected)
  flagged_tbl <- decisions[flagged, c("sample_id", "Sample_Name", "chunk", "failed_metrics_detected"), drop = FALSE]
  colnames(flagged_tbl)[4] <- "failed"
  flagged_tbl$failed <- reason_with_sex(flagged_tbl$failed, flagged_tbl$sample_id, sex_detail)
  
  # ---- probes ----
  attr_files <- vapply(chunks, function(k) {
    study_file(config, "probe_filter_attrition", chunk = k, ext = "csv")
  }, character(1))
  probes <- NULL            # the table of the first chunk
  probe_tables <- NULL      # the tables of all chunks, named by chunk
  probes_differ <- FALSE    # do the chunks have different numbers of probes?
  if (all(file.exists(attr_files))) {
    probe_tables <- lapply(attr_files, read_table_chr)
    names(probe_tables) <- as.character(chunks)
    probes <- probe_tables[[1]]
    probes_differ <- !all(vapply(probe_tables[-1], function(t) {
      identical(t$probes_remaining, probes$probes_remaining)
    }, logical(1)))
  }
  
  # ---- data for the sample QC plots: one row per sample, with its status ----
  qc <- tryCatch({
    d <- bind_chunk_tables(lapply(chunks, function(k) {
      read_table_chr(qc_metrics_file(config, k), colClasses = c(sample_id = "character"))
    }))
    status <- ifelse(decisions$removed_from_pipeline, "removed",
                     ifelse(nzchar(decisions$failed_metrics_detected), "flagged", "passed"))
    d$status <- factor(status[match(as.character(d$sample_id), decisions$sample_id)],
                       levels = names(RP_STATUS))
    d$status[is.na(d$status)] <- "passed"
    d
  }, error = function(e) NULL)
  # ---- FraminghamPACE probe coverage (saved by step 09, one file per chunk) ----
  fram <- tryCatch({
    fs <- vapply(chunks, function(k) study_file(config, "FraminghamPACE_coverage", chunk = k, ext = "csv"), character(1))
    fs <- fs[file.exists(fs)]
    if (length(fs)) do.call(rbind, lapply(fs, read_table_chr)) else NULL
  }, error = function(e) NULL)
  pca_file <- file.path(config$paths$output_data, "control_probe_pca_scores_all_chunks.csv")
  pca <- if (file.exists(pca_file)) tryCatch(read_table_chr(pca_file), error = function(e) NULL) else NULL
  
  # ---- probes removed in all chunks, each counted once ----
  # A probe that is removed in several chunks is counted once, at the first filter that removes it in
  # any chunk, so the steps add up: probes at the start - removed = probes kept in all chunks.
  probe_unique <- tryCatch({
    if (is.null(probe_tables)) stop("no probe tables")
    rf <- vapply(chunks, function(k) study_file(config, "probe_filter_removed", chunk = k, ext = "csv"), character(1))
    if (!all(file.exists(rf))) stop("no removed-probe files")
    rem <- lapply(rf, read_table_chr, colClasses = c(probe = "character"))
    filters <- probe_tables[[1]]$filter[-1]
    seen <- character()
    uniq <- numeric(length(filters))
    for (i in seq_along(filters)) {
      ids <- unique(unlist(lapply(rem, function(r) r$probe[r$filter == filters[i]])))
      uniq[i] <- length(setdiff(ids, seen))
      seen <- union(seen, ids)
    }
    start <- probe_tables[[1]]$probes_remaining[1]
    list(filters = filters, removed = uniq, start = start, total_removed = length(seen),
         shared_kept = start - length(seen))
  }, error = function(e) NULL)
  
  # ---- settings ----
  yn <- function(x) if (isTRUE(x)) "yes" else "no"
  st <- list()
  add_setting <- function(group, setting, value) {
    st[[length(st) + 1]] <<- c(group = group, setting = setting, value = paste(as.character(value), collapse = ", "))
  }
  if (!is.null(probes) && "array_type" %in% colnames(probes)) add_setting("Data", "Array", array_name(probes$array_type[1]))
  add_setting("Data", "Chunks", config$chunking$n_chunks %||% 1)
  add_setting("Data", "Chunk grouping columns", unlist(config$chunking$grouping_cols) %||% "none")
  add_setting("Data", "Random seed", config$chunking$seed %||% 123)
  pf <- config$probe_filtering %||% list()
  for (key in c("remove_snp_probes", "remove_cross_reactive", "remove_polymorphic", "remove_xy",
                "remove_failed_detection", "remove_failed_beadcount")) {
    add_setting("Probe filtering", key, yn(pf[[key]]))
  }
  for (key in c("bead_min_count", "bead_fail_fraction", "detection_p_threshold", "detection_p_fail_fraction")) {
    if (!is.null(pf[[key]])) add_setting("Probe filtering", key, pf[[key]])
  }
  add_setting("Normalization", "Noob", "yes")
  add_setting("Normalization", "BMIQ", "yes")
  nd <- config$noob_downstream %||% list()
  add_setting("Downstream modules", "Cell counts", yn(nd$cell_counts))
  if (isTRUE(nd$cell_counts)) add_setting("Downstream modules", "Cell reference", nd$cell_reference %||% "")
  add_setting("Downstream modules", "PC clocks", yn(nd$pc_clocks))
  add_setting("Downstream modules", "SystemsAge", yn(nd$systems_age))
  add_setting("Downstream modules", "DunedinPACE", yn(nd$dunedinpace))
  add_setting("Downstream modules", "FraminghamPACE (on BMIQ betas)", yn((config$bmiq_downstream %||% list())$framinghampace))
  settings <- as.data.frame(do.call(rbind, st), stringsAsFactors = FALSE)
  
  # ---- software ----
  pkgs <- c("minfi", "wateRmelon", "ewastools", "RPMM", "methylCIPHER", "DunedinPACE", "FraminghamPACE",
            "FlowSorted.Blood.EPIC", "FlowSorted.BloodExtended.EPIC",
            "IlluminaHumanMethylationEPICv2manifest", "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
            "IlluminaHumanMethylationEPICmanifest", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
            "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
            "ggplot2")
  versions <- vapply(pkgs, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  }, character(1))
  software <- data.frame(
    item = c("R", pkgs[!is.na(versions)]),
    version = c(paste(R.version$major, R.version$minor, sep = "."), unname(versions[!is.na(versions)])),
    stringsAsFactors = FALSE
  )
  
  list(overview = overview, checks = checks, sex_results = sex_results, ewas_summary = ewas_summary,
       metrics = metrics, removed = removed_tbl, flagged = flagged_tbl,
       probes = probes, probe_tables = probe_tables, probes_differ = probes_differ,
       settings = settings, software = software, qc = qc, pca = pca, name_col = name_col,
       id_col = config$sample_sheet$sample_id_col %||% "IDATid", framingham = fram,
       probe_unique = probe_unique,
       n_chunks_run = length(chunks), n_chunks_all = as.integer(config$chunking$n_chunks %||% 1))
}

# "EPICv2" instead of "epicv2" (falls back to the stored value)
array_name <- function(x) tryCatch(array_label(x), error = function(e) x)

# One row per chunk: samples, probes at the start, probes kept, probes removed
probe_chunk_summary <- function(nums) {
  tabs <- nums$probe_tables
  ov   <- nums$overview
  do.call(rbind, lapply(names(tabs), function(k) {
    t <- tabs[[k]]
    start <- t$probes_remaining[1]
    end   <- t$probes_remaining[nrow(t)]
    n <- ov$kept[ov$chunk == k]
    data.frame(chunk = k, samples = if (length(n)) n else NA, start = start, kept = end,
               removed = start - end, stringsAsFactors = FALSE)
  }))
}

# The filters as a table with the probes removed in each chunk, and a last row with the probes kept
probe_filters_wide <- function(nums, with_reason = FALSE) {
  tabs  <- nums$probe_tables
  first <- tabs[[1]]
  d     <- first[-1, , drop = FALSE]
  pu    <- nums$probe_unique
  n_row <- nrow(d) + 3                     # probes at the start, the filters, total removed, probes kept
  out   <- data.frame(Filter = c("Probes at the start", sub("^Remove ", "", d$filter), "Total removed", "Probes kept"),
                      stringsAsFactors = FALSE)
  for (k in names(tabs)) {
    t <- tabs[[k]]
    start <- t$probes_remaining[1]
    end   <- t$probes_remaining[nrow(t)]
    out[[paste0("Chunk ", k)]] <- c(rp_fmt(start), rp_fmt(t$probes_removed[match(d$filter, t$filter)]),
                                    rp_fmt(start - end), rp_fmt(end))
  }
  if (!is.null(pu)) {
    out[["All chunks (unique)"]] <- c(rp_fmt(pu$start), rp_fmt(pu$removed), rp_fmt(pu$total_removed), rp_fmt(pu$shared_kept))
  }
  if (with_reason) out$Reason <- c("", d$reason, "", "")
  out
}

# A data frame as aligned text lines
table_lines <- function(df) {
  utils::capture.output(print(df, row.names = FALSE, right = FALSE))
}

# The numbers as the lines of a text file
format_methods_summary <- function(config, nums, chunks) {
  lines <- character()
  add <- function(...) lines <<- c(lines, paste0(...))
  fmt <- function(x) format(x, big.mark = ",", trim = TRUE)
  
  add("METHODS SUMMARY: ", config$study_name %||% "")
  add("Written ", format(Sys.time(), "%Y-%m-%d %H:%M"), " by the GCC DNAm preprocessing pipeline.")
  add("All numbers are read from the output files of the pipeline steps.")
  if (nums$n_chunks_run < nums$n_chunks_all) {
    add("Only chunk(s) ", paste(chunks, collapse = ", "), " of ", nums$n_chunks_all,
        " were included: the numbers cover these chunks only.")
  }
  
  # ---- samples ----
  total <- nums$overview[nums$overview$chunk == "all", ]
  add("")
  add("1. SAMPLES")
  add("   Samples in the sample sheet: ", total$samples, " (", nums$n_chunks_run, " chunk(s))")
  add("   Removed by sample QC: ", total$removed)
  add("   Kept for analysis: ", total$kept)
  add("")
  add("   Samples per chunk:")
  lines <- c(lines, paste0("     ", table_lines(nums$overview)))
  add("")
  add("   Identity and control checks (samples that failed; removed; still in the sample):")
  ck <- nums$checks
  for (i in seq_len(nrow(ck))) {
    if (!ck$available[i]) {
      add("     ", ck$check[i], ": not available (this check was not run)")
    } else if (ck$used[i]) {
      add("     ", ck$check[i], ": ", ck$failed[i], " failed; ", ck$removed[i], " removed")
    } else {
      add("     ", ck$check[i], ": ", ck$failed[i], " failed; NOT used for filtering, so ",
          ck$still_in_sample[i], " still in the sample (please check them)")
    }
  }
  if (!is.null(nums$sex_results)) {
    add("     Sex check result of every sample: ",
        paste(names(nums$sex_results), as.integer(nums$sex_results), sep = " = ", collapse = ", "))
  }
  add("")
  add("   All sample QC checks, over all samples (a sample is removed if it fails a check that is used for filtering):")
  m   <- nums$metrics
  std <- m[!grepl("^EWAS tools", m$check), , drop = FALSE]
  for (i in seq_len(nrow(std))) {
    add("     ", std$check[i], ": ", std$failed[i], " failed; ",
        if (std$used[i]) paste0("used for filtering, ", std$removed_because[i], " removed because of it")
        else "NOT used for filtering")
    add("         rule: ", std$rule[i])
  }
  ew <- m[grepl("^EWAS tools", m$check), , drop = FALSE]
  if (nrow(ew)) {
    ew_fail <- ew[ew$failed > 0, , drop = FALSE]
    add("     EWAS tools control metrics (", nrow(ew), " metrics, ", sum(ew$used), " used for filtering)")
    add("         rule: ", ew$rule[1])
    add("         samples failing each: ",
        if (nrow(ew_fail)) paste0(sub("^EWAS tools: ", "", ew_fail$check), " = ", ew_fail$failed, collapse = "; ")
        else "no sample fails any of them")
  }
  add("")
  rt <- nums$removed
  add("   Removed samples (", nrow(rt), "):")
  if (nrow(rt)) {
    for (i in seq_len(nrow(rt))) {
      add("     ", rt$sample_id[i],
          if (!is.na(rt$Sample_Name[i])) paste0(" (", nums$name_col, " ", rt$Sample_Name[i], ")"),
          " [chunk ", rt$chunk[i], "]: ", rt$reasons[i])
    }
  } else {
    add("     none")
  }
  add("")
  ft <- nums$flagged
  add("   Flagged but kept (failed a check that is not used for filtering) (", nrow(ft), "):")
  if (nrow(ft)) {
    for (i in seq_len(nrow(ft))) {
      add("     ", ft$sample_id[i],
          if (!is.na(ft$Sample_Name[i])) paste0(" (", nums$name_col, " ", ft$Sample_Name[i], ")"),
          " [chunk ", ft$chunk[i], "]: ", ft$failed[i])
    }
  } else {
    add("     none")
  }
  
  # ---- probes ----
  add("")
  add("2. PROBES")
  pr <- nums$probes
  if (is.null(pr)) {
    add("   The probe filtering table was not found (run step 07).")
  } else if (isTRUE(nums$probes_differ)) {
    # the probes were filtered chunk by chunk and the chunks ended up with different probes
    tabs <- nums$probe_tables
    if ("array_type" %in% colnames(pr)) add("   Array: ", array_name(pr$array_type[1]))
    add("   The probes were filtered separately in each chunk. The bead count and detection p criteria are")
    add("   judged on the samples of the chunk, so the chunks can end up with different probes.")
    add("")
    add("   Probes per chunk:")
    cs <- probe_chunk_summary(nums)
    for (i in seq_len(nrow(cs))) {
      add("     Chunk ", cs$chunk[i], " (", cs$samples[i], " samples): ", fmt(cs$start[i]), " probes in the data, ",
          fmt(cs$kept[i]), " kept (", sprintf("%.2f%%", 100 * cs$kept[i] / cs$start[i]), "), ",
          fmt(cs$removed[i]), " removed")
    }
    add("")
    add("   Filters, in the order they were applied (probes removed in each chunk):")
    for (i in seq_len(nrow(pr))[-1]) {
      vals <- vapply(tabs, function(t) fmt(t$probes_removed[match(pr$filter[i], t$filter)]), character(1))
      add(sprintf("     %-46s removed per chunk: %s", pr$filter[i], paste(vals, collapse = " | ")))
      add("         reason: ", pr$reason[i])
    }
    pu <- nums$probe_unique
    if (!is.null(pu)) {
      add("")
      add("   All chunks together (a probe that is removed in several chunks is counted once, at the first filter that removes it):")
      add("     Probes at the start: ", fmt(pu$start))
      for (i in seq_along(pu$filters)) {
        add(sprintf("     %-46s removed %9s", pu$filters[i], fmt(pu$removed[i])))
      }
      add("     Total unique probes removed: ", fmt(pu$total_removed), " (",
          sprintf("%.2f%%", 100 * pu$total_removed / pu$start), " of the starting probes)")
      add("     Probes kept in all chunks: ", fmt(pu$shared_kept))
    }
  } else {
    start <- pr$probes_remaining[1]
    end   <- pr$probes_remaining[nrow(pr)]
    if ("array_type" %in% colnames(pr)) add("   Array: ", array_name(pr$array_type[1]))
    add("   Probes in the data: ", fmt(start))
    add("   Probes kept after all filters: ", fmt(end), " (", sprintf("%.2f%%", 100 * end / start),
        " of the starting probes); removed: ", fmt(start - end))
    add("")
    add("   Filters, in the order they were applied (% = share of the starting probes):")
    for (i in seq_len(nrow(pr))[-1]) {
      add(sprintf("     %-46s removed %9s (%5.2f%%) | remaining %9s",
                  pr$filter[i], fmt(pr$probes_removed[i]), pr$percent_of_start[i], fmt(pr$probes_remaining[i])))
      add("         reason: ", pr$reason[i])
    }
  }
  
  # ---- settings and software ----
  add("")
  add("3. SETTINGS")
  st <- nums$settings
  for (g in unique(st$group)) {
    add("   ", g, ":")
    s <- st[st$group == g, , drop = FALSE]
    for (i in seq_len(nrow(s))) add("     ", s$setting[i], ": ", s$value[i])
  }
  add("")
  add("4. SOFTWARE")
  sw <- nums$software
  for (i in seq_len(nrow(sw))) add("   ", sw$item[i], " ", sw$version[i])
  
  lines
}

write_methods_summary <- function(config, nums, chunks, suffix = "") {
  lines <- format_methods_summary(config, nums, chunks)
  file  <- file.path(config$paths$output_data, paste0(study_prefix(config), "methods_summary", suffix, ".txt"))
  writeLines(lines, file)
  message("Methods summary written: ", file)
  invisible(file)
}

# ---------- REPORT: THE PDF ------------------------------------------------------------
# A4 pages drawn with the grid package (part of base R) and ggplot2, so no extra software
# is needed. Page 1: overview. Page 2: samples. Page 3: probes.

PIPELINE_URL <- "https://github.com/CPRyan/GCC-DNAm-preprocessing-pipeline"

RP <- list(navy = "#16324F", blue = "#1F4E79", teal = "#1B998B", light = "#EEF3F8",
           dark = "#1F2933", grey = "#6B7280", red = "#C8553D", orange = "#E09F3E")
RP_STATUS <- c(passed = "#2F6B9A", flagged = "#E09F3E", removed = "#C8553D")   # colours of the sample status
RP_W  <- 8.27    # page width in inches (A4 portrait)
RP_H  <- 11.69   # page height in inches
RP_X0 <- 0.07    # left and right margins, as a fraction of the page width
RP_X1 <- 0.93

rp_fmt <- function(x) format(x, big.mark = ",", trim = TRUE)

# width in inches of each string at this font size
rp_text_width <- function(s, fontsize, face = "plain") {
  vapply(s, function(x) {
    grid::convertWidth(
      grid::grobWidth(grid::textGrob(x, gp = grid::gpar(fontsize = fontsize, fontface = face))),
      "inches", valueOnly = TRUE
    )
  }, numeric(1), USE.NAMES = FALSE)
}

# New page: a dark header band with the title (and a link and a subtitle below it), a thin
# accent stripe, and a footer with the study and the page number.
# Returns the y position where content starts.
rp_page <- function(title, subtitle, page, n_pages, study, link = NULL) {
  grid::grid.newpage()
  band <- if (is.null(link)) 0.075 else if (is.null(subtitle)) 0.080 else 0.098
  grid::grid.rect(x = 0, y = 1, width = 1, height = band, just = c("left", "top"),
                  gp = grid::gpar(fill = RP$navy, col = NA))
  grid::grid.rect(x = 0, y = 1 - band, width = 1, height = 0.004, just = c("left", "top"),
                  gp = grid::gpar(fill = RP$teal, col = NA))
  grid::grid.text(title, x = RP_X0, y = 0.986, just = c("left", "top"),
                  gp = grid::gpar(fontsize = 17, fontface = "bold", col = "white"))
  y_next <- 0.959
  if (!is.null(link)) {
    grid::grid.text(link, x = RP_X0, y = y_next, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 8.5, col = "#9CC7E0"))
    y_next <- y_next - 0.019
  }
  if (!is.null(subtitle)) {
    grid::grid.text(subtitle, x = RP_X0, y = y_next, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 9, col = "#C9D6E2"))
  }
  grid::grid.lines(x = c(RP_X0, RP_X1), y = 0.040, gp = grid::gpar(col = "grey80", lwd = 0.8))
  grid::grid.text(paste0(study, "   |   GCC DNAm preprocessing pipeline"), x = RP_X0, y = 0.025,
                  just = c("left", "center"), gp = grid::gpar(fontsize = 7, col = RP$grey))
  grid::grid.text(paste0("page ", page, " of ", n_pages), x = RP_X1, y = 0.025,
                  just = c("right", "center"), gp = grid::gpar(fontsize = 7, col = RP$grey))
  1 - band - 0.026
}

rp_heading <- function(text, y, x = RP_X0) {
  grid::grid.text(text, x = x, y = y, just = c("left", "top"),
                  gp = grid::gpar(fontsize = 11, fontface = "bold", col = RP$blue))
  y - 0.026
}

# Lines of plain text, one below the other. Returns the new y.
rp_lines <- function(lines, y, fontsize = 8.5, col = RP$dark) {
  lh <- fontsize * 1.5 / 72 / RP_H
  for (l in lines) {
    grid::grid.text(l, x = RP_X0, y = y, just = c("left", "top"),
                    gp = grid::gpar(fontsize = fontsize, col = col))
    y <- y - lh
  }
  y - lh * 0.4
}

# A row of boxes with a big number. cards: list of list(label, value, note). Returns the new y.
rp_cards <- function(cards, y, height = 0.085) {
  n   <- length(cards)
  gap <- 0.02
  w   <- ((RP_X1 - RP_X0) - gap * (n - 1)) / n
  inner <- w * RP_W - 0.30                 # inches for the big number
  for (i in seq_len(n)) {
    x <- RP_X0 + (i - 1) * (w + gap)
    value <- as.character(cards[[i]]$value)
    size  <- max(9, min(17, 17 * inner / rp_text_width(value, 17, "bold")))
    grid::grid.rect(x = x, y = y, width = w, height = height, just = c("left", "top"),
                    gp = grid::gpar(fill = RP$light, col = NA))
    grid::grid.rect(x = x, y = y, width = 0.006, height = height, just = c("left", "top"),
                    gp = grid::gpar(fill = RP$blue, col = NA))
    grid::grid.text(cards[[i]]$label, x = x + 0.02, y = y - 0.012, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 8, col = RP$grey))
    grid::grid.text(value, x = x + 0.02, y = y - 0.030, just = c("left", "top"),
                    gp = grid::gpar(fontsize = size, fontface = "bold", col = RP$dark))
    if (!is.null(cards[[i]]$note)) {
      grid::grid.text(cards[[i]]$note, x = x + 0.02, y = y - 0.066, just = c("left", "top"),
                      gp = grid::gpar(fontsize = 7.5, col = RP$grey))
    }
  }
  y - height - 0.025
}

# A table with a coloured header and light rows. Columns named in `wrap` are wrapped at that
# many characters (they are narrowed further if the table is wider than the page). At most
# max_rows rows are drawn. Returns the new y.
rp_table <- function(df, y, fontsize = 8, align = NULL, wrap = NULL, max_rows = 30, empty = "none",
                     more_note = "see the tables in the data folder", x0 = RP_X0, x1 = RP_X1) {
  if (!nrow(df)) return(rp_lines(empty, y, fontsize + 0.5))
  more <- 0L
  if (nrow(df) > max_rows) {
    more <- nrow(df) - max_rows
    df <- df[seq_len(max_rows), , drop = FALSE]
  }
  df[] <- lapply(df, as.character)
  df[is.na(df)] <- ""
  nc    <- ncol(df)
  align <- align %||% rep("l", nc)
  avail <- (x1 - x0) * RP_W
  pad   <- 0.07                            # inches left and right of every cell
  lh    <- fontsize * 1.3 / 72 / RP_H      # one text line, as a fraction of the page height
  vpad  <- fontsize * 0.45 / 72 / RP_H     # space above and below the text of a row
  
  wrap_n <- rep(Inf, nc)
  if (!is.null(wrap)) wrap_n[match(names(wrap), colnames(df))] <- wrap
  make_cells <- function(wn) {
    lapply(seq_len(nc), function(j) {
      lapply(df[[j]], function(s) {
        if (!is.finite(wn[j])) return(s)
        w <- strwrap(s, width = wn[j])
        if (length(w)) w else ""
      })
    })
  }
  repeat {
    cells  <- make_cells(wrap_n)
    widths <- vapply(seq_len(nc), function(j) {
      max(c(rp_text_width(colnames(df)[j], fontsize, "bold"),
            rp_text_width(unlist(cells[[j]]), fontsize))) + 2 * pad
    }, numeric(1))
    if (sum(widths) <= avail || !any(is.finite(wrap_n))) break
    j <- which.max(ifelse(is.finite(wrap_n), widths, -Inf))
    if (wrap_n[j] <= 14) break
    wrap_n[j] <- floor(wrap_n[j] * 0.85)
  }
  if (sum(widths) > avail) widths <- widths * avail / sum(widths)
  x_edge <- x0 + c(0, cumsum(widths)) / RP_W
  table_w <- sum(widths) / RP_W
  text_x  <- function(j) if (align[j] == "r") x_edge[j + 1] - pad / RP_W else x_edge[j] + pad / RP_W
  text_h  <- function(j) if (align[j] == "r") "right" else "left"
  
  hh <- lh + 2 * vpad
  grid::grid.rect(x = x0, y = y, width = table_w, height = hh, just = c("left", "top"),
                  gp = grid::gpar(fill = RP$blue, col = NA))
  for (j in seq_len(nc)) {
    grid::grid.text(colnames(df)[j], x = text_x(j), y = y - vpad, just = c(text_h(j), "top"),
                    gp = grid::gpar(fontsize = fontsize, fontface = "bold", col = "white"))
  }
  y <- y - hh
  
  for (i in seq_len(nrow(df))) {
    n_lines <- max(vapply(seq_len(nc), function(j) length(cells[[j]][[i]]), integer(1)))
    rh <- n_lines * lh + 2 * vpad
    if (i %% 2 == 0) {
      grid::grid.rect(x = x0, y = y, width = table_w, height = rh, just = c("left", "top"),
                      gp = grid::gpar(fill = RP$light, col = NA))
    }
    for (j in seq_len(nc)) {
      for (k in seq_along(cells[[j]][[i]])) {
        grid::grid.text(cells[[j]][[i]][k], x = text_x(j), y = y - vpad - (k - 1) * lh,
                        just = c(text_h(j), "top"), gp = grid::gpar(fontsize = fontsize, col = RP$dark))
      }
    }
    y <- y - rh
  }
  grid::grid.lines(x = c(x0, x0 + table_w), y = y, gp = grid::gpar(col = "grey70", lwd = 0.8))
  y <- y - 0.012
  if (more > 0) {
    y <- rp_lines(paste0("... and ", more, " more (", more_note, ")"), y, fontsize, RP$grey)
  }
  y - 0.012
}

# The three identity and control checks as a table
rp_checks_df <- function(nums) {
  ck <- nums$checks
  data.frame(
    Check = ck$check,
    Failed = ifelse(ck$available, as.character(ck$failed), "-"),
    Removed = ifelse(ck$available & ck$used, as.character(ck$removed), "-"),
    `Still in the sample` = ifelse(ck$available & !ck$used, as.character(ck$still_in_sample), "-"),
    Status = ifelse(!ck$available, "not run", ifelse(ck$used, "removed", "flag only")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# The probe filters as a table (the starting row is left out)
rp_probes_df <- function(pr, with_reason = FALSE) {
  d <- pr[-1, , drop = FALSE]
  out <- data.frame(
    Filter = sub("^Remove ", "", d$filter),
    Removed = rp_fmt(d$probes_removed),
    `% of start` = sprintf("%.2f%%", d$percent_of_start),
    Remaining = rp_fmt(d$probes_remaining),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (with_reason) out$Reason <- d$reason
  out
}

# ---- the time of the run ----
# run_all_steps.R puts the start of the run and its steps in the environment. Without it (the
# steps were started one by one), only the time of the report is known.
rp_run_info <- function() {
  now <- Sys.time()
  start_txt <- Sys.getenv("GCC_RUN_START", unset = "")
  start <- if (nzchar(start_txt)) tryCatch(as.POSIXct(start_txt, format = "%Y-%m-%d %H:%M:%S"), error = function(e) NA) else NA
  if (length(start) != 1 || is.na(start)) return(list(known = FALSE, end = now))
  list(known = TRUE, start = start, end = now, steps = Sys.getenv("GCC_RUN_STEPS", unset = ""),
       minutes = as.numeric(difftime(now, start, units = "mins")))
}

rp_format_duration <- function(minutes) {
  if (minutes < 1) return("< 1 min")
  h <- floor(minutes / 60)
  m <- round(minutes - 60 * h)
  if (m == 60) { h <- h + 1; m <- 0 }
  if (h > 0) sprintf("%d h %02d min", h, m) else sprintf("%d min", m)
}

# ---- page 1: overview ----
rp_page_overview <- function(config, nums, page, n_pages, toc = NULL) {
  study <- config$study_name %||% ""
  pr    <- nums$probes
  array_txt <- if (!is.null(pr) && "array_type" %in% colnames(pr)) array_name(pr$array_type[1]) else "unknown"
  run <- rp_run_info()
  y <- rp_page("GCC DNAm preprocessing pipeline - Run Overview", NULL, page, n_pages, study, link = PIPELINE_URL)
  
  # ---- row 1: the run ----
  row1 <- list(
    list(label = "STUDY", value = study, note = paste0("report of ", format(run$end, "%Y-%m-%d"))),
    list(label = "ARRAY", value = array_txt, note = "detected from the data")
  )
  if (run$known) {
    if (nzchar(run$steps)) {
      row1[[length(row1) + 1]] <- list(label = "STEPS RUN", value = run$steps, note = "of the pipeline")
    }
    same_day <- identical(as.Date(run$start), as.Date(run$end))
    row1[[length(row1) + 1]] <- list(
      label = "RUN TIME", value = rp_format_duration(run$minutes),
      note = paste0(format(run$start, "%Y-%m-%d %H:%M"), " - ", format(run$end, if (same_day) "%H:%M" else "%Y-%m-%d %H:%M")))
  } else {
    row1[[length(row1) + 1]] <- list(label = "REPORT CREATED", value = format(run$end, "%Y-%m-%d"),
                                     note = format(run$end, "%H:%M:%S"))
  }
  y <- rp_cards(row1, y)
  
  # ---- row 2: the results ----
  total <- nums$overview[nums$overview$chunk == "all", ]
  probes_card <- NULL
  if (!is.null(nums$probe_tables)) {
    cs <- probe_chunk_summary(nums)
    pu <- nums$probe_unique
    probes_card <- if (!is.null(pu) && nrow(cs) > 1) {
      # the probes that all chunks share: the "all" row of the probes table on the probes page
      list(label = "PROBES KEPT", value = rp_fmt(pu$shared_kept), note = paste0("in all chunks, of ", rp_fmt(pu$start)))
    } else if (nrow(cs) > 1 && length(unique(cs$kept)) > 1) {
      list(label = "PROBES KEPT", value = rp_fmt(min(cs$kept)), note = paste0("fewest of ", nrow(cs), " chunks"))
    } else {
      list(label = "PROBES KEPT", value = rp_fmt(cs$kept[1]), note = paste0("of ", rp_fmt(cs$start[1])))
    }
  }
  row2 <- c(
    list(list(label = "SAMPLES KEPT", value = paste0(total$kept, " of ", total$samples),
              note = paste0(total$removed, " removed by sample QC")),
         list(label = "FLAGGED, KEPT", value = as.character(nrow(nums$flagged)), note = "failed a check, not removed")),
    if (!is.null(probes_card)) list(probes_card),
    list(list(label = "CHUNKS", value = as.character(nums$n_chunks_run), note = paste0("of ", nums$n_chunks_all)))
  )
  y <- rp_cards(row2, y)
  
  # ---- what the next pages show ----
  if (!is.null(toc) && nrow(toc)) {
    y <- rp_heading("Contents", y)
    y <- rp_table(toc, y, fontsize = 8.5, align = c("r", "l", "l"), wrap = c(`What it shows` = 75))
  }
  rp_lines("All numbers are read from the output files of the pipeline steps.", y, fontsize = 7.5, col = RP$grey)
  invisible(NULL)
}

# ---- page 2: samples ----
rp_page_samples <- function(config, nums, page, n_pages) {
  study <- config$study_name %||% ""
  y <- rp_page("Samples", "Which samples were removed, and why", page, n_pages, study)
  
  y <- rp_heading("Samples per chunk", y)
  ov <- nums$overview
  colnames(ov) <- c("Chunk", "Samples", "Removed", "Kept")
  y <- rp_table(ov, y, align = c("l", "r", "r", "r"), x1 = RP_X0 + 0.45)
  
  y <- rp_heading("Identity and control checks", y)
  y <- rp_table(rp_checks_df(nums), y, align = c("l", "r", "r", "r", "l"))
  
  y <- rp_heading("All sample QC checks", y)
  m   <- nums$metrics
  std <- m[!grepl("^EWAS tools", m$check), , drop = FALSE]
  df <- data.frame(Check = std$check, Rule = std$rule, Used = ifelse(std$used, "yes", "no"),
                   Failed = std$failed, Removed = std$removed_because, stringsAsFactors = FALSE)
  ew <- nums$ewas_summary
  if (ew$n > 0) {
    df <- rbind(df, data.frame(
      Check = paste0("EWAS tools control metrics (", ew$n, ")"),
      Rule = "Illumina control probe metrics; thresholds are set by the ewastools package",
      Used = paste0(ew$n_used, " of ", ew$n), Failed = ew$failed, Removed = ew$removed,
      stringsAsFactors = FALSE))
  }
  y <- rp_table(df, y, align = c("l", "l", "l", "r", "r"), wrap = c(Rule = 62))
  y <- rp_lines(c("Failed = samples failing the check. Removed = samples removed because of it. For the EWAS tools metrics,",
                  "these are samples failing at least one of the metrics.",
                  paste0("EWAS tools metrics with failing samples: ", ew$detail)),
                y, fontsize = 7.5, col = RP$grey)
  
  rt <- nums$removed
  y <- rp_heading(paste0("Removed samples (", nrow(rt), ")"), y)
  if (nrow(rt)) {
    rd <- data.frame(`Sample ID` = rt$sample_id, Name = rt$Sample_Name, Chunk = rt$chunk,
                     Reason = rt$reasons, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(rd)[1:2] <- c(nums$id_col, nums$name_col)
    y <- rp_table(rd, y, align = c("l", "l", "r", "l"), wrap = c(Reason = 55), max_rows = 25)
  } else {
    y <- rp_lines("No sample was removed.", y)
  }
  
  ft <- nums$flagged
  y <- rp_heading(paste0("Flagged but kept (", nrow(ft), ")"), y)
  if (nrow(ft)) {
    fd <- data.frame(`Sample ID` = ft$sample_id, Name = ft$Sample_Name, Chunk = ft$chunk,
                     `Failed check` = ft$failed, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(fd)[1:2] <- c(nums$id_col, nums$name_col)
    y <- rp_table(fd, y, align = c("l", "l", "r", "l"), wrap = c(`Failed check` = 55), max_rows = 15)
    rp_lines("These samples failed a check that is not used for filtering. Please check them.",
             y, fontsize = 7.5, col = RP$grey)
  } else {
    rp_lines("No sample failed a check without being removed.", y)
  }
  invisible(NULL)
}

# ---- page 3: probes ----
rp_page_probes <- function(config, nums, page, n_pages) {
  study <- config$study_name %||% ""
  y <- rp_page("Probes", "Which probes were removed, and why", page, n_pages, study)
  pr <- nums$probes
  if (is.null(pr)) {
    rp_lines("The probe filtering table was not found (run step 07).", y)
    return(invisible(NULL))
  }
  tabs   <- nums$probe_tables
  differ <- isTRUE(nums$probes_differ)
  
  if (differ) {
    cs <- probe_chunk_summary(nums)
    y <- rp_lines(c("The probes were filtered separately in each chunk. The bead count and detection p criteria",
                    "are judged on the samples of the chunk, so the chunks can end up with different probes."),
                  y, fontsize = 9)
    if (length(tabs) > 6) {
      y <- rp_heading("Probes per chunk", y)
      cd <- data.frame(Chunk = cs$chunk, Samples = cs$samples, `Probes at start` = rp_fmt(cs$start),
                       `Probes kept` = rp_fmt(cs$kept), Removed = rp_fmt(cs$removed),
                       check.names = FALSE, stringsAsFactors = FALSE)
      rp_table(cd, y, align = c("l", "r", "r", "r", "r"), max_rows = 40)
      return(invisible(NULL))
    }
    y <- rp_heading("Probes removed by each filter, in each chunk", y)
    y <- rp_table(probe_filters_wide(nums, with_reason = TRUE), y,
                  align = c("l", rep("r", length(tabs) + !is.null(nums$probe_unique)), "l"), wrap = c(Reason = 45))
    if (!is.null(nums$probe_unique)) {
      y <- rp_lines(strwrap(paste("All chunks (unique): a probe that is removed in several chunks is counted once,",
                                  "at the first filter that removes it, so the steps add up. The probes kept in all chunks",
                                  "are the probes that the chunks share."), width = 140),
                    y, fontsize = 7, col = RP$grey)
    }
    long <- do.call(rbind, lapply(names(tabs), function(k) {
      d <- tabs[[k]][-1, , drop = FALSE]
      data.frame(filter = sub("^Remove ", "", d$filter), chunk = k, pct = d$percent_of_start,
                 stringsAsFactors = FALSE)
    }))
    long$filter <- factor(long$filter, levels = rev(unique(long$filter)))
    long$chunk  <- factor(long$chunk, levels = names(tabs))
    fills <- rep_len(c("#5C7FA3", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7"), nlevels(long$chunk))
    chart <- ggplot2::ggplot(long, ggplot2::aes(x = pct, y = filter, fill = chunk)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
      ggplot2::scale_fill_manual(values = fills, name = "Chunk") +
      ggplot2::labs(x = "Probes removed (% of the starting probes)", y = NULL,
                    title = "Probes removed by each filter") +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(face = "bold"))
  } else {
    start <- pr$probes_remaining[1]
    end   <- pr$probes_remaining[nrow(pr)]
    y <- rp_lines(c(paste0("Probes in the data: ", rp_fmt(start)),
                    paste0("Probes kept after all filters: ", rp_fmt(end), " (",
                           sprintf("%.2f%%", 100 * end / start), " of the starting probes); removed: ",
                           rp_fmt(start - end))),
                  y, fontsize = 9)
    
    y <- rp_heading("Filters, in the order they were applied", y)
    y <- rp_table(rp_probes_df(pr, with_reason = TRUE), y, align = c("l", "r", "r", "r", "l"),
                  wrap = c(Reason = 50))
    
    d <- pr[-1, , drop = FALSE]
    d$filter <- factor(sub("^Remove ", "", d$filter), levels = rev(sub("^Remove ", "", d$filter)))
    d$label  <- paste0(rp_fmt(d$probes_removed), "  (", sprintf("%.2f%%", d$percent_of_start), ")")
    chart <- ggplot2::ggplot(d, ggplot2::aes(x = percent_of_start, y = filter)) +
      ggplot2::geom_col(fill = RP$blue, width = 0.65) +
      ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.08, size = 3) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.4))) +
      ggplot2::labs(x = "Probes removed (% of the starting probes)", y = NULL,
                    title = "Probes removed by each filter") +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(face = "bold"))
  }
  chart_h <- 0.30
  print(chart, vp = grid::viewport(x = 0.5, y = max(y - 0.02, chart_h + 0.05), width = 0.86, height = chart_h,
                                   just = c("center", "top")))
  invisible(NULL)
}

# ---- plots ----

# A saved plot (list(plot, width, height)) at the full content width, as high as its shape needs
# but not higher than max_height (a fraction of the page). Returns the new y.
rp_place_plot <- function(fig, y, max_height = 0.8) {
  w_npc <- RP_X1 - RP_X0
  h_npc <- (fig$height / fig$width) * (w_npc * RP_W) / RP_H
  if (h_npc > max_height) {
    w_npc <- w_npc * max_height / h_npc
    h_npc <- max_height
  }
  print(fig$plot, vp = grid::viewport(x = 0.5, y = y, width = w_npc, height = h_npc, just = c("center", "top")))
  y - h_npc - 0.015
}

# A page with one to three saved plots below each other
rp_page_figures <- function(config, figs, keys, title, subtitle, page, n_pages, reserve = 0, extra = NULL) {
  y <- rp_page(title, subtitle, page, n_pages, config$study_name %||% "")
  keys <- keys[keys %in% names(figs)]
  room <- y - 0.06 - reserve
  for (k in keys) {
    y <- rp_place_plot(figs[[k]], y, max_height = (room - 0.015 * (length(keys) - 1)) / length(keys))
  }
  if (!is.null(extra)) extra(y)
  invisible(NULL)
}

# The colours of passed, flagged and removed samples, as a legend line. Returns the new y.
rp_status_legend <- function(y, keep = names(RP_STATUS)) {
  labs <- c(passed = "passed all checks", flagged = "flagged, kept in the analysis", removed = "removed")[keep]
  x <- RP_X0
  for (nm in names(labs)) {
    grid::grid.points(x = x + 0.006, y = y - 0.006, pch = 16, size = grid::unit(5, "pt"),
                      gp = grid::gpar(col = RP_STATUS[[nm]]))
    grid::grid.text(labs[[nm]], x = x + 0.017, y = y - 0.006, just = c("left", "center"),
                    gp = grid::gpar(fontsize = 8, col = RP$dark))
    x <- x + 0.017 + 0.0072 * nchar(labs[[nm]]) + 0.035
  }
  y - 0.025
}

rp_qc_theme <- function() {
  ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", size = 8.5),
                   axis.title = ggplot2::element_text(size = 7.5),
                   legend.position = "none",
                   plot.margin = ggplot2::margin(4, 8, 4, 6))
}

# One value per sample, sorted, coloured by status, with an optional threshold line
rp_qc_dot <- function(d, col, title, ylab, thr = NULL, subtitle = NULL, plain_labels = FALSE) {
  if (is.null(d) || !col %in% colnames(d)) return(NULL)
  d$val <- suppressWarnings(as.numeric(d[[col]]))
  d <- d[!is.na(d$val), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d <- d[order(d$val), , drop = FALSE]
  d$idx <- seq_len(nrow(d))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = idx, y = val, colour = status)) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_colour_manual(values = RP_STATUS, drop = FALSE) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Samples (sorted)", y = ylab) +
    rp_qc_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 6.5, colour = RP$grey))
  if (plain_labels) {
    p <- p + ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = FALSE, trim = TRUE))
  }
  if (!is.null(thr)) {
    p <- p + ggplot2::geom_hline(yintercept = thr, linetype = "dashed", colour = RP$grey, linewidth = 0.4)
  }
  p
}

rp_qc_signal <- function(d, thr) {
  if (is.null(d) || !all(c("MQC", "UQC") %in% colnames(d))) return(NULL)
  d$MQC <- suppressWarnings(as.numeric(d$MQC))
  d$UQC <- suppressWarnings(as.numeric(d$UQC))
  ggplot2::ggplot(d, ggplot2::aes(x = UQC, y = MQC, colour = status)) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_manual(values = RP_STATUS, drop = FALSE) +
    ggplot2::geom_hline(yintercept = thr, linetype = "dashed", colour = RP$grey, linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = thr, linetype = "dashed", colour = RP$grey, linewidth = 0.4) +
    ggplot2::labs(title = "Signal intensity", x = "Unmethylated signal (median log2)",
                  y = "Methylated signal (median log2)") +
    rp_qc_theme()
}

# X against Y signal, coloured by the reported sex; samples that are not "passed" get a ring
rp_qc_sex <- function(d) {
  if (is.null(d) || !all(c("sex_x_intensity", "sex_y_intensity") %in% colnames(d))) return(NULL)
  d$xval <- suppressWarnings(as.numeric(d$sex_x_intensity))
  d$yval <- suppressWarnings(as.numeric(d$sex_y_intensity))
  d <- d[!is.na(d$xval) & !is.na(d$yval), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  rep_sex <- if ("reported_sex" %in% colnames(d)) d$reported_sex else rep(NA_character_, nrow(d))
  d$sex <- ifelse(rep_sex %in% "f", "reported female", ifelse(rep_sex %in% "m", "reported male", "unknown"))
  odd <- d[d$status != "passed", , drop = FALSE]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = xval, y = yval, colour = sex)) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_manual(
      values = c("reported female" = "#CC79A7", "reported male" = "#0072B2", unknown = "#9CA3AF"), name = NULL) +
    ggplot2::labs(title = "Sex check", x = "X chromosome signal (normalised)",
                  y = "Y chromosome signal (normalised)") +
    rp_qc_theme() +
    ggplot2::theme(legend.position = "bottom", legend.key.size = grid::unit(8, "pt"),
                   legend.text = ggplot2::element_text(size = 7))
  if (nrow(odd)) {
    p <- p + ggplot2::geom_point(data = odd, ggplot2::aes(x = xval, y = yval), inherit.aes = FALSE,
                                 shape = 1, size = 4, stroke = 0.8, colour = RP$red)
  }
  p
}

# The first two principal components of the control probes
rp_qc_pca <- function(pcs, d) {
  if (is.null(pcs) || is.null(d)) return(NULL)
  id <- intersect(c("sample_id", "Sample", "sample", "Basename"), colnames(pcs))
  pc <- grep("^PC[0-9]+$", colnames(pcs), value = TRUE)
  if (!length(id) || length(pc) < 2) return(NULL)
  st <- as.character(d$status[match(as.character(pcs[[id[1]]]), as.character(d$sample_id))])
  pcs$status <- factor(ifelse(is.na(st), "passed", st), levels = names(RP_STATUS))
  pcs <- pcs[pcs$status != "removed", , drop = FALSE]   # removed samples are not shown
  pcs$xval <- suppressWarnings(as.numeric(pcs[[pc[1]]]))
  pcs$yval <- suppressWarnings(as.numeric(pcs[[pc[2]]]))
  ggplot2::ggplot(pcs, ggplot2::aes(x = xval, y = yval, colour = status)) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_manual(values = RP_STATUS, drop = FALSE) +
    ggplot2::labs(title = "Control-probe PCA", x = pc[1], y = pc[2]) +
    rp_qc_theme()
}

# ---- FraminghamPACE: probes and predictive coverage ----
# One row per model. If the chunks have the same numbers, they are one row for all chunks.
rp_framingham_df <- function(fc) {
  if (is.null(fc) || !nrow(fc)) return(NULL)
  cols <- c("model", "chunk", "background_probes_present", "background_probes_total", "scoring_probes_present",
            "scoring_probes_total", "predictive_coverage_min", "predictive_coverage_max", "samples_scored", "samples_total")
  if (!all(cols %in% colnames(fc))) return(NULL)
  fc$key <- paste(fc$model, fc$background_probes_present, fc$background_probes_total, fc$scoring_probes_present,
                  fc$scoring_probes_total, fc$predictive_coverage_min, fc$predictive_coverage_max)
  rows <- list()
  for (m in unique(fc$model)) {
    d <- fc[fc$model == m, , drop = FALSE]
    groups <- if (length(unique(d$key)) == 1) list(d) else split(d, d$chunk)
    for (g in groups) {
      cov <- if (g$predictive_coverage_min[1] == g$predictive_coverage_max[1]) {
        sprintf("%.1f%%", g$predictive_coverage_min[1])
      } else {
        sprintf("%.1f-%.1f%%", g$predictive_coverage_min[1], g$predictive_coverage_max[1])
      }
      rows[[length(rows) + 1]] <- data.frame(
        Model = m,
        Chunk = if (nrow(g) > 1) "all" else as.character(g$chunk[1]),
        `Scoring probes present` = sprintf("%s / %s", rp_fmt(g$scoring_probes_present[1]), rp_fmt(g$scoring_probes_total[1])),
        `Background probes present` = sprintf("%s / %s", rp_fmt(g$background_probes_present[1]), rp_fmt(g$background_probes_total[1])),
        `Predictive coverage` = cov,
        `Samples scored` = sprintf("%d / %d", sum(g$samples_scored), sum(g$samples_total)),
        check.names = FALSE, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

rp_page_framingham_table <- function(nums, y) {
  df <- rp_framingham_df(nums$framingham)
  if (is.null(df)) return(invisible(y))
  y <- rp_heading("FraminghamPACE: probes and predictive coverage", y)
  y <- rp_table(df, y, fontsize = 7.5, align = c("l", "l", "r", "r", "r", "r"))
  note <- paste("Scoring probes: probes of the model that are in the data. Background probes: probes of the",
                "reference set that are in the data. Predictive coverage is reported by the FraminghamPACE",
                "package; a sample is scored only if its coverage is at least 75%.")
  rp_lines(strwrap(note, width = 140), y, fontsize = 7, col = RP$grey)
  invisible(NULL)
}

# ---- page: sample QC plots (a 3 x 2 grid) ----
rp_page_qcplots <- function(config, nums, page, n_pages) {
  study <- config$study_name %||% ""
  y <- rp_page("Sample QC plots",
               "One point per sample. Removed samples are not shown. Dashed lines are the thresholds.",
               page, n_pages, study)
  y <- rp_status_legend(y, keep = c("passed", "flagged"))
  d_all <- nums$qc
  d <- d_all[d_all$status != "removed", , drop = FALSE]   # removed samples are not shown
  q <- config$qc_thresholds %||% list()
  panels <- list(
    rp_qc_signal(d, q$signal_intensity_threshold %||% 10.5),
    rp_qc_dot(d, "average_det_pval", "Detection p-value", "Mean detection p-value per sample",
              subtitle = paste0("A sample fails if more than ", 100 * (q$detection_p_fail_fraction %||% 0.01),
                                "% of its probes have a p-value above ", q$detection_p_threshold %||% 0.05),
              plain_labels = TRUE),
    rp_qc_dot(d, "mean_beadcount", "Bead count", "Mean bead count per probe", q$mean_bead_threshold %||% 2),
    rp_qc_dot(d, "snp_outlier_log_odds", "SNP outliers", "SNP outlier score", -(q$snp_outlier_threshold %||% 4)),
    rp_qc_sex(d),
    rp_qc_pca(nums$pca, d_all)
  )
  grid::pushViewport(grid::viewport(x = 0.5, y = y, width = RP_X1 - RP_X0, height = y - 0.06,
                                    just = c("center", "top"), layout = grid::grid.layout(3, 2)))
  for (i in seq_along(panels)) {
    if (is.null(panels[[i]])) next
    print(panels[[i]], vp = grid::viewport(layout.pos.row = (i - 1) %/% 2 + 1, layout.pos.col = (i - 1) %% 2 + 1))
  }
  grid::popViewport()
  invisible(NULL)
}

# ---- page: settings and software ----
rp_page_appendix <- function(config, nums, page, n_pages) {
  study <- config$study_name %||% ""
  y <- rp_page("Settings and software", "What was used for this run", page, n_pages, study)
  y <- rp_heading("Settings", y)
  st <- nums$settings
  sd <- data.frame(Group = ifelse(duplicated(st$group), "", st$group), Setting = st$setting,
                   Value = st$value, stringsAsFactors = FALSE)
  y <- rp_table(sd, y, fontsize = 7.5, wrap = c(Value = 60), max_rows = 45)
  y <- rp_heading("Software", y)
  sw  <- nums$software
  txt <- paste0(sw$item, " ", sw$version)
  n_col <- 3L
  n_row <- ceiling(length(txt) / n_col)
  txt <- c(txt, rep("", n_row * n_col - length(txt)))
  swd <- as.data.frame(matrix(txt, nrow = n_row, ncol = n_col), stringsAsFactors = FALSE)
  colnames(swd) <- c("R and packages", " ", "  ")
  rp_table(swd, y, fontsize = 7.5)
  invisible(NULL)
}

# Runs one page builder; if it fails, the page is replaced by a note and a warning is given
rp_safe_page <- function(label, fun) {
  tryCatch(fun(), error = function(e) {
    warning("Report page '", label, "' could not be created: ", conditionMessage(e), call. = FALSE)
    try(grid::upViewport(0), silent = TRUE)
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
    grid::grid.text(paste0("The page '", label, "' could not be created."), y = 0.52,
                    gp = grid::gpar(fontsize = 12, fontface = "bold", col = RP$red))
    grid::grid.text(paste(strwrap(conditionMessage(e), 100), collapse = "\n"), y = 0.47,
                    gp = grid::gpar(fontsize = 8, col = RP$grey))
  })
}

# The PDF report. figs: the saved plots of step 09 (list of list(plot, width, height)); pages
# for plots that were not made (switched off in the config) are left out.
write_report_pdf <- function(config, nums, chunks, suffix = "", figs = list()) {
  require_package("ggplot2")
  file <- derived_file(config, "report", suffix, "pdf")
  grDevices::pdf(file, width = RP_W, height = RP_H, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  pages <- list()
  add_page <- function(label, shows, fun) {
    pages[[length(pages) + 1]] <<- list(label = label, shows = shows, fun = fun)
  }
  have <- function(keys) any(keys %in% names(figs))
  
  add_page("Run Overview", "", function(i, n, toc) rp_page_overview(config, nums, i, n, toc))
  add_page("Samples", "Samples per chunk, the identity and control checks, all sample QC checks, removed and flagged samples",
           function(i, n, toc) rp_page_samples(config, nums, i, n))
  if (!is.null(nums$qc)) {
    add_page("Sample QC plots", "Signal intensity, detection p-value, bead count, SNP outliers, sex check, control-probe PCA",
             function(i, n, toc) rp_page_qcplots(config, nums, i, n))
  }
  add_page("Probes", "Probes removed by each filter, in each chunk and in all chunks together",
           function(i, n, toc) rp_page_probes(config, nums, i, n))
  if (have(c("beta", "M"))) {
    add_page("Normalization", "Distribution of the BMIQ beta values and the M values", function(i, n, toc) {
      rp_page_figures(config, figs, c("beta", "M"), "Normalization",
                      "BMIQ-normalised beta values and M values of every sample", i, n)
    })
  }
  cc <- grep("^cell_counts_", names(figs), value = TRUE)
  if (length(cc)) {
    add_page("Cell composition", "Estimated proportions of the cell types", function(i, n, toc) {
      rp_page_figures(config, figs, utils::head(cc, 2), "Cell composition",
                      "Estimated proportions of cell types in each sample", i, n)
    })
  }
  if (have("clocks_vs_age")) {
    add_page("Epigenetic clocks", "PC clocks against chronological age", function(i, n, toc) {
      rp_page_figures(config, figs, "clocks_vs_age", "Epigenetic clocks",
                      "PC clocks against chronological age", i, n)
    })
  }
  if (have("systemsage_vs_age")) {
    add_page("SystemsAge", "Scores of the SystemsAge systems against chronological age", function(i, n, toc) {
      rp_page_figures(config, figs, "systemsage_vs_age", "SystemsAge",
                      "Scores of the SystemsAge systems against chronological age", i, n)
    })
  }
  if (have(c("pace_distribution", "pace_by_sex"))) {
    add_page("Pace of aging", "DunedinPACE and FraminghamPACE scores, and the FraminghamPACE probe coverage",
             function(i, n, toc) {
               rp_page_figures(config, figs, c("pace_distribution", "pace_by_sex"), "Pace of aging",
                               "DunedinPACE and FraminghamPACE scores", i, n,
                               reserve = if (!is.null(nums$framingham)) 0.20 else 0,
                               extra = function(y) rp_page_framingham_table(nums, y))
             })
  }
  add_page("Settings and software", "The settings of the run and the versions of R and the main packages",
           function(i, n, toc) rp_page_appendix(config, nums, i, n))
  
  n_pages <- length(pages)
  toc <- data.frame(Page = seq_len(n_pages)[-1],
                    Section = vapply(pages[-1], function(p) p$label, character(1)),
                    `What it shows` = vapply(pages[-1], function(p) p$shows, character(1)),
                    check.names = FALSE, stringsAsFactors = FALSE)
  for (i in seq_len(n_pages)) {
    rp_safe_page(pages[[i]]$label, function() pages[[i]]$fun(i, n_pages, toc))
  }
  message("PDF report written: ", file, " (", n_pages, " pages)")
  invisible(file)
}