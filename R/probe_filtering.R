# ------------------------------------------------------------------------------
# Probe filtering. Works for any array that detect_array_type() recognises
# (EPICv1, EPICv2, 450k). Every filter prints how many probes it removed and
# why; the same numbers are saved in the attrition table (CSV).
#
# run_probe_filtering() is the entry point. The probes are filtered chunk by chunk:
# every chunk is counted, judged and filtered with its own samples, so a chunk can be
# run alone (--chunk=N) and the chunks are independent of each other.
#
#   1. count, for every probe, the samples of the chunk that fail the criteria
#   2. decide which probes fail ("in more than X% of the samples of the chunk")
#   3. apply all filters and save
# ------------------------------------------------------------------------------

read_probe_vector <- function(file, column = NULL, header = TRUE) {
  if (!file.exists(file)) {
    warning("Probe list not found: ", file, call. = FALSE)
    return(character())
  }
  if (grepl("[.]txt$", file, ignore.case = TRUE) && is.null(column)) {
    return(unique(trimws(readLines(file, warn = FALSE))))
  }
  x <- utils::read.csv(file, stringsAsFactors = FALSE, header = header, check.names = FALSE)
  if (!nrow(x)) {
    return(character())
  }
  if (!is.null(column) && column %in% colnames(x)) {
    return(unique(as.character(x[[column]])))
  }
  unique(as.character(x[[1]]))
}

# array_type: "epicv1", "epicv2" or "450k" (from detect_array_type(), i.e. detected from the data)
# Returns the unique probe IDs from the lists.
cross_reactive_probe_ids <- function(config, array_type) {
  opts <- config$probe_filtering %||% list()
  files <- opts$cross_reactive_files %||% list()
  probe_dir <- config_path(config, config$paths$probe_list_dir %||% "data/probe_lists")
  
  # EPICv1 and EPICv2 share the same lists (plain cg IDs)
  if (array_type %in% c("epicv1", "epicv2")) {
    lists <- list(
      list(file = files$epic_pidsley %||% "13059_2016_1066_MOESM1_ESM.csv"),
      list(file = files$epic_mccartney_cpg %||% "1-s2.0-S221359601630071X-mmc2.txt", header = FALSE),
      list(file = files$epic_mccartney_noncpg %||% "1-s2.0-S221359601630071X-mmc3.txt", header = FALSE)
    )
  } else if (array_type == "450k") {
    lists <- list(
      list(file = files$k450_chen %||% "48639-non-specific-probes-Illumina450k.csv", column = "TargetID"),
      list(file = files$k450_benton %||% "HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt", header = FALSE)
    )
  } else {
    warning("No cross-reactive probe lists are defined for array type '", array_type,
            "'; this filter is skipped.", call. = FALSE)
    return(character())
  }
  
  ids <- character()
  for (l in lists) {
    ids <- c(ids, read_probe_vector(file.path(probe_dir, l$file), column = l$column, header = l$header %||% TRUE))
  }
  unique(ids)
}

# ---------- STEP 1: count failing samples per probe ----------------------------
# Returns, for every probe of the chunk, the number of its samples with a detection p
# above the threshold and the number with a bead count below the minimum.
probe_qc_counts_chunk <- function(config, chunk) {
  require_package("minfi")
  opts     <- config$probe_filtering %||% list()
  bead_min <- opts$bead_min_count %||% 3
  det_thr  <- opts$detection_p_threshold %||% 0.05
  
  grset <- readRDS(noob_grset_file(config, chunk))
  probes <- rownames(grset)
  grset_samples <- colnames(grset)
  rm(grset)
  
  # judge probes only on the samples that are in the analysis
  rgset  <- readRDS(good_rgset_file(config, chunk))
  common <- intersect(grset_samples, colnames(rgset))
  if (!length(common)) stop("No shared sample names between GRset and RGset for chunk ", chunk)
  message("Probe QC counts (chunk ", chunk, "): ", length(common), " samples (GRset ",
          length(grset_samples), ", RGset ", ncol(rgset), ")")
  rgset <- rgset[, common]
  
  det_p  <- minfi::detectionP(rgset)
  nbeads <- extract_nbeads_matrix(rgset)
  rm(rgset); gc()
  
  # probes whose row names are not in these matrices cannot be judged (counted as NA)
  if (sum(probes %in% rownames(det_p)) < 0.5 * length(probes)) {
    warning("Detection p-value matrix matches only ", sum(probes %in% rownames(det_p)), " of ",
            length(probes), " probes; check that its rownames match the GRset.", call. = FALSE)
  }
  if (sum(probes %in% rownames(nbeads)) < 0.5 * length(probes)) {
    warning("Bead count matrix matches only ", sum(probes %in% rownames(nbeads)), " of ",
            length(probes), " probes; check that its rownames match the GRset ",
            "(EPICv2 rownames carry suffixes such as cgXXXX_TC11).", call. = FALSE)
  }
  bead_ok <- !all(is.na(nbeads))
  if (!bead_ok) {
    warning("Bead counts are all missing; the bead-count filter will be skipped.", call. = FALSE)
  }
  
  list(
    chunk = chunk,
    samples = common,
    bead_counts_available = bead_ok,
    counts = data.frame(
      probe = probes,
      n_det_fail = unname(rowSums(is.na(det_p) | det_p > det_thr)[probes]),
      n_bead_fail = unname(rowSums(is.na(nbeads) | nbeads < bead_min)[probes]),
      stringsAsFactors = FALSE
    )
  )
}

# ---------- STEP 2: decide which probes fail in this chunk ---------------------
# counts: the step 1 result of one chunk
probe_qc_decide <- function(config, counts) {
  opts <- config$probe_filtering %||% list()
  det_thr       <- opts$detection_p_threshold %||% 0.05
  bead_min      <- opts$bead_min_count %||% 3
  det_frac_max  <- opts$detection_p_fail_fraction %||% 0.01
  bead_frac_max <- opts$bead_fail_fraction %||% 0.05
  
  probes    <- counts$counts$probe
  n_samples <- length(counts$samples)
  det_fail  <- counts$counts$n_det_fail
  bead_fail <- counts$counts$n_bead_fail
  bead_ok   <- counts$bead_counts_available
  
  # probes that could not be judged (NA) fail
  good_detection <- probes[!is.na(det_fail) & det_fail / n_samples <= det_frac_max]
  good_beadcount <- if (bead_ok) probes[!is.na(bead_fail) & bead_fail / n_samples <= bead_frac_max] else probes
  
  message("Probe QC (chunk ", counts$chunk, "): ", n_samples, " samples")
  message("  detection p > ", det_thr, " in more than ", format(100 * det_frac_max), "% of the samples: ",
          format(length(probes) - length(good_detection), big.mark = ","), " of ",
          format(length(probes), big.mark = ","), " probes (before the other filters)")
  if (det_frac_max * n_samples < 1) {
    message("  (with ", n_samples, " samples, one failing sample is already more than ",
            format(100 * det_frac_max), "%)")
  }
  if (bead_ok) {
    message("  bead count < ", bead_min, " in more than ", format(100 * bead_frac_max), "% of the samples: ",
            format(length(probes) - length(good_beadcount), big.mark = ","), " of ",
            format(length(probes), big.mark = ","), " probes (before the other filters)")
    if (bead_frac_max * n_samples < 1) {
      message("  (with ", n_samples, " samples, one failing sample is already more than ",
              format(100 * bead_frac_max), "%)")
    }
  } else {
    message("  bead counts are missing: the bead-count filter will be skipped")
  }
  
  list(
    n_samples = n_samples,
    bead_counts_available = bead_ok,
    good_detection = good_detection,
    good_beadcount = good_beadcount
  )
}

# Saves the attrition table (one row per filter, with the reason) as a CSV.
# There is no figure.
record_probe_attrition <- function(config, chunk, attrition) {
  out_csv <- study_file(config, "probe_filter_attrition", chunk = chunk, ext = "csv")
  utils::write.csv(attrition, out_csv, row.names = FALSE)
  invisible(out_csv)
}

# Adds one row to the attrition table and prints what the filter did and why.
# The % is relative to the probes at the start (first row of the table), so the
# percentages of all filters add up to the total removed.
add_attrition_step <- function(attrition, filter_name, before, after, reason = "") {
  removed <- before - after
  start   <- attrition$probes_remaining[1]
  pct     <- if (start > 0) 100 * removed / start else 0
  message(sprintf("  %-46s removed %9s (%5.2f%%) | remaining %9s | %s",
                  filter_name, format(removed, big.mark = ","), pct,
                  format(after, big.mark = ","), reason))
  if (after == 0) stop("Filter '", filter_name, "' removed all probes", call. = FALSE)
  rbind(
    attrition,
    data.frame(
      filter = filter_name,
      probes_remaining = after,
      probes_removed = removed,
      percent_of_start = round(pct, 2),
      reason = reason,
      stringsAsFactors = FALSE
    )
  )
}

# ---------- STEP 3: apply all filters to the chunk -----------------------------
# global:  the decision of probe_qc_decide() for this chunk (NULL if bead count and
#          detection p filtering are both off)
# samples: samples counted in step 1 (NULL: use all samples of the GRset)
probe_filter_chunk <- function(config, chunk, global = NULL, samples = NULL) {
  require_package("minfi")
  opts  <- config$probe_filtering %||% list()
  grset <- readRDS(noob_grset_file(config, chunk))
  
  # Array type is detected from the data (nothing to set in the config)
  array_type <- detect_array_type(grset)
  
  # the samples counted in step 1 are the samples in the analysis
  if (!is.null(samples)) {
    common <- intersect(colnames(grset), samples)
    if (!length(common)) stop("No shared sample names between GRset and the counted samples for chunk ", chunk)
    n_drop <- ncol(grset) - length(common)
    if (n_drop > 0) {
      message("  ", n_drop, " sample(s) in the Noob GRset are not in the good-sample set and are dropped")
    }
    grset <- grset[, common]
  }
  message("Probe filtering: ", ncol(grset), " samples")
  
  # The probes that each filter removes are recorded (by name), so that the probes removed in
  # all chunks can be counted once in the report: <study>_probe_filter_removed_chunkN.csv
  removed_ids <- list()
  prev_ids    <- rownames(grset)
  record_and_add <- function(attrition, filter_name, before, after, reason = "") {
    out <- add_attrition_step(attrition, filter_name, before, after, reason)
    now <- rownames(grset)
    removed_ids[[filter_name]] <<- setdiff(prev_ids, now)
    prev_ids <<- now
    out
  }
  
  n_start <- nrow(grset)   # all probes of the Noob-normalized data
  message("Probe filtering (chunk ", chunk, ", ", array_type, "): ",
          format(n_start, big.mark = ","), " probes at the start")
  message("  (% = share of the ", format(n_start, big.mark = ","), " probes at the start)")
  attrition <- data.frame(
    filter = "Starting probes",
    probes_remaining = n_start,
    probes_removed = 0,
    percent_of_start = 0,
    reason = "probes in the Noob-normalized data",
    stringsAsFactors = FALSE
  )
  
  # minfi::mapToGenome() drops the probes that have no genomic position in the annotation
  grset <- minfi::mapToGenome(grset)
  attrition <- record_and_add(attrition, "Remove probes without a genomic position", n_start, nrow(grset),
                              "no genomic position in the array annotation (dropped by minfi::mapToGenome)")
  
  if (isTRUE(opts$remove_snp_probes)) {
    before <- nrow(grset)
    # dropLociWithSnps needs SNP columns in the array annotation; check they are filled
    si <- as.data.frame(minfi::getSnpInfo(grset))
    snp_cols <- intersect(c("CpG_rs", "CpG_maf", "SBE_rs", "SBE_maf"), colnames(si))
    if (!length(snp_cols) || all(colSums(!is.na(si[, snp_cols, drop = FALSE])) == 0)) {
      warning("Array annotation has no SNP information; SNP probe filter skipped. ",
              "Use an external SNP/probe mask instead.", call. = FALSE)
      reason <- "skipped: no SNP information in the array annotation"
    } else {
      grset <- minfi::dropLociWithSnps(grset, snps = c("CpG", "SBE"), maf = 0, snpAnno = NULL)
      reason <- "SNP at the CpG or single-base-extension site (any minor allele frequency, minfi annotation)"
    }
    attrition <- record_and_add(attrition, "Remove CpG/SBE SNP probes", before, nrow(grset), reason)
  }
  
  # The probe lists combine cross-reactive and polymorphic probes, so these two
  # settings are one filter: it runs if either one is true.
  if (isTRUE(opts$remove_cross_reactive) || isTRUE(opts$remove_polymorphic)) {
    before <- nrow(grset)
    cr_probes <- cross_reactive_probe_ids(config, array_type)
    if (!length(cr_probes)) {
      reason <- "skipped: no cross-reactive/polymorphic probe lists available"
    } else {
      # v2 rownames carry replicate suffixes (cgXXXX_TC11); the lists use plain IDs
      base_id <- sub("_.*$", "", rownames(grset))
      grset <- grset[!base_id %in% cr_probes, ]
      reason <- "cross-reactive or polymorphic probes (from the probe lists)"
    }
    attrition <- record_and_add(attrition, "Remove cross-reactive/polymorphic probes",
                                before, nrow(grset), reason)
  }
  
  if (isTRUE(opts$remove_xy)) {
    before <- nrow(grset)
    annotation <- minfi::getAnnotation(grset)
    chr_col <- intersect(c("chr", "seqnames"), colnames(annotation))
    if (!length(chr_col)) {
      warning("Array annotation has no chromosome column; XY filter skipped.", call. = FALSE)
      reason <- "skipped: no chromosome information in the array annotation"
    } else {
      xy_probes <- rownames(annotation)[as.character(annotation[[chr_col[[1]]]]) %in% c("chrX", "chrY", "X", "Y")]
      grset <- grset[!rownames(grset) %in% xy_probes, ]
      reason <- "probes on chrX or chrY"
    }
    attrition <- record_and_add(attrition, "Remove XY probes", before, nrow(grset), reason)
  } else {
    message("  XY probes are kept (remove_xy: false)")
  }
  
  # Bead count and detection p: decided on the samples of this chunk (step 2)
  if (isTRUE(opts$remove_failed_beadcount)) {
    before <- nrow(grset)
    bead_min  <- opts$bead_min_count %||% 3
    bead_frac <- opts$bead_fail_fraction %||% 0.05
    if (!global$bead_counts_available) {
      reason <- "skipped: bead counts are missing"
    } else {
      grset <- grset[rownames(grset) %in% global$good_beadcount, ]
      reason <- paste0("bead count < ", bead_min, " in more than ", format(100 * bead_frac),
                       "% of the samples of the chunk")
    }
    attrition <- record_and_add(attrition, "Remove beadcount-failed probes", before, nrow(grset), reason)
  }
  
  if (isTRUE(opts$remove_failed_detection)) {
    before <- nrow(grset)
    det_thr  <- opts$detection_p_threshold %||% 0.05
    det_frac <- opts$detection_p_fail_fraction %||% 0.01
    grset <- grset[rownames(grset) %in% global$good_detection, ]
    attrition <- record_and_add(
      attrition, "Remove detection-P-failed probes", before, nrow(grset),
      paste0("detection p > ", det_thr, " in more than ", format(100 * det_frac),
             "% of the samples of the chunk"))
  }
  
  total_removed <- n_start - nrow(grset)
  message(sprintf("Probe filtering (chunk %s, %s) done: %s -> %s probes (%s removed, %.2f%%); %d samples",
                  chunk, array_type, format(n_start, big.mark = ","),
                  format(nrow(grset), big.mark = ","), format(total_removed, big.mark = ","),
                  100 * total_removed / n_start, ncol(grset)))
  
  attrition$array_type <- array_type   # record which array the filtering was run for
  record_probe_attrition(config, chunk, attrition)
  removed_tbl <- do.call(rbind, lapply(names(removed_ids)[lengths(removed_ids) > 0], function(f) {
    data.frame(probe = removed_ids[[f]], filter = f, stringsAsFactors = FALSE)
  }))
  if (is.null(removed_tbl)) removed_tbl <- data.frame(probe = character(), filter = character())
  utils::write.csv(removed_tbl, study_file(config, "probe_filter_removed", chunk = chunk, ext = "csv"),
                   row.names = FALSE)
  saveRDS(grset, probe_filtered_grset_file(config, chunk))
  invisible(probe_filtered_grset_file(config, chunk))
}

# ---------- ENTRY POINT --------------------------------------------------------
# Runs steps 1-3 for each chunk, one after the other. Results are saved per chunk.
run_probe_filtering <- function(config) {
  opts   <- config$probe_filtering %||% list()
  chunks <- get_chunks_to_run(config)
  
  need_counts <- isTRUE(opts$remove_failed_beadcount) || isTRUE(opts$remove_failed_detection)
  
  for (chunk in chunks) {
    message("Probe filtering for chunk ", chunk)
    counts <- NULL
    global <- NULL
    if (need_counts) {
      counts <- probe_qc_counts_chunk(config, chunk)
      global <- probe_qc_decide(config, counts)
    }
    probe_filter_chunk(config, chunk, global = global, samples = counts$samples)
    rm(counts, global)
    gc()
  }
  invisible(NULL)
}