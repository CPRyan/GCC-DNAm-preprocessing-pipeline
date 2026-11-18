# BMIQ is an intra-sample normalisation: each sample is normalised on its own, so
# chunking does not change the results (unlike the cohort-wide probe criteria).
bmiq_chunk <- function(config, chunk) {
  require_package("minfi")
  require_package("wateRmelon")
  require_package("RPMM")   # wateRmelon::BMIQ stops without it
  
  n_cores <- as.integer(config$bmiq$n_cores %||% 1)
  seed    <- as.integer(config$chunking$seed %||% 123)
  n_chunks <- as.integer(config$chunking$n_chunks %||% 1)
  name_col <- name_column(config)
  names_of <- sample_name_lookup(config)   # sample_id -> name, shown next to the IDATid
  verbose <- isTRUE(config$bmiq$verbose)   # TRUE: show BMIQ's own progress output
  if (n_cores > 1 && .Platform$OS.type == "windows") {
    warning("Parallel BMIQ is not available on Windows; using 1 core.", call. = FALSE)
    n_cores <- 1L
  }
  
  grset <- readRDS(probe_filtered_grset_file(config, chunk))
  beta <- minfi::getBeta(grset)
  annotation <- minfi::getAnnotation(grset)
  if (!"Type" %in% colnames(annotation)) {
    stop("Annotation is missing probe Type required for BMIQ.", call. = FALSE)
  }
  type <- as.character(annotation[rownames(beta), "Type"])
  if (!all(type %in% c("I", "II"))) {
    stop("Probe Type must be I or II for every probe; found: ",
         paste(unique(type[!type %in% c("I", "II")]), collapse = ", "), call. = FALSE)
  }
  design <- ifelse(type == "I", 1, 2)
  rm(grset, annotation)
  
  message("BMIQ (chunk ", chunk, "): ", ncol(beta), " samples, ", format(nrow(beta), big.mark = ","),
          " probes (", format(sum(design == 1), big.mark = ","), " type I, ",
          format(sum(design == 2), big.mark = ","), " type II), ",
          round(100 * mean(is.na(beta)), 2), "% missing beta values, ", n_cores, " core(s)")
  
  # BMIQ picks probes at random for its model fit, so each sample gets its own seed
  normalize_sample <- function(i) {
    set.seed(seed + i)
    # BMIQ labels its output with sampleID (default 1); use the position and the real ID
    nm <- if (!is.null(names_of)) unname(names_of[colnames(beta)[i]]) else NA_character_
    sample_label <- paste0(i, "/", ncol(beta), " ", colnames(beta)[i],
                           if (!is.na(nm)) paste0(" (", name_col, " ", nm, ")"))
    label <- paste0("chunk ", chunk, "/", n_chunks, " - sample ", sample_label)
    # header printed before BMIQ's own output (same stream, so the order is kept)
    if (verbose) cat("\n--- BMIQ - ", label, " ---\n", sep = "")
    # R's lgamma() warns "underflow occurred in 'lgammacor'" when a shape parameter of the
    # beta mixture fit gets huge. It is a harmless numerical note, but it can be repeated
    # thousands of times. It is not hidden: the log records it with a note, and here it is
    # counted per sample, so that the samples with many of them can be found.
    n_underflow <- 0L
    tryCatch({
      nbeta <- withCallingHandlers(
        wateRmelon::BMIQ(beta[, i], design.v = design, plots = FALSE,
                         pri = verbose, sampleID = sample_label)$nbeta,
        warning = function(w) {
          if (grepl("underflow occurred in 'lgammacor'", conditionMessage(w), fixed = TRUE)) {
            n_underflow <<- n_underflow + 1L
          }
        }
      )
      if (!verbose || n_underflow > 0) {
        message("  BMIQ done for ", label,
                if (n_underflow > 0) paste0(" (", format(n_underflow, big.mark = ","),
                                            " numerical underflow warnings in the model fit)"))
      }
      nbeta
    }, error = function(e) {
      message("BMIQ failed for ", label, ": ", conditionMessage(e))
      NULL
    })
  }
  
  if (verbose && n_cores > 1) {
    message("Note: with several cores the verbose output of different samples is interleaved; ",
            "use n_cores: 1 for a readable log.")
  }
  
  t0  <- Sys.time()
  res <- parallel::mclapply(seq_len(ncol(beta)), normalize_sample, mc.cores = n_cores)
  
  # a sample that failed gets NA values and is reported
  bad <- vapply(res, function(r) is.null(r) || length(r) != nrow(beta), logical(1))
  if (all(bad)) stop("BMIQ failed for all samples in chunk ", chunk, call. = FALSE)
  if (any(bad)) {
    warning("BMIQ failed for ", sum(bad), " sample(s) in chunk ", chunk, " (values set to NA): ",
            paste(colnames(beta)[bad], collapse = ", "), call. = FALSE)
    res[bad] <- list(rep(NA_real_, nrow(beta)))
  }
  
  beta_bmiq <- do.call(cbind, res)
  dimnames(beta_bmiq) <- dimnames(beta)
  
  message("BMIQ (chunk ", chunk, ") done: ", ncol(beta_bmiq) - sum(bad), " of ", ncol(beta_bmiq),
          " samples normalised in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  
  saveRDS(beta_bmiq, bmiq_betas_file(config, chunk))
  invisible(bmiq_betas_file(config, chunk))
}