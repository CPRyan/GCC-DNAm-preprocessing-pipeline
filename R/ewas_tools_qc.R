run_ewas_tools_qc <- function(config) {
  if (!isTRUE(config$ewas_tools_qc$run)) {
    stop("EWAS tools QC is disabled in config, but this pipeline requires it.", call. = FALSE)
  }
  
  manifest <- if (file.exists(chunks_file(config))) {
    read_chunk_assignments(config)
  } else {
    load_and_validate_manifest(config)
  }
  
  qc <- control_metric_qc_from_manifest(
    manifest,
    sex_col = config$sample_sheet$sex_col %||% "Female",
    snp_outlier_threshold = config$qc_thresholds$snp_outlier_threshold %||% 4
  )
  saveRDS(qc, ewas_qc_file(config))
  qc
}

# Reported sex from the sample sheet as "f", "m" or NA. Only unambiguous codes are read:
# 1 / yes / true / female = female, 0 / no / false / male = male. Letters such as "F" or "M"
# are NOT read, because "F" can mean female or false; those samples count as unknown.
reported_sex_from <- function(x) {
  v   <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(v))
  out[v %in% c("1", "yes", "y", "true", "female")] <- "f"
  out[v %in% c("0", "no", "n", "false", "male")]   <- "m"
  out
}

# Normalised average total intensity of the X and Y chromosome probes, calculated like
# ewastools::check_sex() but accepting the chromosome labels "X" and "chrX" (the EPICv2
# annotation writes "chrX"; check_sex() only finds "X" and then returns NaN for every sample).
sex_intensities <- function(raw) {
  man <- raw$manifest
  chr <- sub("^chr", "", as.character(man$chr))
  x <- man$index[which(chr == "X")]
  y <- man$index[which(chr == "Y")]
  a <- man$index[which(!chr %in% c("X", "Y"))]
  message("  X chromosome probes: ", length(x), "; Y chromosome probes: ", length(y))
  if (!length(x) || !length(y)) {
    stop("no X or Y chromosome probes were found in the ewastools manifest (chromosome labels: ",
         paste(utils::head(unique(as.character(man$chr)), 6), collapse = ", "), ")")
  }
  total <- function(idx) colMeans(raw$M[idx, , drop = FALSE] + raw$U[idx, , drop = FALSE], na.rm = TRUE)
  auto  <- total(a)
  list(X = total(x) / auto, Y = total(y) / auto)
}

# Two sample identity checks of ewastools, run on the data that was read for the control
# metrics:
#  - Sex check: the sex predicted from the X and Y chromosome probes is compared with the
#    reported sex. A mismatch suggests a swapped or mislabelled sample. It uses dye-bias
#    corrected intensities BEFORE masking undetected probes: masking removes most of the
#    Y probes of females, which hides the difference between the sexes (ewastools workflow).
#  - SNP outliers: the SNP probes of the array have three genotype clusters (AA, AB, BB);
#    a high score means the sample fits them badly, which suggests contamination or a
#    failed sample. Samples with a score above -snp_outlier_threshold (default -4, the
#    cut-off recommended by the ewastools authors) are flagged. Detection p-values, masking
#    and dye-bias correction come first, as in the ewastools workflow.
# Each check runs on its own. If one cannot be run, it is left out with a warning and no
# sample is flagged because of it.
# Returns list(pass = <TRUE = passed>, values = <the measured values>), one entry per
# metric, each with one value per sample.
identity_checks <- function(data, reported_sex, snp_outlier_threshold = 4) {
  n      <- length(reported_sex)
  pass   <- list()
  values <- list()
  
  # ---- sex check ----
  tryCatch({
    known_m <- which(reported_sex == "m")
    known_f <- which(reported_sex == "f")
    if (!length(known_m) || !length(known_f)) {
      stop("it needs reported males and females; found ", length(known_m), " male and ",
           length(known_f), " female sample(s)")
    }
    xy <- sex_intensities(ewastools::correct_dye_bias(data))   # not masked
    values$sex_x_intensity <- xy$X   # kept in the QC metrics table, to plot Y against X
    values$sex_y_intensity <- xy$Y
    message("  normalised X / Y intensity (median): females ",
            sprintf("%.2f / %.2f", stats::median(xy$X[known_f], na.rm = TRUE), stats::median(xy$Y[known_f], na.rm = TRUE)),
            "; males ",
            sprintf("%.2f / %.2f", stats::median(xy$X[known_m], na.rm = TRUE), stats::median(xy$Y[known_m], na.rm = TRUE)))
    predicted <- as.character(ewastools::predict_sex(xy$X, xy$Y, known_m, known_f))
    status <- ifelse(is.na(reported_sex), "reported sex missing",
                     ifelse(is.na(predicted), "inconclusive",
                            ifelse(predicted == reported_sex, "match", "mismatch")))
    tab <- table(status)
    message("Sex check (ewastools): ", paste(names(tab), tab, sep = " = ", collapse = ", "))
    if (mean(status == "inconclusive") > 0.2) {
      stop("the sex could not be predicted for more than 20% of the samples, so the check is not reliable")
    }
    if (mean(status == "mismatch") > 0.5) {
      stop("more than half of the samples have a sex mismatch; the coding of the sex column is probably inverted")
    }
    # only a proven mismatch fails; samples with a missing or inconclusive sex are counted
    # separately in the log
    pass$passed_sex_mismatch <- status != "mismatch"
    values$reported_sex  <- reported_sex
    values$predicted_sex <- predicted
    values$sex_check     <- status
  }, error = function(e) {
    warning("The sex check could not be run: ", conditionMessage(e), call. = FALSE)
  })
  invisible(gc())
  
  # ---- SNP outliers ----
  tryCatch({
    prepared <- ewastools::correct_dye_bias(ewastools::mask(ewastools::detectionP(data), 0.01))
    beta <- ewastools::dont_normalize(prepared)
    man  <- prepared$manifest
    idx  <- man$index[man$probe_type == "rs"]
    if (!length(idx)) stop("there are no SNP (rs) probes in the ewastools manifest of this array")
    genotypes <- ewastools::call_genotypes(beta[idx, , drop = FALSE], learn = FALSE)
    score <- as.numeric(ewastools::snp_outliers(genotypes))
    if (length(score) != n || all(is.na(score))) stop("the outlier scores could not be calculated")
    pass$passed_snp_outlier_bl4 <- score <= -snp_outlier_threshold
    values$snp_outlier_log_odds <- score
    message("SNP outliers (ewastools): ", length(idx), " SNP probes; ",
            sum(!pass$passed_snp_outlier_bl4, na.rm = TRUE), " of ", n,
            " sample(s) above the cut-off of ", -snp_outlier_threshold)
  }, error = function(e) {
    warning("The SNP outlier check could not be run: ", conditionMessage(e), call. = FALSE)
  })
  
  list(pass = pass, values = values)
}

control_metric_qc_from_manifest <- function(manifest, sex_col = "Female", snp_outlier_threshold = 4) {
  require_package("ewastools")
  
  if (!"Basename" %in% colnames(manifest)) {
    stop("Manifest must contain Basename for EWAS tools QC.", call. = FALSE)
  }
  if (!"sample_id" %in% colnames(manifest)) {
    stop("Manifest must contain sample_id for EWAS tools QC.", call. = FALSE)
  }
  
  basenames <- manifest$Basename
  sample_ids <- manifest$sample_id
  
  message("Reading the IDATs of ", length(basenames), " samples (EWAS tools) ...")
  invisible(utils::capture.output(
    data_for_control_metrics <- ewastools::read_idats(basenames, quiet = TRUE)
  ))
  message("  done")
  control_test_values <- ewastools::control_metrics(data_for_control_metrics)
  
  passed_qc <- data.frame(matrix(ncol = length(control_test_values), nrow = length(sample_ids)))
  qc_result_values <- data.frame(matrix(ncol = length(control_test_values), nrow = length(sample_ids)))
  rownames(passed_qc) <- sample_ids
  rownames(qc_result_values) <- sample_ids
  colnames(passed_qc) <- names(control_test_values)
  colnames(qc_result_values) <- names(control_test_values)
  
  for (metric_name in names(control_test_values)) {
    metric_value <- as.numeric(control_test_values[[metric_name]])
    threshold <- attr(control_test_values[[metric_name]], "threshold")
    if (length(metric_value) != length(sample_ids)) {
      stop(
        "EWAS tools metric length did not match manifest sample count for metric: ",
        metric_name,
        call. = FALSE
      )
    }
    passed_qc[[metric_name]] <- metric_value >= threshold
    qc_result_values[[metric_name]] <- metric_value
    attr(qc_result_values[[metric_name]], "threshold") <- threshold
  }
  
  # sample identity checks: SNP outliers and sex
  if (!is.null(sex_col) && sex_col %in% colnames(manifest)) {
    reported_sex <- reported_sex_from(manifest[[sex_col]])
    message("Reported sex (column '", sex_col, "'): ", sum(reported_sex == "f", na.rm = TRUE),
            " female, ", sum(reported_sex == "m", na.rm = TRUE), " male, ",
            sum(is.na(reported_sex)), " unknown")
  } else {
    message("The sample sheet has no sex column '", sex_col, "': the sex check is skipped.")
    reported_sex <- rep(NA_character_, length(sample_ids))
  }
  checks <- identity_checks(data_for_control_metrics, reported_sex, snp_outlier_threshold)
  for (name in names(checks$pass))   passed_qc[[name]]        <- checks$pass[[name]]
  for (name in names(checks$values)) qc_result_values[[name]] <- checks$values[[name]]
  
  list(
    passed_QC = passed_qc,
    QC_result_values = qc_result_values
  )
}

# Uses the saved EWAS tools QC if it exists; otherwise calculates it (slow: reads every IDAT) and saves it.
load_or_run_ewas_tools_qc <- function(config) {
  file <- ewas_qc_file(config)
  if (file.exists(file)) {
    return(readRDS(file))
  }
  run_ewas_tools_qc(config)
}

as_ewas_table <- function(x, sample_id_col = "sample_id") {
  if (is.null(x)) {
    return(NULL)
  }
  out <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  out[[sample_id_col]] <- rownames(out)
  rownames(out) <- NULL
  out
}