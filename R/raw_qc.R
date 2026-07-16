ensure_manifest_package <- function(rgset) {
  require_package("minfi")
  array_name <- minfi::annotation(rgset)[["array"]]
  if (is.null(array_name) || is.na(array_name) || !nzchar(array_name)) {
    return(invisible(TRUE))
  }
  manifest_package <- paste0(array_name, "manifest")
  if (!requireNamespace(manifest_package, quietly = TRUE)) {
    stop("Required array manifest package is not installed: ", manifest_package, call. = FALSE)
  }
  invisible(TRUE)
}

calculate_control_probe_pca <- function(rgset, variance_threshold = 0.90, pseudocount = 100) {
  require_package("minfi")
  ensure_manifest_package(rgset)
  control_probes <- minfi::getControlAddress(rgset)
  control_probes <- intersect(control_probes, rownames(minfi::getRed(rgset)))
  if (!length(control_probes)) {
    stop("No control probes were found in the RGset.", call. = FALSE)
  }

  red <- minfi::getRed(rgset)[control_probes, , drop = FALSE]
  green <- minfi::getGreen(rgset)[control_probes, , drop = FALSE]
  control_beta <- red / (red + green + pseudocount)
  keep_probe <- apply(control_beta, 1, function(x) all(is.finite(x)) && stats::sd(x) > 0)
  control_beta <- control_beta[keep_probe, , drop = FALSE]
  pca_result <- stats::prcomp(t(control_beta), center = TRUE, scale. = TRUE)
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cumulative_variance_explained <- cumsum(variance_explained)
  cutoff <- which(cumulative_variance_explained >= variance_threshold)[[1]]

  list(
    scores = data.frame(sample_id = rownames(pca_result$x), pca_result$x[, seq_len(cutoff), drop = FALSE], check.names = FALSE),
    variance = data.frame(
      pc = paste0("PC", seq_along(variance_explained)),
      variance_explained = variance_explained,
      cumulative_variance_explained = cumulative_variance_explained
    ),
    pca = pca_result,
    n_control_probes_used = nrow(control_beta),
    n_pcs_to_threshold = cutoff
  )
}

extract_nbeads_matrix <- function(rgset) {
  require_package("Biobase")
  if ("getNBeads" %in% getNamespaceExports("minfi")) {
    nbeads <- try(minfi::getNBeads(rgset), silent = TRUE)
    if (!inherits(nbeads, "try-error")) {
      return(nbeads)
    }
  }
  assay_elements <- Biobase::assayDataElementNames(rgset)
  if ("NBeads" %in% assay_elements) {
    return(Biobase::assayDataElement(rgset, "NBeads"))
  }
  if (requireNamespace("wateRmelon", quietly = TRUE)) {
    return(wateRmelon::beadcount(rgset))
  }
  stop("Could not extract bead counts. Load IDATs with extended = TRUE, or install wateRmelon.", call. = FALSE)
}

calculate_bead_metrics <- function(rgset, bead_min_count = 5, bead_fail_fraction = 0.05, mean_bead_threshold = 2) {
  nbeads <- extract_nbeads_matrix(rgset)
  mean_bead_count <- colMeans(nbeads, na.rm = TRUE)
  fraction_low <- colMeans(nbeads < bead_min_count, na.rm = TRUE)
  data.frame(
    sample_id = names(mean_bead_count),
    mean_beadcount = as.numeric(mean_bead_count),
    fraction_probes_low_beadcount = as.numeric(fraction_low),
    bead_min_count = bead_min_count,
    bead_fail_fraction = bead_fail_fraction,
    mean_bead_threshold = mean_bead_threshold,
    passed_beadcount = ifelse(fraction_low > bead_fail_fraction | mean_bead_count <= mean_bead_threshold, "fail", "pass"),
    stringsAsFactors = FALSE
  )
}

calculate_detection_metrics <- function(rgset, detection_p_threshold = 0.05, detection_p_fail_fraction = 0.01) {
  require_package("minfi")
  ensure_manifest_package(rgset)
  det_p <- minfi::detectionP(rgset)
  average_det_pval <- colMeans(det_p, na.rm = TRUE)
  fraction_failed <- colMeans(det_p > detection_p_threshold, na.rm = TRUE)
  data.frame(
    sample_id = names(average_det_pval),
    average_det_pval = as.numeric(average_det_pval),
    fraction_detection_p_failed = as.numeric(fraction_failed),
    detection_p_threshold = detection_p_threshold,
    detection_p_fail_fraction = detection_p_fail_fraction,
    passed_detection_cutoff = ifelse(fraction_failed > detection_p_fail_fraction, "fail", "pass"),
    stringsAsFactors = FALSE
  )
}

calculate_signal_metrics <- function(rgset, signal_intensity_threshold = 10.5) {
  require_package("minfi")
  mset <- minfi::preprocessRaw(rgset)
  meth <- minfi::getMeth(mset)
  unmeth <- minfi::getUnmeth(mset)
  mqc <- log2(apply(meth, 2, stats::median, na.rm = TRUE))
  uqc <- log2(apply(unmeth, 2, stats::median, na.rm = TRUE))
  data.frame(
    sample_id = names(mqc),
    MQC = as.numeric(mqc),
    UQC = as.numeric(uqc),
    signal_intensity_threshold = signal_intensity_threshold,
    passed_signal_quality = ifelse(mqc < signal_intensity_threshold | uqc < signal_intensity_threshold, "fail", "pass"),
    stringsAsFactors = FALSE
  )
}

run_raw_qc_for_rgset <- function(config, rgset) {
  thresholds <- config$qc_thresholds
  signal <- calculate_signal_metrics(rgset, thresholds$signal_intensity_threshold %||% 10.5)
  beads <- calculate_bead_metrics(
    rgset,
    bead_min_count = thresholds$bead_min_count %||% 5,
    bead_fail_fraction = thresholds$bead_fail_fraction %||% 0.05,
    mean_bead_threshold = thresholds$mean_bead_threshold %||% 2
  )
  detection <- calculate_detection_metrics(
    rgset,
    detection_p_threshold = thresholds$detection_p_threshold %||% 0.05,
    detection_p_fail_fraction = thresholds$detection_p_fail_fraction %||% 0.01
  )
  control_pca <- calculate_control_probe_pca(rgset, thresholds$control_pca_variance_threshold %||% 0.90)

  list(signal = signal, beads = beads, detection = detection, control_pca = control_pca)
}
