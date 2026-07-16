bmiq_chunk <- function(config, chunk) {
  require_package("minfi")
  require_package("wateRmelon")
  grset <- readRDS(probe_filtered_grset_file(config, chunk))
  beta <- minfi::getBeta(grset)
  annotation <- minfi::getAnnotation(grset)
  if (!"Type" %in% colnames(annotation)) {
    stop("Annotation is missing probe Type required for BMIQ.", call. = FALSE)
  }
  design <- ifelse(annotation[rownames(beta), "Type"] == "I", 1, 2)
  beta_bmiq <- apply(beta, 2, function(x) {
    wateRmelon::BMIQ(x, design.v = design, plots = FALSE)$nbeta
  })
  saveRDS(beta_bmiq, bmiq_betas_file(config, chunk))
  invisible(bmiq_betas_file(config, chunk))
}
