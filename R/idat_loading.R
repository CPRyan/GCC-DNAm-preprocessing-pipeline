load_rgset_from_manifest <- function(manifest, sample_id_col = "sample_id", extended = TRUE) {
  require_package("minfi")
  if (!"Basename" %in% colnames(manifest)) {
    stop("Manifest must contain Basename.", call. = FALSE)
  }
  if (!sample_id_col %in% colnames(manifest)) {
    stop("Manifest is missing sample ID column: ", sample_id_col, call. = FALSE)
  }
  if (any(!manifest$idat_pair_exists)) {
    stop("Manifest contains samples with missing IDAT pairs.", call. = FALSE)
  }

  rgset <- minfi::read.metharray(
    basenames = manifest$Basename,
    extended = extended,
    verbose = TRUE
  )
  colnames(rgset) <- make.unique(as.character(manifest[[sample_id_col]]))
  rgset
}

load_and_save_rgset_chunk <- function(config, chunk) {
  manifest <- manifest_for_chunk(config, chunk)
  rgset <- load_rgset_from_manifest(manifest, sample_id_col = "sample_id", extended = TRUE)
  saveRDS(rgset, raw_rgset_file(config, chunk))
  invisible(raw_rgset_file(config, chunk))
}
