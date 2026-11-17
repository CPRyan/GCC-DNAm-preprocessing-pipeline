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
  
  # force = TRUE: IDATs from the SAME array type can differ slightly in probe
  # count (e.g. EPICv1 chips scanned with different DMAP files); minfi then keeps
  # the shared probes instead of stopping. IDATs from different array types
  # (450k / EPICv1 / EPICv2) still fail, and we say so in plain words.
  rgset <- tryCatch(
    minfi::read.metharray(
      basenames = manifest$Basename,
      extended = extended,
      verbose = TRUE,
      force = TRUE
    ),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("different", msg, fixed = TRUE)) {
        stop(msg, "\nThe IDATs in this chunk seem to come from different array types ",
             "(450k, EPICv1, EPICv2). Split the samples by array type and run them as ",
             "separate studies.", call. = FALSE)
      }
      stop(e)
    }
  )
  colnames(rgset) <- make.unique(as.character(manifest[[sample_id_col]]))
  
  # Log the detected array at the start of the run; set minfi's annotation if
  # an older minfi left it unknown (no change if it was recognised)
  array_type <- detect_array_type(rgset)
  ensure_array_annotation(rgset, array_type)
}

load_and_save_rgset_chunk <- function(config, chunk) {
  manifest <- manifest_for_chunk(config, chunk)
  rgset <- load_rgset_from_manifest(manifest, sample_id_col = "sample_id", extended = TRUE)
  saveRDS(rgset, raw_rgset_file(config, chunk))
  invisible(raw_rgset_file(config, chunk))
}