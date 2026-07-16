noob_normalize_chunk <- function(config, chunk) {
  require_package("minfi")
  rgset <- readRDS(good_rgset_file(config, chunk))
  grset <- minfi::preprocessNoob(rgset)
  saveRDS(grset, noob_grset_file(config, chunk))
  invisible(noob_grset_file(config, chunk))
}
