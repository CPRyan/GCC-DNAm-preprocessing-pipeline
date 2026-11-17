# ---------------------------------------------------------------------------
# Array type: detected from the data, nothing to set in the config.
#
# Each stage that needs the array calls detect_array_type() on the object it
# has just loaded (RGChannelSet, MethylSet or GenomicRatioSet), so stages do
# not depend on each other. It logs what it found and returns one of
# "epicv1", "epicv2" or "450k".
# ---------------------------------------------------------------------------

# Display name for messages
array_label <- function(x) c(epicv1 = "EPICv1", epicv2 = "EPICv2", `450k` = "450k")[[x]]

# Order of evidence:
#   1. the array name minfi stores in the object's annotation
#   2. replicate-probe suffixes in the rownames (cgXXXX_TC11, _BC21, ...): EPICv2 only
#   3. probe count of a raw RGChannelSet (EPICv2 has 1,105,209 rows; EPICv1 ~1.05M)
detect_array_type <- function(x, verbose = TRUE) {
  requireNamespace("minfi", quietly = TRUE)
  
  ann <- tryCatch(as.character(BiocGenerics::annotation(x)[["array"]]),
                  error = function(e) NA_character_)
  known <- c(IlluminaHumanMethylationEPICv2 = "epicv2",
             IlluminaHumanMethylationEPIC   = "epicv1",
             IlluminaHumanMethylation450k   = "450k")
  
  type <- NA_character_
  if (length(ann) == 1L && !is.na(ann) && ann %in% names(known)) {
    type <- unname(known[[ann]])
  } else if (any(grepl("_[TB]C[0-9]+$", rownames(x)))) {
    type <- "epicv2"
  } else if (inherits(x, "RGChannelSet")) {
    n <- nrow(x)
    if      (n > 1.08e6) type <- "epicv2"
    else if (n > 1.00e6) type <- "epicv1"
    else if (n > 5.00e5) type <- "450k"
  }
  
  if (is.na(type)) {
    stop("Cannot determine the array type from the data (minfi annotation: '",
         paste(ann, collapse = ""), "', ", nrow(x), " rows). ",
         "Supported arrays: 450k, EPICv1, EPICv2.", call. = FALSE)
  }
  if (verbose) message("[array] ", array_label(type), " (detected from the data)")
  type
}

# Older minfi versions leave the annotation "Unknown" for EPICv2 IDATs, and then
# detectionP / Noob / getBeta cannot find the right manifest. If the annotation is
# not one of the three known arrays, set it from the detected array. If minfi
# already recognised the array, the object is returned unchanged.
ensure_array_annotation <- function(x, array_type) {
  known <- c("IlluminaHumanMethylationEPICv2", "IlluminaHumanMethylationEPIC",
             "IlluminaHumanMethylation450k")
  have <- tryCatch(BiocGenerics::annotation(x), error = function(e) NULL)
  if (!is.null(have) && isTRUE(unname(have["array"]) %in% known)) return(x)
  
  want <- switch(array_type,
                 epicv2 = c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38"),
                 epicv1 = c(array = "IlluminaHumanMethylationEPIC",   annotation = "ilm10b4.hg19"),
                 "450k" = c(array = "IlluminaHumanMethylation450k",   annotation = "ilmn12.hg19"))
  
  tryCatch({
    x <- BiocGenerics::`annotation<-`(x, value = want)
    message("[array] minfi annotation was not set; set to ", want[["array"]])
    x
  }, error = function(e) {
    warning("[array] Could not set the minfi annotation to ", want[["array"]], ": ",
            conditionMessage(e), call. = FALSE)
    x
  })
}