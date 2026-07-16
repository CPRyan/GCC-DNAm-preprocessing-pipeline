source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "_common.R"))

manifest <- load_and_validate_manifest(config)
utils::write.csv(manifest, validated_manifest_file(config), row.names = FALSE)

missing <- manifest[!manifest$idat_pair_exists, , drop = FALSE]
if (nrow(missing)) {
  missing_file <- study_file(config, "missing_idats", ext = "csv")
  utils::write.csv(missing, missing_file, row.names = FALSE)
  stop("Missing IDAT pairs for ", nrow(missing), " sample(s). Details: ", missing_file, call. = FALSE)
}

message("Validated ", nrow(manifest), " samples.")
message("Wrote: ", validated_manifest_file(config))
