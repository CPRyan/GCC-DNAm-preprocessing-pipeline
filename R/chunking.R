create_chunk_assignments <- function(manifest, n_chunks = 1, grouping_cols = NULL, seed = 123) {
  if (n_chunks <= 1) {
    manifest$chunk <- 1L
    return(manifest)
  }

  set.seed(seed)
  if (is.null(grouping_cols) || !length(grouping_cols)) {
    groups <- manifest$sample_id
  } else {
    missing_cols <- setdiff(grouping_cols, colnames(manifest))
    if (length(missing_cols)) {
      stop("Chunk grouping column(s) missing from manifest: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    groups <- do.call(paste, c(manifest[grouping_cols], sep = "||"))
  }

  split_idx <- split(seq_len(nrow(manifest)), groups)
  split_idx <- sample(split_idx)
  chunk_sizes <- integer(n_chunks)
  chunk_assignments <- integer(nrow(manifest))

  for (idx in split_idx) {
    best_chunk <- which.min(chunk_sizes)
    chunk_assignments[idx] <- best_chunk
    chunk_sizes[best_chunk] <- chunk_sizes[best_chunk] + length(idx)
  }

  manifest$chunk <- chunk_assignments
  manifest
}

read_chunk_assignments <- function(config) {
  file <- chunks_file(config)
  if (!file.exists(file)) {
    stop("Chunk assignment file does not exist. Run scripts/01_create_targets_and_chunks.R first: ", file, call. = FALSE)
  }
  utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
}

manifest_for_chunk <- function(config, chunk) {
  manifest <- read_chunk_assignments(config)
  manifest[manifest$chunk == chunk, , drop = FALSE]
}
