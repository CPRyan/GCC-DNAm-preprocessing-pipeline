# ------------------------------------------------------------------------------
# Cell counts: reference in the config, based on access:
#   - FlowSorted.Blood.EPIC          (6 cell types, public)
#   - FlowSorted.BloodExtended.EPIC  (12 cell types, restricted access)
#
# Clocks: PCClocks and SystemsAge (methylCIPHER), DunedinPACE
# ------------------------------------------------------------------------------
# IDOL CpGs are projected (CP/QP) onto reference profiles built from the
# Noob-processed v1 reference library.
#
# Array type (EPICv1 / EPICv2 / 450k): detected from the data, nothing to set
# in the config (see R/array_utils.R). On EPICv2, replicate probes are
# collapsed by lowest detection p-value for the cell counts; clocks /
# DunedinPACE use betas with replicates averaged.
# ------------------------------------------------------------------------------

# ---------- SETTINGS ----------------------------------------------------------
BLOOD_EXTENDED_CELLTYPES <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv",
                              "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")
BLOOD6_CELLTYPES <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")

# One entry per reference. `library` is both the reference name for
# libraryDataGet() and the package holding the IDOL CpG list(s).
CELL_REFERENCES <- list(
  blood6   = list(label        = "Blood6",
                  library      = "FlowSorted.Blood.EPIC",
                  celltypes    = BLOOD6_CELLTYPES,
                  idol_default = "IDOLOptimizedCpGs",
                  idol_450k    = "IDOLOptimizedCpGs450klegacy",
                  out_name     = "Blood6_cell_counts"),
  extended = list(label        = "BloodExtended",
                  library      = "FlowSorted.BloodExtended.EPIC",
                  celltypes    = BLOOD_EXTENDED_CELLTYPES,
                  idol_default = "IDOLOptimizedCpGsBloodExtended",
                  idol_450k    = NULL,
                  out_name     = "BloodExtended_cell_counts")
)

# Read Config value -> "FlowSorted.Blood.EPIC" (default) or "FlowSorted.BloodExtended.EPIC"
normalize_cell_reference <- function(x) {
  if (is.null(x)) {
    message("Cell counts: noob_downstream.cell_reference is not set; ",
            "using FlowSorted.Blood.EPIC (6 cell types).")
    return("blood6")
  }
  x <- tolower(trimws(x))
  if (x %in% c("extended", "flowsorted.bloodextended.epic")) {
    "extended"
  } else if (x %in% c("blood6", "flowsorted.blood.epic")) {
    "blood6"
  } else {
    stop("Unknown noob_downstream.cell_reference: '", x, "'. ",
         "Use FlowSorted.BloodExtended.EPIC or FlowSorted.Blood.EPIC.", call. = FALSE)
  }
}

# Betas with plain cg IDs (replicate probes cgXXXX_TC11, _BC21, ... averaged)
get_betas <- function(grset) {
  betas <- minfi::getBeta(grset)
  ids <- sub("_.*$", "", rownames(betas))
  ok  <- !is.na(betas)
  betas[!ok] <- 0
  out <- rowsum(betas, ids, reorder = FALSE) / rowsum(ok * 1, ids, reorder = FALSE)
  out[is.nan(out)] <- NA
  out
}

# This applies to EPICv2: Per sample, keep the replicate with the lowest detection p-value
# (CpG x sample matrix for the requested plain cg IDs)
collapse_lowest_detp <- function(betas, detp, cpgs) {
  ids <- sub("_.*$", "", rownames(betas))
  sel <- ids %in% cpgs
  b   <- betas[sel, , drop = FALSE]
  d   <- detp[rownames(b), colnames(b), drop = FALSE]
  d[is.na(d)] <- Inf
  by_id <- split(seq_len(nrow(b)), ids[sel])
  cols  <- seq_len(ncol(b))
  out <- t(vapply(by_id, function(rows) {
    # for each sample (column): the replicate with the lowest detection p-value.
    # max.col() works per ROW, so the matrix is transposed: one row per sample.
    best <- if (length(rows) == 1L) rep(1L, length(cols)) else max.col(t(-d[rows, , drop = FALSE]), ties.method = "first")
    b[cbind(rows[best], cols)]
  }, numeric(ncol(b))))
  colnames(out) <- colnames(b)
  out
}

# ---------- BLOOD CELL COMPOSITION --------------------------------------------
# reference: "extended" or "blood6" (see CELL_REFERENCES)
# The code detect array type and handle the deconvolution accordingly using projectCellType_CP()
estimate_cell_counts <- function(config, chunk, reference) {
  ref   <- CELL_REFERENCES[[reference]]
  label <- ref$label
  require_package("minfi")
  require_package("FlowSorted.Blood.EPIC")   # libraryDataGet() and projectCellType_CP()
  require_package(ref$library)
  
  grset <- readRDS(noob_grset_file(config, chunk))
  
  # the array is detected from the data (nothing to set in the config)
  array_type <- detect_array_type(grset)
  
  # IDOL CpG list for this array
  if (array_type == "450k" && is.null(ref$idol_450k)) {
    message(label, " [450k]: no 450k-specific IDOL list; using the EPIC list")
  }
  idol_name <- if (array_type == "450k" && !is.null(ref$idol_450k)) ref$idol_450k else ref$idol_default
  env <- new.env()
  utils::data(list = idol_name, package = ref$library, envir = env)
  idol <- env[[idol_name]]
  if (!is.character(idol)) {
    stop(idol_name, " is not a character vector of CpG IDs; check str(", idol_name, ")")
  }
  
  betas    <- minfi::getBeta(grset)
  has_reps <- anyDuplicated(sub("_.*$", "", rownames(betas))) > 0
  
  if (has_reps) {
    # EPICv2: needs detection p-values to choose among replicate probes
    rgset  <- readRDS(raw_rgset_file(config, chunk))
    # the RGset may hold more samples (or a different order) than the Noob GRset
    common <- intersect(colnames(grset), colnames(rgset))
    if (length(common) == 0) {
      stop("No shared sample names. GRset: ", paste(head(colnames(grset), 3), collapse = ", "),
           " | RGset: ", paste(head(colnames(rgset), 3), collapse = ", "))
    }
    message(label, " [", array_type, "]: ", length(common), " samples in common (GRset ",
            ncol(grset), ", RGset ", ncol(rgset), ")")
    betas <- betas[, common, drop = FALSE]
    detp  <- minfi::detectionP(rgset[, common])
    rm(rgset)
  }
  rm(grset); gc()
  
  # reference profiles (v1 library, Noob)
  refset <- FlowSorted.Blood.EPIC::libraryDataGet(ref$library)
  refb   <- minfi::getBeta(minfi::preprocessNoob(refset))
  ct     <- as.character(refset$CellType)
  rm(refset); gc()
  
  types <- intersect(ref$celltypes, ct)
  if (length(types) < length(ref$celltypes)) {
    stop("Reference ", ref$library, " is missing cell types: ",
         paste(setdiff(ref$celltypes, ct), collapse = ", "))
  }
  
  wanted <- intersect(idol, rownames(refb))
  if (has_reps) {
    y <- collapse_lowest_detp(betas, detp, wanted)
    rm(detp)
  } else {
    y <- betas[rownames(betas) %in% wanted, , drop = FALSE]   # EPICv1 / 450k: IDs are already plain
  }
  rm(betas); gc()
  
  # keep CpGs without missing values
  cpgs <- rownames(y)[rowSums(is.na(y)) == 0 & rowSums(is.na(refb[rownames(y), , drop = FALSE])) == 0]
  y <- y[cpgs, , drop = FALSE]
  message(label, " [", array_type, "]: using ", length(cpgs), " of ", length(idol), " IDOL CpGs")
  if (length(cpgs) < 0.8 * length(idol)) {
    warning(label, " [", array_type, "]: only ", length(cpgs), " of ", length(idol),
            " IDOL CpGs are available; cell proportions may be less reliable.", call. = FALSE)
  }
  
  coef <- sapply(types, function(k) rowMeans(refb[cpgs, ct == k, drop = FALSE]))
  prop <- FlowSorted.Blood.EPIC::projectCellType_CP(y, coef)
  
  message(label, ": ", ncol(y), " samples, ", length(types), " cell types, ",
          length(cpgs), " of ", length(idol), " IDOL CpGs used (", array_type, ")")
  
  prop <- as.data.frame(prop)
  prop <- cbind(sample_id = rownames(prop), prop)
  
  file <- study_file(config, ref$out_name, chunk = chunk, ext = "csv")
  utils::write.csv(prop, file, row.names = FALSE)
  invisible(file)
}

# ---------- PC CLOCKS ---------------------------------------------------------
run_pc_clocks_noob <- function(config, chunk) {
  require_package("minfi")
  require_package("methylCIPHER")
  
  # download PCClocks_data.qs2 (2.1 GB) if missing or incomplete
  pc_file <- config_path(config, config$paths$PCClock_data)
  message("PC clocks data: ", pc_file)
  md5 <- "0f7fa3a89b6559e98478a9e986d28db0"
  if (!file.exists(pc_file) || unname(tools::md5sum(pc_file)) != md5) {
    dir.create(dirname(pc_file), recursive = TRUE, showWarnings = FALSE)
    options(timeout = 4 * 3600)
    download.file("https://zenodo.org/records/19455622/files/PCClocks_data.qs2?download=1",
                  destfile = pc_file, mode = "wb")
    stopifnot(unname(tools::md5sum(pc_file)) == md5)
  }
  
  betas <- get_betas(readRDS(noob_grset_file(config, chunk)))
  message("Input: ", ncol(betas), " samples, ", nrow(betas), " CpGs, ",
          round(100 * mean(is.na(betas)), 2), "% missing beta values")
  manifest <- manifest_for_chunk(config, chunk)
  pheno <- manifest[, c("sample_id", "Age", "Female")]
  pheno$Age <- as.numeric(pheno$Age)
  pheno$Female <- as.numeric(pheno$Female)
  rownames(pheno) <- pheno$sample_id
  
  stopifnot(all(colnames(betas) %in% rownames(pheno)))
  
  pheno <- pheno[colnames(betas), ]
  pc_clocks <- methylCIPHER::calcPCClocks(t(betas), pheno, RData = pc_file)
  message("methylCIPHER ", as.character(packageVersion("methylCIPHER")), ": ",
          length(setdiff(names(pc_clocks), c("sample_id", "Age", "Female"))),
          " PC clock outputs for ", nrow(pc_clocks), " samples")
  
  file <- study_file(config, "PCClocks_noob", chunk = chunk, ext = "csv")
  utils::write.csv(as.data.frame(pc_clocks), file, row.names = FALSE)
  invisible(file)
}

# ---------- SystemsAge --------------------------------------------------------
run_systems_age_noob <- function(config, chunk) {
  require_package("minfi")
  require_package("methylCIPHER")
  
  # download SystemsAge_data.qs2 (3.6 GB) if missing or incomplete
  sa_file <- config_path(config, config$paths$SystemsAge_data)
  md5     <- "c923ca6c1d62391ce5b54662b98a864c"
  if (!file.exists(sa_file) || unname(tools::md5sum(sa_file)) != md5) {
    dir.create(dirname(sa_file), recursive = TRUE, showWarnings = FALSE)
    options(timeout = 8 * 3600)
    download.file("https://zenodo.org/records/19455622/files/SystemsAge_data.qs2?download=1",
                  destfile = sa_file, mode = "wb")
    stopifnot(unname(tools::md5sum(sa_file)) == md5)
  }
  
  betas <- get_betas(readRDS(noob_grset_file(config, chunk)))
  message("Input: ", ncol(betas), " samples, ", nrow(betas), " CpGs, ",
          round(100 * mean(is.na(betas)), 2), "% missing beta values")
  manifest <- manifest_for_chunk(config, chunk)
  pheno <- manifest[, c("sample_id", "Age", "Female")]
  pheno$Age <- as.numeric(pheno$Age)
  pheno$Female <- as.numeric(pheno$Female)
  rownames(pheno) <- pheno$sample_id
  stopifnot(all(colnames(betas) %in% rownames(pheno)))
  pheno <- pheno[colnames(betas), ]
  systems_age <- methylCIPHER::calcSystemsAge(t(betas), pheno, RData = sa_file)
  message("methylCIPHER ", as.character(packageVersion("methylCIPHER")), ": ",
          length(setdiff(names(systems_age), c("sample_id", "Age", "Female"))),
          " SystemsAge outputs for ", nrow(systems_age), " samples")
  
  file <- study_file(config, "SystemsAge_noob", chunk = chunk, ext = "csv")
  utils::write.csv(as.data.frame(systems_age), file, row.names = FALSE)
  invisible(file)
}

# ---------- DunedinPACE -------------------------------------------------------
run_dunedinpace_noob <- function(config, chunk) {
  require_package("DunedinPACE")
  require_package("minfi")
  
  betas <- get_betas(readRDS(noob_grset_file(config, chunk)))
  message("Input: ", ncol(betas), " samples, ", nrow(betas), " CpGs, ",
          round(100 * mean(is.na(betas)), 2), "% missing beta values")
  pace <- DunedinPACE::PACEProjector(betas)
  pace <- as.data.frame(pace)
  message("DunedinPACE ", as.character(packageVersion("DunedinPACE")), ": scores for ", nrow(pace), " samples")
  pace <- cbind(sample_id = rownames(pace), pace)
  
  file <- study_file(config, "DunedinPACE_noob", chunk = chunk, ext = "csv")
  utils::write.csv(pace, file, row.names = FALSE)
  invisible(file)
}

# ---------- RUN NOOB DOWNSTREAM -----------------------------------------------
run_noob_downstream <- function(config, chunk) {
  opts <- config$noob_downstream %||% list()
  
  # a failing module warns instead of stopping the others
  safe <- function(name, expr) {
    message("\n======= Running ", name, " (chunk ", chunk, ") =======")
    tryCatch({
      expr
      message("======= ", name, " completed correctly =======\n")
    }, error = function(e) {
      message("======= ", name, " FAILED: ", conditionMessage(e), " =======\n")
      warning("Chunk ", chunk, ": ", name, " failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
  }
  
  run_cells       <- isTRUE(opts$cell_counts)
  run_PCclocks    <- isTRUE(opts$pc_clocks)
  run_systemsage  <- isTRUE(opts$systems_age)
  run_dunedinpace <- isTRUE(opts$dunedinpace)
  
  # an invalid cell_reference is a config error: stop here rather than warn
  cell_ref <- if (run_cells) normalize_cell_reference(opts$cell_reference) else NULL
  
  if (run_cells)       safe(CELL_REFERENCES[[cell_ref]]$label, estimate_cell_counts(config, chunk, cell_ref))
  if (run_PCclocks)    safe("PC clocks",    run_pc_clocks_noob(config, chunk))
  if (run_systemsage)  safe("SystemsAge",   run_systems_age_noob(config, chunk))
  if (run_dunedinpace) safe("DunedinPACE",  run_dunedinpace_noob(config, chunk))
  invisible(NULL)
}