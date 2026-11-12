# ---------------------------------------------------------------------------
# Noob downstream modules for EPIC v2: cell counts, PC clocks, DunedinPACE
#
# Array type: EPICv2
# Tissue: skeletal muscle (cell conuts and clocks are optimized from blood!)
# Cell counts methods: FlowSorted.Blood.EPIC, EpiDISH (blood), EpiDISH (hiearchical)
# Clocks methods: PC_Clocks, DunedinPACE  
# 
# Nte: replicate probes are collapsed by lowest detection p-value, 
# IDOL CpGs from the v1 library are projected (CP/QP) onto reference profiles 
# built from the Noob-processed v1 reference. 
# Clocks / DunedinPACE use betas with replicates averaged.
# ---------------------------------------------------------------------------

BLOOD_EXTENDED_CELLTYPES <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv",
                              "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")

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

# Per sample, keep the replicate with the lowest detection p-value
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
    best <- if (length(rows) == 1L) 1L else max.col(-d[rows, , drop = FALSE], ties.method = "first")
    b[cbind(rows[best], cols)]
  }, numeric(ncol(b))))
  colnames(out) <- colnames(b)
  out
}

# ---------- FlowSorted.BloodExtended.EPIC CELL COMPOSITION --------------------
estimate_cell_counts_blood_extended <- function(config, chunk) {
  require_package("minfi")
  require_package("FlowSorted.Blood.EPIC")
  require_package("FlowSorted.BloodExtended.EPIC")
  utils::data("IDOLOptimizedCpGsBloodExtended", package = "FlowSorted.BloodExtended.EPIC")
  
  grset <- readRDS(noob_grset_file(config, chunk))
  rgset <- readRDS(raw_rgset_file(config, chunk))
  
  # the RGset may hold more samples (or a different order) than the Noob GRset
  common <- intersect(colnames(grset), colnames(rgset))
  if (length(common) == 0) {
    stop("No shared sample names. GRset: ", paste(head(colnames(grset), 3), collapse = ", "),
         " | RGset: ", paste(head(colnames(rgset), 3), collapse = ", "))
  }
  message("BloodExtended: ", length(common), " samples in common (GRset ",
          ncol(grset), ", RGset ", ncol(rgset), ")")
  
  detp  <- minfi::detectionP(rgset[, common])
  betas <- minfi::getBeta(grset)[, common, drop = FALSE]
  rm(rgset, grset); gc()
  
  # reference profiles (v1 library, Noob)
  ref  <- FlowSorted.Blood.EPIC::libraryDataGet("FlowSorted.BloodExtended.EPIC")
  refb <- minfi::getBeta(minfi::preprocessNoob(ref))
  ct   <- as.character(ref$CellType)
  rm(ref); gc()
  
  y <- collapse_lowest_detp(betas, detp,
                            intersect(IDOLOptimizedCpGsBloodExtended, rownames(refb)))
  rm(betas, detp); gc()
  
  # keep CpGs without missing values
  cpgs <- rownames(y)[rowSums(is.na(y)) == 0 & rowSums(is.na(refb[rownames(y), ])) == 0]
  y <- y[cpgs, , drop = FALSE]
  message("BloodExtended: using ", length(cpgs), " of ",
          length(IDOLOptimizedCpGsBloodExtended), " IDOL CpGs")
  
  types <- intersect(BLOOD_EXTENDED_CELLTYPES, ct)
  coef  <- sapply(types, function(k) rowMeans(refb[cpgs, ct == k, drop = FALSE]))
  prop  <- FlowSorted.Blood.EPIC::projectCellType_CP(y, coef)
  
  file <- study_file(config, "BloodExtended_cell_counts", chunk = chunk, ext = "csv")
  utils::write.csv(as.data.frame(prop), file, row.names = TRUE)
  invisible(file)
}

# ---------- EPIDISH CELL COMPOSITION ------------------------------------------
estimate_cell_counts_hepidish_noob <- function(config, chunk, method = "RPC") {
  require_package("EpiDISH")
  require_package("minfi")
  
  grset <- readRDS(noob_grset_file(config, chunk))
  betas <- get_betas(grset)      # Betas with plain cg IDs (replicate probes cgXXXX_TC11, _BC21, ... averaged)
  data("centEpiFibIC.m", package = "EpiDISH") 
  data("centBloodSub.m", package = "EpiDISH") 
  frac <- EpiDISH::hepidish(
    beta.m = betas,
    ref1.m = centEpiFibIC.m,
    ref2.m = centBloodSub.m, # 
    h.CT.idx = 3,
    method = method)
  
  file <- study_file(config, paste0("HEpiDISH_", method, "_cell_counts"), chunk = chunk, ext = "csv")
  utils::write.csv(as.data.frame(frac), file, row.names = TRUE)
  invisible(file)
}

# ---------- PC CLOCKS ---------------------------------------------------------
run_pc_clocks_noob <- function(config, chunk) {
  require_package("minfi")
  require_package("methylCIPHER")
  
  betas <- get_betas(readRDS(noob_grset_file(config, chunk)))
  manifest <- manifest_for_chunk(config, chunk)
  pheno <- manifest[, c("sample_id", "Age", "Female")]
  pheno$Age <- as.numeric(pheno$Age)
  pheno$Female <- as.numeric(pheno$Female)
  rownames(pheno) <- pheno$sample_id
  pheno <- pheno[colnames(betas), ]
  pc_clocks <- methylCIPHER::calcPCClocks(t(betas), pheno, RData = config$paths$PCClocks_data)
  
  file <- study_file(config, "PCClocks_noob", chunk = chunk, ext = "csv")
  utils::write.csv(as.data.frame(pc_clocks), file, row.names = TRUE)
  invisible(file)
}

# ---------- DunedinPACE -------------------------------------------------------
run_dunedinpace_noob <- function(config, chunk) {
  require_package("DunedinPACE")
  require_package("minfi")
  
  betas <- get_betas(readRDS(noob_grset_file(config, chunk)))
  pace <- DunedinPACE::PACEProjector(betas)
  
  file <- study_file(config, "DunedinPACE_noob", chunk = chunk, ext = "csv")
  utils::write.csv(as.data.frame(pace), file, row.names = TRUE)
  invisible(file)
}

# ---------- PLOTS -------------------------------------------------------
# Reads the CSV files written by the modules above (same chunk) and saves 3 PDFs
# in the pipeline figure folder. The code inside is plain step-by-step; the
# function wrapper is only there so the driver can run it once per chunk.
plot_noob_downstream <- function(config, chunk) {
  require_package("ggplot2")
  library(ggplot2)
  
  # ---- 1. Read the result files ----
  be <- read.csv(study_file(config, "BloodExtended_cell_counts", chunk = chunk, ext = "csv"), row.names = 1)
  he <- read.csv(study_file(config, "HEpiDISH_RPC_cell_counts",  chunk = chunk, ext = "csv"), row.names = 1)
  pc <- read.csv(study_file(config, "PCClocks_noob",             chunk = chunk, ext = "csv"), row.names = 1)
  
  # ---- 2. Boxplots of the cell counts, one panel per method ----
  long <- rbind(
    data.frame(method = "BloodExtended", cell = rep(colnames(be), each = nrow(be)), prop = unlist(be)),
    data.frame(method = "HEpiDISH",      cell = rep(colnames(he), each = nrow(he)), prop = unlist(he)))
  
  p1 <- ggplot(long, aes(cell, prop)) +
    geom_boxplot(outlier.size = 0.8) +
    facet_wrap(~method, scales = "free_x") +
    labs(x = NULL, y = "Estimated proportion") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # ---- 3. Agreement between the two methods on the shared cell types ----
  # BloodExtended subtypes are summed to match HEpiDISH:
  #   B = Bmem + Bnv,  CD4T = CD4mem + CD4nv + Treg,  CD8T = CD8mem + CD8nv
  # Basophils have no HEpiDISH counterpart, so they are left out.
  be_c <- data.frame(B = be$Bmem + be$Bnv, NK = be$NK,
                     CD4T = be$CD4mem + be$CD4nv + be$Treg, CD8T = be$CD8mem + be$CD8nv,
                     Mono = be$Mono, Neu = be$Neu, Eos = be$Eos, row.names = rownames(be))
  he_c <- data.frame(B = he$B, NK = he$NK, CD4T = he$CD4T, CD8T = he$CD8T,
                     Mono = he$Mono, Neu = he$Neutro, Eos = he$Eosino, row.names = rownames(he))
  
  same <- intersect(rownames(be_c), rownames(he_c))
  be_c <- be_c[same, ]
  he_c <- he_c[same, ]
  
  # rescale each sample so the shared groups sum to 1 (makes the methods comparable)
  be_c <- be_c / rowSums(be_c)
  he_c <- he_c / rowSums(he_c)
  
  # correlation per cell type, shown in each panel title
  cells  <- colnames(be_c)
  r      <- sapply(cells, function(k) suppressWarnings(cor(be_c[[k]], he_c[[k]], use = "complete.obs")))
  labels <- paste0(cells, " (r = ", round(r, 2), ")")
  
  both <- data.frame(cell = factor(rep(labels, each = nrow(be_c)), levels = labels),
                     BloodExtended = unlist(be_c), HEpiDISH = unlist(he_c))
  
  p2 <- ggplot(both, aes(BloodExtended, HEpiDISH)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "firebrick", linewidth = 0.6) +
    facet_wrap(~cell, scales = "free") +
    labs(x = "BloodExtended (relative proportion)", y = "HEpiDISH (relative proportion)") +
    theme_bw()
  
  # ---- 4. Clocks vs chronological age ----
  clocks <- intersect(c("PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge", "PCGrimAge"), colnames(pc))
  
  clk <- NULL
  for (k in clocks) {
    ok  <- !is.na(pc$Age) & !is.na(pc[[k]])
    r   <- cor(pc$Age[ok], pc[[k]][ok])
    err <- median(abs(pc[[k]][ok] - pc$Age[ok]))
    clk <- rbind(clk, data.frame(clock = sprintf("%s\nr = %.2f, median error = %.1f y", k, r, err),
                                 Age = pc$Age[ok], estimate = pc[[k]][ok]))
  }
  clk$clock <- factor(clk$clock, levels = unique(clk$clock))
  
  p3 <- ggplot(clk, aes(Age, estimate)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "firebrick", linewidth = 0.6) +
    facet_wrap(~clock, scales = "free_y") +
    labs(x = "Chronological age (years)", y = "Clock estimate (years)") +
    theme_bw()
  
  # ---- 5. Save the three figures ----
  fig_dir <- config$paths$figure_dir %||% output_figure_dir(config)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig <- function(name) study_file(config, name, chunk = chunk, ext = "pdf", dir = fig_dir)
  ggsave(fig("cell_counts_boxplots"),           p1, width = 12, height = 5)
  suppressWarnings(ggsave(fig("cell_counts_method_correlation"), p2, width = 12, height = 8))
  ggsave(fig("clocks_vs_age"),                  p3, width = 12, height = 7)
  
  invisible(NULL)
}

# ---------- RUN NOOB DOWNSTREAM -----------------------------------------------
run_noob_downstream <- function(config, chunk) {
  opts <- config$noob_downstream %||% list()
  
  # a failing module warns instead of stopping the others
  safe <- function(name, expr) {
    tryCatch(expr, error = function(e) {
      warning("Chunk ", chunk, ": ", name, " failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
  }
  
  # Existing config keys keep working: `cell_counts: true` runs both cell-count
  # methods, and plots are made when cell counts and PC clocks are both on.
  # Optional keys (blood_extended, hepidish, plots) override this if you add them.
  run_cells  <- isTRUE(opts$cell_counts)
  run_be     <- isTRUE(opts$blood_extended %||% run_cells)
  run_he     <- isTRUE(opts$hepidish       %||% run_cells)
  run_clocks <- isTRUE(opts$pc_clocks)
  run_plots  <- isTRUE(opts$plots %||% (run_cells && run_clocks))
  
  if (run_be)                   safe("BloodExtended", estimate_cell_counts_blood_extended(config, chunk))
  if (run_he)                   safe("HEpiDISH",      estimate_cell_counts_hepidish_noob(config, chunk))
  if (run_clocks)               safe("PC clocks",     run_pc_clocks_noob(config, chunk))
  if (isTRUE(opts$dunedinpace)) safe("DunedinPACE",   run_dunedinpace_noob(config, chunk))
  if (run_plots)                safe("Plots",         plot_noob_downstream(config, chunk))
  invisible(NULL)
}