# ------------------------------------------------------------------------------
# BMIQ downstream. All chunks by default; --chunk=N runs one chunk, and its merged
# table and figures then cover that chunk only (file names end in _chunkN).
#   0. M values from the BMIQ betas          (per chunk, saved per chunk)
#   1. FraminghamPACE on the BMIQ betas      (per chunk, saved per chunk)
#   2. merge of all downstream outputs into ONE table (cell_composition_and_clocks, in Derived_Variables)
#   3. figures: cell counts, PC clocks vs age, pace-score distributions and by sex,
#      beta and M value distributions
# ------------------------------------------------------------------------------

# ---------- 0. Settings -------------------------------------------------------
DOWNSTREAM_OUTPUTS <- c("Blood6_cell_counts", "BloodExtended_cell_counts",
                        "PCClocks_noob", "SystemsAge_noob", "DunedinPACE_noob",
                        "FraminghamPACE_BMIQ")

# Prefix of each output's columns in the joined tables, so names never collide when
# tools are added. A column that already starts with the prefix is left as it is
# (PCHorvath1, DunedinPACE, FraminghamPACE_Ridge, SystemsAge).
# For a new tool: add its output name to DOWNSTREAM_OUTPUTS and a line here.
# Without a line, the prefix is the output name without _noob, _BMIQ or _cell_counts.
OUTPUT_PREFIX <- c(
  Blood6_cell_counts        = "Blood6",
  BloodExtended_cell_counts = "BloodExtended",
  PCClocks_noob             = "PC",
  SystemsAge_noob           = "SystemsAge",
  DunedinPACE_noob          = "DunedinPACE",
  FraminghamPACE_BMIQ       = "FraminghamPACE"
)

# Cell types grouped by lineage for the cell-count figure. The order of this vector
# is the order on the x axis. Cell types that are not listed are shown last, in grey.
CELL_LINEAGE <- c(
  Neu = "Granulocytes", Eos = "Granulocytes", Bas = "Granulocytes",
  Mono = "Monocytes",
  NK = "NK cells",
  Bnv = "B cells", Bmem = "B cells", Bcell = "B cells",
  CD4nv = "CD4 T cells", CD4mem = "CD4 T cells", CD4T = "CD4 T cells", Treg = "CD4 T cells",
  CD8nv = "CD8 T cells", CD8mem = "CD8 T cells", CD8T = "CD8 T cells"
)
# colour-blind-safe palette (Okabe-Ito)
LINEAGE_COLOURS <- c(
  "Granulocytes" = "#E69F00", "Monocytes" = "#56B4E9", "NK cells" = "#009E73",
  "B cells" = "#CC79A7", "CD4 T cells" = "#0072B2", "CD8 T cells" = "#D55E00",
  "Other" = "grey60"
)

# colour-blind-safe colours (Okabe-Ito) for chunks and other groups
PUB_COLOURS <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#999999")

# fill of the pace-score histograms (muted slate blue)
PACE_FILL <- "#5C7FA3"

# One theme for all figures: classic axes, light grid, black text, legend at the bottom
theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text         = ggplot2::element_text(colour = "black"),
      panel.grid.major  = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      plot.title        = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle     = ggplot2::element_text(size = base_size - 2, colour = "grey30"),
      strip.background  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(size = base_size - 1),
      legend.position   = "bottom",
      legend.text       = ggplot2::element_text(size = base_size - 2),
      legend.key.size   = grid::unit(0.4, "cm")
    )
}

# Runs one step with start/finish banners; an error becomes a warning
run_safely <- function(name, expr) {
  message("\n======= Running ", name, " =======")
  tryCatch({
    out <- expr
    message("======= ", name, " completed correctly =======\n")
    out
  }, error = function(e) {
    message("======= ", name, " FAILED: ", conditionMessage(e), " =======\n")
    warning(name, " failed: ", conditionMessage(e), call. = FALSE)
    NULL
  })
}

# Pearson correlation and p-value of two vectors (complete pairs only);
# NA if there are fewer than 3 pairs
cor_summary <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 3) return(c(r = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- cor.test(x[ok], y[ok])
  c(r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

# "<study_name>_": the prefix of every output file, with the same cleaning of the
# study name as study_file() (safe_study_name() in paths.R)
study_prefix <- function(config) {
  paste0(safe_study_name(config), "_")
}

# Path of a figure: <figures>/<study_name>_<name><suffix>.pdf
figure_file <- function(config, name, suffix = "") {
  file.path(config$paths$output_figures, paste0(study_prefix(config), name, suffix, ".pdf"))
}

# Is this plot switched on? `plots: false` switches all off; `plots: true` or a
# missing key means on; a list switches single plots on or off.
plot_wanted <- function(config, key) {
  plot_opt <- config$bmiq_downstream$plots %||% TRUE
  if (is.list(plot_opt)) isTRUE(plot_opt[[key]] %||% TRUE) else isTRUE(plot_opt)
}

# "p < 1e-4" or "p = 0.032" (for panel titles)
p_text <- function(p) {
  p <- unname(p)
  if (is.na(p)) "p = NA" else if (p < 1e-4) "p < 1e-4" else paste0("p = ", format(signif(p, 2), scientific = FALSE))
}

# Joins module tables by sample_id (one row per sample). Age and Female are kept
# once; every other column gets the tool prefix from OUTPUT_PREFIX.
# Returns NULL for an empty list.
join_modules <- function(tables) {
  joined    <- NULL
  have_demo <- FALSE
  for (name in names(tables)) {
    d    <- tables[[name]]
    demo <- intersect(c("Age", "Female"), names(d))
    if (length(demo)) {
      if (have_demo) d <- d[, setdiff(names(d), demo), drop = FALSE] else have_demo <- TRUE
    }
    prefix <- if (name %in% names(OUTPUT_PREFIX)) OUTPUT_PREFIX[[name]] else
      sub("_(noob|BMIQ|cell_counts)$", "", name)
    cols <- !names(d) %in% c("sample_id", "Age", "Female")
    old_names <- names(d)[cols]
    names(d)[cols] <- ifelse(startsWith(old_names, prefix), old_names, paste0(prefix, "_", old_names))
    joined <- if (is.null(joined)) d else merge(joined, d, by = "sample_id", all = TRUE)
  }
  joined
}

# ---------- 0. M VALUES (BMIQ betas) -------------------------------------------
# M = log2(beta / (1 - beta)). Same rows and columns as the BMIQ betas. An infinite M
# (a beta of exactly 0 or 1) is set to NA and counted. The log reports how long the
# tails are: values near 0 or 1 give large |M|.
calculate_m_values <- function(config, chunk) {
  betas <- readRDS(bmiq_betas_file(config, chunk))
  if (any(betas < 0 | betas > 1, na.rm = TRUE)) {
    stop("Values outside 0-1: these do not look like beta values")
  }
  
  m <- log2(betas / (1 - betas))
  dimnames(m) <- dimnames(betas)
  n_inf <- sum(is.infinite(m))
  m[is.infinite(m)] <- NA
  
  q <- stats::quantile(m, c(0.001, 0.999), na.rm = TRUE, names = FALSE)
  message("M values (chunk ", chunk, "): ", format(nrow(m), big.mark = ","), " CpGs x ", ncol(m),
          " samples; range ", round(min(m, na.rm = TRUE), 1), " to ", round(max(m, na.rm = TRUE), 1))
  message("  tails: 0.1% / 99.9% quantiles ", round(q[1], 1), " / ", round(q[2], 1),
          "; |M| > 5: ", round(100 * mean(abs(m) > 5, na.rm = TRUE), 2), "%",
          "; |M| > 10: ", round(100 * mean(abs(m) > 10, na.rm = TRUE), 3), "%",
          "; infinite (beta of exactly 0 or 1) set to NA: ", format(n_inf, big.mark = ","))
  
  file <- study_file(config, "BMIQ_Mvalues", chunk = chunk, ext = "rds")
  saveRDS(m, file)
  invisible(file)
}

# ---------- 1. FRAMINGHAMPACE (BMIQ betas) ------------------------------------
# Models: ElasticNet and Ridge (package FraminghamPACE).
# Probe coverage is printed as messages.
run_framinghampace_bmiq <- function(config, chunk) {
  require_package("FraminghamPACE")
  
  models <- c("ElasticNet", "Ridge")
  
  # betas: rows = CpGs, columns = samples
  betas <- readRDS(bmiq_betas_file(config, chunk))
  if (is.null(rownames(betas)) || is.null(colnames(betas))) {
    stop("BMIQ betas have no row or column names; check that bmiq_chunk() keeps the dimnames")
  }
  
  # EPICv2: average replicate probes (cgXXXX_TC11, _BC21, ...) to plain cg IDs
  ids <- sub("_.*$", "", rownames(betas))
  if (anyDuplicated(ids)) {
    ok <- !is.na(betas)
    betas[!ok] <- 0
    betas <- rowsum(betas, ids, reorder = FALSE) / rowsum(ok * 1, ids, reorder = FALSE)
    betas[is.nan(betas)] <- NA
    message("Replicate probes averaged to plain cg IDs")
  }
  
  stopifnot(
    "Column names (sample IDs) must be unique"  = !anyDuplicated(colnames(betas)),
    "Row names don't look like CpG IDs (cg...)" = mean(grepl("^cg", rownames(betas))) > 0.5
  )
  if (any(betas < 0 | betas > 1, na.rm = TRUE)) {
    stop("Values outside 0-1: these do not look like beta values")
  }
  
  message("Input: ", nrow(betas), " CpGs x ", ncol(betas), " samples, ",
          round(100 * mean(is.na(betas)), 2), "% missing beta values")
  
  probe_types <- table(substr(rownames(betas), 1, 2))
  message("Probe types in input (the package uses only cg probes): ",
          paste(names(probe_types), probe_types, sep = " = ", collapse = ", "))
  
  # ---- probe coverage before running (numbers for the methods section) ----
  bg <- FraminghamPACE::getFraminghamPACEProbes(backgroundList = TRUE)
  bg_n <- c(present = sum(bg %in% rownames(betas)), total = length(bg))
  message(sprintf("Background probes present: %d / %d", bg_n[["present"]], bg_n[["total"]]))
  
  scoring_n <- list()
  for (m in models) {
    sp <- FraminghamPACE::getFraminghamPACEProbes(model_type = m)
    scoring_n[[m]] <- c(present = sum(sp %in% rownames(betas)), total = length(sp))
    message(sprintf("%s scoring probes present: %d / %d", m, scoring_n[[m]][["present"]], scoring_n[[m]][["total"]]))
  }
  coverage_rows <- list()   # saved as a table, for the report
  
  # ---- run each model ----
  results <- data.frame(sample_id = colnames(betas), stringsAsFactors = FALSE)
  
  for (m in models) {
    message("--- FraminghamPACE ", m, " ---")
    out <- FraminghamPACE::getFraminghamPACE(betas = betas, model_type = m, return_coverage = TRUE)
    
    stopifnot(all(results$sample_id %in% rownames(out$pace_score)))
    results[[paste0("FraminghamPACE_", m)]] <- out$pace_score[results$sample_id, 1]
    
    # one coverage value per model (a range is shown only if samples differ)
    cov <- round(range(out$coverage$sample, na.rm = TRUE), 1)
    cov_txt <- if (cov[1] == cov[2]) sprintf("%.1f", cov[1]) else sprintf("%.1f-%.1f", cov[1], cov[2])
    n_scored <- sum(!is.na(results[[paste0("FraminghamPACE_", m)]]))
    message(sprintf("%s: probe coverage %s%% | samples scored: %d / %d", m, cov_txt, n_scored, nrow(results)))
    coverage_rows[[m]] <- data.frame(
      chunk = chunk, model = m,
      background_probes_present = bg_n[["present"]], background_probes_total = bg_n[["total"]],
      scoring_probes_present = scoring_n[[m]][["present"]], scoring_probes_total = scoring_n[[m]][["total"]],
      predictive_coverage_min = cov[1], predictive_coverage_max = cov[2],
      samples_scored = n_scored, samples_total = nrow(results), stringsAsFactors = FALSE)
  }
  utils::write.csv(do.call(rbind, coverage_rows),
                   study_file(config, "FraminghamPACE_coverage", chunk = chunk, ext = "csv"), row.names = FALSE)
  
  file <- study_file(config, "FraminghamPACE_BMIQ", chunk = chunk, ext = "csv")
  utils::write.csv(results, file, row.names = FALSE)
  invisible(file)
}

# ---------- 2. MERGE OUTPUTS ---------------------------------------------------
merge_downstream_outputs <- function(config, chunks, suffix = "") {
  tables  <- list()
  
  # ---- A. stack the chunks of each module ----
  for (name in DOWNSTREAM_OUTPUTS) {
    files <- vapply(chunks,
                    function(k) study_file(config, name, chunk = k, ext = "csv"),
                    character(1))
    found <- file.exists(files)
    if (!any(found)) next
    if (!all(found)) {
      warning(name, ": ", sum(!found), " of ", length(chunks),
              " chunk file(s) missing; merging the rest", call. = FALSE)
    }
    
    d <- do.call(rbind, lapply(files[found], read.csv,
                               check.names = FALSE, stringsAsFactors = FALSE))
    stopifnot(!anyDuplicated(d$sample_id))
    tables[[name]] <- d
    
    message(name, ": ", nrow(d), " samples from ", sum(found), " chunk(s)")
  }
  if (!length(tables)) stop("No downstream output files were found to merge")
  
  # ---- B. all modules in one table ----
  all_tbl  <- join_modules(tables)
  all_file <- derived_file(config, "cell_composition_and_clocks", suffix)
  utils::write.csv(add_name_column(config, all_tbl), all_file, row.names = FALSE)
  message("All modules merged (", paste(names(tables), collapse = ", "), "): ",
          nrow(all_tbl), " samples x ", ncol(all_tbl) - 1, " variables -> ", basename(all_file))
  
  invisible(list(tables = tables, all = all_tbl))
}

# ---------- 3. PLOTS -----------------------------------------------------------
# Reads the merged tables and saves PDFs (vector, at most 7.2 in wide) in the figures
# folder; file names get the study name and, for a single chunk, the chunk suffix.
# All figures share theme_pub(). Each plot can be switched off in the config (see the
# top of this file). `plots: false` switches all plots off; `plots: true` or a
# missing key means on. A plot is also skipped when its module produced no output.
# Correlations are Pearson; the p-value is from cor.test().
plot_downstream_results <- function(config, res, suffix = "") {
  require_package("ggplot2")
  library(ggplot2)
  
  tables  <- res$tables
  all_tbl <- res$all
  want    <- function(key) plot_wanted(config, key)
  figs    <- list()   # the plots, with their size; the report uses them too
  
  fig_dir <- config$paths$output_figures
  id_cols <- c("sample_id", "Age", "Female")
  
  # ---- cell counts: boxplots (one figure per reference) ----
  # Boxes by lineage, every sample as a point, proportions in %.
  if (want("cell_counts")) {
    for (name in grep("_cell_counts$", names(tables), value = TRUE)) {
      cc    <- tables[[name]][, setdiff(names(tables[[name]]), id_cols), drop = FALSE]
      label <- sub("_cell_counts$", "", name)
      
      cells <- c(intersect(names(CELL_LINEAGE), colnames(cc)), setdiff(colnames(cc), names(CELL_LINEAGE)))
      long  <- data.frame(cell = factor(rep(cells, each = nrow(cc)), levels = cells),
                          prop = unlist(cc[cells], use.names = FALSE))
      long$lineage <- unname(CELL_LINEAGE[as.character(long$cell)])
      long$lineage[is.na(long$lineage)] <- "Other"
      long$lineage <- factor(long$lineage, levels = c(unique(CELL_LINEAGE), "Other"))
      
      rotate <- length(cells) > 12
      p1 <- ggplot(long, aes(cell, prop, fill = lineage)) +
        geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4,
                     colour = "grey20", alpha = 0.7) +
        geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1),
                   size = 0.8, alpha = 0.6, colour = "grey20", show.legend = FALSE) +
        scale_fill_manual(values = LINEAGE_COLOURS, name = NULL) +
        scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                           expand = expansion(mult = c(0.02, 0.05))) +
        guides(fill = guide_legend(nrow = 1)) +
        labs(title = sprintf("Blood cell composition (%s reference)", label),
             subtitle = sprintf("n = %d samples", nrow(cc)),
             x = NULL, y = "Estimated proportion of cells") +
        theme_pub() +
        theme(axis.text.x = element_text(angle = if (rotate) 45 else 0, hjust = if (rotate) 1 else 0.5))
      fig_w <- max(6.5, min(7.2, 0.55 * length(cells) + 1)); fig_h <- 4.2
      ggsave(figure_file(config, paste0("cell_counts_boxplots_", label), suffix), p1,
             width = fig_w, height = fig_h)
      figs[[paste0("cell_counts_", label)]] <- list(plot = p1, width = fig_w, height = fig_h)
    }
  }
  
  # ---- PC clocks vs chronological age ----
  pc <- tables[["PCClocks_noob"]]
  if (want("pc_clocks") && !is.null(pc)) {
    clocks <- intersect(c("PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge", "PCGrimAge"), colnames(pc))
    
    clk <- NULL
    for (k in clocks) {
      cs  <- cor_summary(pc$Age, pc[[k]])
      ok  <- !is.na(pc$Age) & !is.na(pc[[k]])
      err <- median(abs(pc[[k]][ok] - pc$Age[ok]))
      clk <- rbind(clk, data.frame(
        clock = sprintf("%s\nr = %.2f, %s\nmedian error = %.1f y", k, cs["r"], p_text(cs["p"]), err),
        Age = pc$Age[ok], estimate = pc[[k]][ok]))
    }
    
    if (!is.null(clk)) {
      clk$clock <- factor(clk$clock, levels = unique(clk$clock))
      p2 <- ggplot(clk, aes(Age, estimate)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
        geom_point(alpha = 0.7, size = 1.3, colour = "#0072B2") +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#D55E00", linewidth = 0.7) +
        facet_wrap(~clock, scales = "free_y", ncol = 3) +
        labs(title = "PC clocks vs chronological age",
             subtitle = sprintf("n = %d samples; dashed line: identity, solid line: linear fit", nrow(pc)),
             x = "Chronological age (years)", y = "Epigenetic age estimate (years)") +
        theme_pub()
      fig_h <- 2.6 * ceiling(length(clocks) / 3) + 1
      ggsave(figure_file(config, "clocks_vs_age", suffix), p2, width = 7.2, height = fig_h)
      figs[["clocks_vs_age"]] <- list(plot = p2, width = 7.2, height = fig_h)
    }
  }
  
  # ---- SystemsAge: the overall score and the system scores vs chronological age ----
  # The scores are on the scale of the SystemsAge model (not necessarily years), so
  # there is no identity line; r and p show how each score relates to age.
  sa <- tables[["SystemsAge_noob"]]
  if (want("systems_age") && !is.null(sa)) {
    sa_order   <- c("SystemsAge", "Heart", "Lung", "Kidney", "Liver", "Brain", "Immune",
                    "Inflammation", "Blood", "MusculoSkeletal", "Hormone", "Metabolic")
    score_cols <- setdiff(names(sa), id_cols)
    score_cols <- c(intersect(sa_order, score_cols), setdiff(score_cols, sa_order))
    
    sa_long <- NULL
    for (k in score_cols) {
      cs <- cor_summary(sa$Age, sa[[k]])
      ok <- !is.na(sa$Age) & !is.na(sa[[k]])
      sa_long <- rbind(sa_long, data.frame(
        panel = sprintf("%s\nr = %.2f, %s", k, cs["r"], p_text(cs["p"])),
        Age = sa$Age[ok], score = sa[[k]][ok]))
    }
    
    if (!is.null(sa_long)) {
      sa_long$panel <- factor(sa_long$panel, levels = unique(sa_long$panel))
      p6 <- ggplot(sa_long, aes(Age, score)) +
        geom_point(alpha = 0.7, size = 1.1, colour = "#0072B2") +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#D55E00", linewidth = 0.6) +
        facet_wrap(~panel, scales = "free_y", ncol = 4) +
        labs(title = "SystemsAge scores vs chronological age",
             subtitle = sprintf("n = %d samples; solid line: linear fit", nrow(sa)),
             x = "Chronological age (years)", y = "Score") +
        theme_pub()
      fig_h <- 1.9 * ceiling(length(score_cols) / 4) + 1
      ggsave(figure_file(config, "systemsage_vs_age", suffix), p6, width = 7.2, height = fig_h)
      figs[["systemsage_vs_age"]] <- list(plot = p6, width = 7.2, height = fig_h)
    }
  }
  
  # ---- distribution of the pace scores: histograms, and by sex ----
  # (DunedinPACE and FraminghamPACE are rates, not ages: no plot against age)
  fr_cols   <- grep("^FraminghamPACE_", names(all_tbl), value = TRUE)
  pace_cols <- intersect(c("DunedinPACE", fr_cols), names(all_tbl))
  if (length(pace_cols) && (want("pace_distribution") || want("pace_by_sex"))) {
    pace_names <- sub("^FraminghamPACE_", "FraminghamPACE ", pace_cols)
    has_female <- "Female" %in% names(all_tbl)
    pace_long <- data.frame(
      score_name = factor(rep(pace_names, each = nrow(all_tbl)), levels = pace_names),
      score      = unlist(all_tbl[pace_cols], use.names = FALSE),
      Female     = if (has_female) rep(all_tbl$Female, times = length(pace_cols)) else NA_real_
    )
    pace_long  <- pace_long[!is.na(pace_long$score), ]
    pace_width <- max(4, min(7.2, 2.6 * length(pace_cols) + 0.5))
    
    if (want("pace_distribution")) {
      p4 <- ggplot(pace_long, aes(score)) +
        geom_histogram(bins = 15, fill = PACE_FILL, colour = "grey20", linewidth = 0.3)
      if ("DunedinPACE" %in% pace_cols) {   # 1 = aging at the calendar rate
        ref <- data.frame(score_name = factor("DunedinPACE", levels = pace_names), x = 1)
        p4  <- p4 + geom_vline(data = ref, aes(xintercept = x), linetype = "dashed", colour = "grey30")
      }
      p4 <- p4 +
        facet_wrap(~score_name, scales = "free") +
        labs(title = "Distribution of the pace-of-aging scores",
             subtitle = sprintf("n = %d samples; dashed line: DunedinPACE = 1 (aging at the calendar rate)",
                                nrow(all_tbl)),
             x = "Score", y = "Samples") +
        theme_pub()
      ggsave(figure_file(config, "pace_scores_distribution", suffix), p4,
             width = pace_width, height = 3.6)
      figs[["pace_distribution"]] <- list(plot = p4, width = pace_width, height = 3.6)
    }
    
    if (want("pace_by_sex")) {
      d <- pace_long[!is.na(pace_long$Female), ]
      if (nrow(d) == 0) {
        message("Pace scores by sex skipped: no Female column (run PC clocks or SystemsAge, which carry it)")
      } else {
        d$Sex <- factor(d$Female, levels = c(0, 1), labels = c("Male", "Female"))
        if (nlevels(droplevels(d$Sex)) < 2) {
          message("Pace scores by sex skipped: only one sex in the data")
        } else {
          # Wilcoxon rank-sum p-value in each panel title
          panel <- vapply(levels(d$score_name), function(s) {
            x <- d[d$score_name == s, ]
            p <- tryCatch(suppressWarnings(wilcox.test(score ~ Sex, data = x)$p.value),
                          error = function(e) NA_real_)
            sprintf("%s\nWilcoxon %s", s, p_text(p))
          }, character(1))
          d$panel <- factor(as.character(d$score_name), levels = levels(d$score_name), labels = panel)
          n_sex <- table(d$Sex[d$score_name == levels(d$score_name)[1]])
          
          p5 <- ggplot(d, aes(Sex, score, fill = Sex)) +
            geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.4, colour = "grey20", alpha = 0.7) +
            geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1),
                       size = 0.9, alpha = 0.6, colour = "grey20") +
            scale_fill_manual(values = c(Male = "#0072B2", Female = "#D55E00"), guide = "none") +
            facet_wrap(~panel, scales = "free_y") +
            labs(title = "Pace-of-aging scores by sex",
                 subtitle = sprintf("Male n = %d, Female n = %d", n_sex[["Male"]], n_sex[["Female"]]),
                 x = NULL, y = "Score") +
            theme_pub()
          ggsave(figure_file(config, "pace_scores_by_sex", suffix), p5,
                 width = pace_width, height = 3.8)
          figs[["pace_by_sex"]] <- list(plot = p5, width = pace_width, height = 3.8)
        }
      }
    }
  }
  
  message("Figures saved in: ", fig_dir)
  invisible(figs)
}

# ---------- 3b. BETA AND M VALUE DISTRIBUTIONS ---------------------------------
# One density curve per sample, read one chunk at a time (nothing large is kept).
# kind: "beta" (BMIQ betas, 0 to 1) or "M" (M values, drawn from -20 to 20).
# Returns a long table (chunk, sample_id, x, density), or NULL if there is no data;
# attr(, "value_range") is the observed minimum and maximum of the values.
sample_densities <- function(config, chunks, kind = c("beta", "M")) {
  kind <- match.arg(kind)
  rng  <- if (kind == "beta") c(0, 1) else c(-20, 20)
  out  <- list()
  vmin <- Inf   # observed range of the values, all samples
  vmax <- -Inf
  
  for (chunk in chunks) {
    file <- if (kind == "beta") bmiq_betas_file(config, chunk) else
      study_file(config, "BMIQ_Mvalues", chunk = chunk, ext = "rds")
    if (!file.exists(file)) {
      message("  ", kind, " values of chunk ", chunk, " not found; skipped")
      next
    }
    mat  <- readRDS(file)
    vmin <- min(vmin, min(mat, na.rm = TRUE))
    vmax <- max(vmax, max(mat, na.rm = TRUE))
    for (j in seq_len(ncol(mat))) {
      x <- mat[, j]
      x <- x[!is.na(x)]
      if (length(x) < 2) next   # e.g. a sample where BMIQ failed
      dn <- stats::density(x, from = rng[1], to = rng[2], n = 256)
      out[[length(out) + 1]] <- data.frame(chunk = chunk, sample_id = colnames(mat)[j],
                                           x = dn$x, density = dn$y)
    }
    rm(mat); gc()
  }
  if (!length(out)) return(NULL)
  if (vmin < rng[1] || vmax > rng[2]) {
    message("  ", kind, " values outside ", rng[1], " to ", rng[2], " are not drawn")
  }
  res <- do.call(rbind, out)
  attr(res, "value_range") <- c(vmin, vmax)
  res
}

# All samples of the chunks that were run: thin lines are samples (coloured by
# chunk when there are several), the thick black line is the mean of all samples.
plot_value_distributions <- function(config, chunks, suffix = "") {
  require_package("ggplot2")
  library(ggplot2)
  
  kinds <- list(
    beta = list(key = "beta_distribution", file = "bmiq_beta_distribution",
                xlab = "Beta value (BMIQ)", title = "Distribution of BMIQ beta values"),
    M    = list(key = "m_distribution", file = "bmiq_mvalue_distribution",
                xlab = "M value (from BMIQ betas)", title = "Distribution of M values")
  )
  
  figs <- list()   # the plots, with their size; the report uses them too
  for (kind in names(kinds)) {
    k <- kinds[[kind]]
    if (!plot_wanted(config, k$key)) next
    
    d <- sample_densities(config, chunks, kind)
    if (is.null(d)) {
      message(k$title, " skipped: no data found")
      next
    }
    vr <- attr(d, "value_range")   # observed minimum and maximum of the values
    d$chunk <- factor(d$chunk)
    mean_curve <- aggregate(density ~ x, data = d, FUN = mean)
    n_samples  <- length(unique(d$sample_id))
    
    # M: the x axis ends at the observed range (plus half a unit). The curves are smoothed,
    # so they would otherwise run on past the last real value.
    x_lim <- if (kind == "M") c(floor(vr[1]) - 0.5, ceiling(vr[2]) + 0.5) else NULL
    digits <- if (kind == "beta") 3 else 1
    
    p <- ggplot(d, aes(x, density, group = sample_id)) +
      geom_line(aes(colour = chunk), alpha = 0.6, linewidth = 0.35) +
      geom_line(data = mean_curve, aes(x, density), inherit.aes = FALSE,
                colour = "black", linewidth = 0.9) +
      scale_colour_manual(values = rep_len(PUB_COLOURS, nlevels(d$chunk)), name = "Chunk") +
      guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, linewidth = 1))) +
      labs(title = k$title,
           subtitle = sprintf("%d samples, %d chunk(s); observed range %s to %s; black line: mean of all samples",
                              n_samples, nlevels(d$chunk),
                              format(round(vr[1], digits)), format(round(vr[2], digits))),
           x = k$xlab, y = "Density") +
      theme_pub()
    # a legend only helps with a few chunks
    if (nlevels(d$chunk) == 1 || nlevels(d$chunk) > 8) p <- p + theme(legend.position = "none")
    if (!is.null(x_lim)) p <- p + coord_cartesian(xlim = x_lim)
    ggsave(figure_file(config, k$file, suffix), p,
           width = 7.2, height = 4.2)
    figs[[kind]] <- list(plot = p, width = 7.2, height = 4.2)
  }
  
  invisible(figs)
}

# ---------- DRIVERS ------------------------------------------------------------
# After the per-chunk steps: merge outputs, then plots
run_downstream_summary <- function(config, chunks, suffix = "") {
  opts <- config$bmiq_downstream %||% list()
  res  <- NULL
  
  if (isTRUE(opts$merge %||% TRUE)) {
    res <- run_safely("Merge of downstream outputs", merge_downstream_outputs(config, chunks, suffix))
  }
  figs <- list()   # the plots are also used in the report
  if (!isFALSE(opts$plots)) {
    if (!is.null(res)) {
      f <- run_safely("Summary plots", plot_downstream_results(config, res, suffix))
      if (is.list(f)) figs <- c(figs, f)
    }
    # the distributions read the per-chunk BMIQ and M value files, not the merged table
    f <- run_safely("Beta and M value distribution plots", plot_value_distributions(config, chunks, suffix))
    if (is.list(f)) figs <- c(figs, f)
  }
  # the sample QC tables of all chunks, in Derived_Variables (R/derived_variables.R)
  run_safely("Sample QC tables (all chunks)", export_qc_tables(config, chunks, suffix))
  # the control-probe PCA scores, in Derived_Variables (R/derived_variables.R)
  run_safely("Control-probe PCA scores", export_control_pca_scores(config, chunks, suffix))
  # the BMIQ beta and M value matrices, in Derived_Variables; switch off with bmiq_downstream: derived_matrices: false
  if (isTRUE(opts$derived_matrices %||% TRUE)) {
    run_safely("BMIQ beta and M value matrices (all chunks)", export_bmiq_matrices(config, chunks, suffix))
  }
  # the report (R/report.R): a text file with the numbers and a PDF; switch off with bmiq_downstream: report: false
  if (isTRUE(opts$report %||% TRUE)) {
    nums <- run_safely("Report numbers", collect_report_numbers(config, chunks))
    if (!is.null(nums)) {
      run_safely("Methods summary (text file)", write_methods_summary(config, nums, chunks, suffix))
      run_safely("PDF report", write_report_pdf(config, nums, chunks, suffix, figs))
    }
  }
  invisible(res)
}

# Entry point: all chunks, or the one chunk selected with --chunk=N
run_bmiq_downstream <- function(config) {
  opts     <- config$bmiq_downstream %||% list()
  chunks   <- get_chunks_to_run(config)
  n_chunks <- as.integer(config$chunking$n_chunks %||% 1)
  
  # one chunk selected out of several: its outputs and plots carry a _chunkN suffix
  suffix <- if (length(chunks) < n_chunks) paste0("_chunk", paste(chunks, collapse = "-")) else ""
  if (nzchar(suffix)) {
    message("Only chunk ", paste(chunks, collapse = ", "), " of ", n_chunks,
            ": merged tables and plots will cover this chunk only")
  }
  
  # 0 and 1: per chunk, one chunk in memory at a time; outputs are saved per chunk
  for (chunk in chunks) {
    message("BMIQ downstream for chunk ", chunk)
    run_safely(paste0("M values (chunk ", chunk, ")"), calculate_m_values(config, chunk))
    if (isTRUE(opts$framinghampace)) {
      run_safely(paste0("FraminghamPACE (chunk ", chunk, ")"), run_framinghampace_bmiq(config, chunk))
    }
    gc()
  }
  
  # 2 and 3: the chunks that were run
  run_downstream_summary(config, chunks, suffix)
  invisible(NULL)
}