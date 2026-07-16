read_probe_vector <- function(file, column = NULL, header = TRUE) {
  if (!file.exists(file)) {
    return(character())
  }
  if (grepl("[.]txt$", file, ignore.case = TRUE) && is.null(column)) {
    return(unique(trimws(readLines(file, warn = FALSE))))
  }
  x <- utils::read.csv(file, stringsAsFactors = FALSE, header = header, check.names = FALSE)
  if (!nrow(x)) {
    return(character())
  }
  if (!is.null(column) && column %in% colnames(x)) {
    return(unique(as.character(x[[column]])))
  }
  unique(as.character(x[[1]]))
}

cross_reactive_probe_ids <- function(config) {
  opts <- config$probe_filtering %||% list()
  array_type <- tolower(opts$array_type %||% "EPICv1")
  files <- opts$cross_reactive_files %||% list()
  probe_dir <- config_path(config, config$paths$probe_list_dir %||% "data-raw/probe_lists")

  if (array_type %in% c("epic", "epicv1", "850k", "epic_v1")) {
    unique(c(
      read_probe_vector(file.path(probe_dir, files$epic_pidsley %||% "13059_2016_1066_MOESM1_ESM.csv")),
      read_probe_vector(file.path(probe_dir, files$epic_mccartney_cpg %||% "1-s2.0-S221359601630071X-mmc2.txt"), header = FALSE),
      read_probe_vector(file.path(probe_dir, files$epic_mccartney_noncpg %||% "1-s2.0-S221359601630071X-mmc3.txt"), header = FALSE)
    ))
  } else if (array_type %in% c("450k", "k450", "hm450")) {
    unique(c(
      read_probe_vector(file.path(probe_dir, files$k450_chen %||% "48639-non-specific-probes-Illumina450k.csv"), column = "TargetID"),
      read_probe_vector(file.path(probe_dir, files$k450_benton %||% "HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt"), header = FALSE)
    ))
  } else {
    stop("Unsupported probe_filtering.array_type: ", opts$array_type, call. = FALSE)
  }
}

filter_good_beadcount_probes <- function(rgset, probes, bead_min_count = 3, bead_fail_fraction = 0.05) {
  nbeads <- extract_nbeads_matrix(rgset)
  probes <- intersect(probes, rownames(nbeads))
  nbeads <- nbeads[probes, , drop = FALSE]
  keep <- rowMeans(is.na(nbeads) | nbeads < bead_min_count, na.rm = TRUE) <= bead_fail_fraction
  names(keep)[keep]
}

filter_good_detection_probes <- function(rgset, probes, detection_p_threshold = 0.05, detection_p_fail_fraction = 0.01) {
  require_package("minfi")
  det_p <- minfi::detectionP(rgset)
  probes <- intersect(probes, rownames(det_p))
  det_p <- det_p[probes, , drop = FALSE]
  keep <- rowMeans(is.na(det_p) | det_p > detection_p_threshold, na.rm = TRUE) <= detection_p_fail_fraction
  names(keep)[keep]
}

record_probe_attrition <- function(config, chunk, attrition) {
  out_csv <- study_file(config, "probe_filter_attrition", chunk = chunk, ext = "csv")
  utils::write.csv(attrition, out_csv, row.names = FALSE)

  if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("scales", quietly = TRUE)) {
    attrition$filter <- factor(attrition$filter, levels = rev(attrition$filter))
    p <- ggplot2::ggplot(attrition) +
      ggplot2::geom_col(ggplot2::aes(filter, -probes_remaining), fill = "grey70", color = "black") +
      ggplot2::geom_col(ggplot2::aes(filter, probes_removed), fill = "darkred", color = "black") +
      ggplot2::geom_text(ggplot2::aes(filter, -min(probes_remaining) / 2, label = scales::comma(probes_remaining))) +
      ggplot2::geom_text(ggplot2::aes(filter, max(probes_removed) / 1.5, label = scales::comma(probes_removed))) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::ylab("") +
      ggplot2::xlab("") +
      ggplot2::scale_x_discrete(position = "top")
    ggplot2::ggsave(
      study_file(config, "probe_filter_attrition", chunk = chunk, ext = "pdf", dir = output_figure_dir(config)),
      plot = p,
      width = 14,
      height = 8
    )
  }

  invisible(out_csv)
}

add_attrition_step <- function(attrition, filter_name, before, after) {
  rbind(
    attrition,
    data.frame(
      filter = filter_name,
      probes_remaining = after,
      probes_removed = before - after,
      stringsAsFactors = FALSE
    )
  )
}

probe_filter_chunk <- function(config, chunk) {
  require_package("minfi")
  opts <- config$probe_filtering %||% list()
  grset <- readRDS(noob_grset_file(config, chunk))
  rgset <- readRDS(good_rgset_file(config, chunk))

  grset <- minfi::mapToGenome(grset)
  annotation <- minfi::getAnnotation(grset)
  current_n <- nrow(grset)
  attrition <- data.frame(
    filter = "Starting probes",
    probes_remaining = current_n,
    probes_removed = 0,
    stringsAsFactors = FALSE
  )

  if (isTRUE(opts$remove_snp_probes)) {
    before <- nrow(grset)
    grset <- minfi::dropLociWithSnps(grset, snps = c("CpG", "SBE"), maf = 0, snpAnno = NULL)
    attrition <- add_attrition_step(attrition, "Remove CpG/SBE SNP probes", before, nrow(grset))
  }

  if (isTRUE(opts$remove_cross_reactive) || isTRUE(opts$remove_polymorphic)) {
    before <- nrow(grset)
    cr_probes <- cross_reactive_probe_ids(config)
    if (!length(cr_probes)) {
      warning("No cross-reactive/polymorphic probe list files were found; skipping this filter.", call. = FALSE)
    } else {
      grset <- grset[!rownames(grset) %in% cr_probes, ]
    }
    attrition <- add_attrition_step(attrition, "Remove cross-reactive/polymorphic probes", before, nrow(grset))
  }

  if (isTRUE(opts$remove_xy)) {
    before <- nrow(grset)
    annotation <- minfi::getAnnotation(grset)
    chr_col <- intersect(c("chr", "seqnames"), colnames(annotation))[[1]]
    xy_probes <- rownames(annotation)[as.character(annotation[[chr_col]]) %in% c("chrX", "chrY", "X", "Y")]
    grset <- grset[!rownames(grset) %in% xy_probes, ]
    attrition <- add_attrition_step(attrition, "Remove XY probes", before, nrow(grset))
  }

  if (isTRUE(opts$remove_failed_beadcount)) {
    before <- nrow(grset)
    good <- filter_good_beadcount_probes(
      rgset,
      rownames(grset),
      bead_min_count = opts$bead_min_count %||% 3,
      bead_fail_fraction = opts$bead_fail_fraction %||% 0.05
    )
    grset <- grset[rownames(grset) %in% good, ]
    attrition <- add_attrition_step(attrition, "Remove beadcount-failed probes", before, nrow(grset))
  }

  if (isTRUE(opts$remove_failed_detection)) {
    before <- nrow(grset)
    good <- filter_good_detection_probes(
      rgset,
      rownames(grset),
      detection_p_threshold = opts$detection_p_threshold %||% 0.05,
      detection_p_fail_fraction = opts$detection_p_fail_fraction %||% 0.01
    )
    grset <- grset[rownames(grset) %in% good, ]
    attrition <- add_attrition_step(attrition, "Remove detection-P-failed probes", before, nrow(grset))
  }

  record_probe_attrition(config, chunk, attrition)
  saveRDS(grset, probe_filtered_grset_file(config, chunk))
  invisible(probe_filtered_grset_file(config, chunk))
}
