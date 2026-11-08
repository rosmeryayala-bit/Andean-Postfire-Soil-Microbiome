# Rarefaction and extrapolation
# Sample-coverage diversity curves using Hill numbers at genus level.
# Input: project-specific phyloseq objects, site metadata, or soil data tables.
# Output: analysis tables and figures as defined by the project directory settings.


rm(list = ls())


packages_cran <- c(
  "ggplot2", "dplyr", "tidyr", "RColorBrewer",
  "scales", "iNEXT", "cowplot"
)

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("phyloseq", quietly = TRUE)) {
  stop("The phyloseq package is required to run this script.")
}

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(scales)
library(iNEXT)
library(cowplot)

cat("\n=====================================================\n")
cat("RAREFACTION CURVES - GENUS - q = 0, 1, 2\n")
cat("RAREFACTION AND EXTRAPOLATION\n")
cat("=====================================================\n")
cat("Hora inicio:", as.character(Sys.time()), "\n")


find_existing_path <- function(paths) {
  paths <- path.expand(paths)
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) return(NA_character_)
  return(hit)
}

find_col <- function(names_vec, candidates) {
  hit <- candidates[candidates %in% names_vec][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  return(hit)
}

fix_duplicate_names <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- make.unique(colnames(df))
  return(df)
}

normalize_region <- function(x) {
  x <- trimws(as.character(x))
  ifelse(
    grepl("amazon", x, ignore.case = TRUE),
    "Amazonas",
    ifelse(
      grepl("cusco|cuzco", x, ignore.case = TRUE),
      "Cusco",
      x
    )
  )
}

normalize_treatment <- function(x) {
  x <- trimws(as.character(x))
  xl <- tolower(x)
  
  ifelse(
    xl %in% c("control", "ctrl", "unburned", "unburnt", "no burned", "no burn"),
    "Control",
    ifelse(
      xl %in% c("burned", "burnt", "wildfire", "fire", "incendio", "incendiado", "quemado"),
      "Wildfire",
      x
    )
  )
}

normalize_year <- function(x) {
  x <- trimws(as.character(x))
  xn <- suppressWarnings(as.numeric(x))
  ifelse(is.na(xn), x, as.character(as.integer(round(xn))))
}

clean_sample_label <- function(x) {
  x <- as.character(x)
  x <- sub("\\.raw$", "", x)
  x <- sub("_(BAC|ARQ|FUN)$", "", x)
  x <- sub("_(BAC|ARQ|FUN)\\.raw$", "", x)
  x <- gsub("_", "", x)
  return(x)
}

clean_genus <- function(x) {
  x <- trimws(as.character(x))
  
  x <- gsub("^g__", "", x)
  x <- gsub("^D_[0-9]+__", "", x)
  
  x[
    is.na(x) |
      x == "" |
      x == "NA" |
      x == "__" |
      grepl(
        "unclassified|uncultured|unknown|unidentified|metagenome|norank|ambiguous|incertae|chloroplast|mitochondria",
        x,
        ignore.case = TRUE
      )
  ] <- NA
  
  return(x)
}

make_empty_zero_table <- function() {
  data.frame(
    Region = character(),
    SampleCode = character(),
    Marker = character(),
    stringsAsFactors = FALSE
  )
}


order_amazonas <- c(
  "ACYU1", "ACYU2", "ACYU3",
  "AYU1", "AYU2", "AYU3",
  "ACKU1", "ACKU2", "ACKU3",
  "AKU1", "AKU2", "AKU3",
  "ACLE1", "ACLE2", "ACLE3",
  "ALE1", "ALE2", "ALE3",
  "ACTI1", "ACTI2", "ACTI3",
  "ATI1", "ATI2", "ATI3"
)

order_cusco <- c(
  "CCPAT1", "CCPAT2", "CCPAT3",
  "CPAT1", "CPAT2", "CPAT3",
  "CCPUC1", "CCPUC2", "CCPUC3",
  "CPUC1", "CPUC2", "CPUC3",
  "CCMAT1", "CCMAT2", "CCMAT3",
  "CMAT1", "CMAT2", "CMAT3",
  "CCHUC1", "CCHUC2", "CCHUC3",
  "CHUC1", "CHUC2", "CHUC3",
  "CCHAC1", "CCHAC2", "CCHAC3",
  "CHAC1", "CHAC2", "CHAC3"
)

sample_order_fixed <- c(order_amazonas, order_cusco)


extract_inext_coverage <- function(inext_obj) {
  
  df <- NULL
  
  if (is.list(inext_obj$iNextEst) &&
      "coverage_based" %in% names(inext_obj$iNextEst)) {
    df <- inext_obj$iNextEst$coverage_based
  }
  
  if (is.null(df) &&
      is.list(inext_obj$iNextEst) &&
      "size_based" %in% names(inext_obj$iNextEst)) {
    df <- inext_obj$iNextEst$size_based
  }
  
  if (is.null(df) && is.data.frame(inext_obj$iNextEst)) {
    df <- inext_obj$iNextEst
  }
  
  if (is.null(df)) {
    stop("No se pudo extraer la tabla de estimaciones de iNEXT.")
  }
  
  df <- fix_duplicate_names(df)
  names_df <- colnames(df)
  
  sample_col <- if ("Assemblage" %in% names_df) {
    "Assemblage"
  } else if ("site" %in% names_df) {
    "site"
  } else {
    stop("Sample identifier column not found in iNEXT.")
  }
  
  q_col <- if ("Order.q" %in% names_df) {
    "Order.q"
  } else if ("order" %in% names_df) {
    "order"
  } else if ("q" %in% names_df) {
    "q"
  } else {
    stop("Diversity order column not found in iNEXT.")
  }
  
  sc_col <- if ("SC" %in% names_df) {
    "SC"
  } else if ("SC.hat" %in% names_df) {
    "SC.hat"
  } else {
    stop("Sample coverage column not found in iNEXT.")
  }
  
  qd_col <- if ("qD" %in% names_df) {
    "qD"
  } else if ("qD.hat" %in% names_df) {
    "qD.hat"
  } else {
    stop("Diversity estimate column not found in iNEXT.")
  }
  
  method_col <- if ("Method" %in% names_df) {
    "Method"
  } else if ("method" %in% names_df) {
    "method"
  } else {
    NA_character_
  }
  
  out <- data.frame(
    SampleCode = as.character(df[[sample_col]]),
    q = as.numeric(as.character(df[[q_col]])),
    SC = as.numeric(df[[sc_col]]),
    qD = as.numeric(df[[qd_col]]),
    stringsAsFactors = FALSE
  )
  
  if (!is.na(method_col)) {
    out$Method <- as.character(df[[method_col]])
  } else {
    out$Method <- "iNEXT"
  }
  
  out <- out %>%
    filter(
      q %in% c(0, 1, 2),
      is.finite(SC),
      is.finite(qD)
    ) %>%
    mutate(
      SC = pmin(pmax(SC, 0), 1),
      q_lab = factor(
        paste0("q = ", q),
        levels = c("q = 0", "q = 1", "q = 2")
      )
    )
  
  return(out)
}


plot_region_rarefaction <- function(genus_mat,
                                    meta_df,
                                    region_value,
                                    region_order,
                                    marker_label) {
  
  meta_df <- fix_duplicate_names(meta_df)
  
  region_meta <- meta_df %>%
    filter(RegionPlot == region_value) %>%
    mutate(
      SampleLabel_chr = as.character(SampleLabel),
      OrderFixed = match(SampleLabel_chr, region_order)
    )
  
  missing_order <- region_meta %>%
    filter(is.na(OrderFixed)) %>%
    pull(SampleLabel_chr) %>%
    unique()
  
  if (length(missing_order) > 0) {
    cat("\nNotice:", region_value, "- samples not found in the specified order:\n")
    print(missing_order)
  }
  
  region_meta <- region_meta %>%
    mutate(
      OrderFixed2 = ifelse(
        is.na(OrderFixed),
        length(region_order) + row_number(),
        OrderFixed
      )
    ) %>%
    arrange(OrderFixed2)
  
  region_meta <- region_meta[
    region_meta$SampleID %in% colnames(genus_mat),
    ,
    drop = FALSE
  ]
  
  sample_ids <- region_meta$SampleID
  sample_codes <- region_meta$SampleLabel_chr
  
  assemblages_all <- list()
  
  for (i in seq_along(sample_ids)) {
    
    sid <- sample_ids[i]
    scode <- sample_codes[i]
    
    tmp <- genus_mat[, sid, drop = TRUE]
    tmp <- tmp[is.finite(tmp) & tmp > 0]
    tmp <- round(as.numeric(tmp))
    
    if (length(tmp) > 0) {
      names(tmp) <- rownames(genus_mat)[genus_mat[, sid, drop = TRUE] > 0]
    }
    
    assemblages_all[[scode]] <- tmp
  }
  
  total_reads_identified <- vapply(
    assemblages_all,
    function(x) sum(x, na.rm = TRUE),
    numeric(1)
  )
  
  zero_codes <- names(total_reads_identified)[total_reads_identified == 0]
  valid_assemblages <- assemblages_all[total_reads_identified > 0]
  
  if (length(zero_codes) > 0) {
    cat("\n", region_value, " - ", marker_label,
        ": samples without an identified genus; shown as a zero line:\n", sep = "")
    print(zero_codes)
  }
  
  est_df <- data.frame()
  
  if (length(valid_assemblages) > 0) {
    
    inext_obj <- iNEXT(
      valid_assemblages,
      q = c(0, 1, 2),
      datatype = "abundance",
      knots = 40,
      se = FALSE,
      nboot = 0
    )
    
    est_df <- extract_inext_coverage(inext_obj)
    
  } else {
    warning("No samples with identified genera for ", region_value, " - ", marker_label)
  }
  
  if (length(zero_codes) > 0) {
    
    zero_df <- expand.grid(
      SampleCode = zero_codes,
      q = c(0, 1, 2),
      SC = c(0, 1),
      stringsAsFactors = FALSE
    )
    
    zero_df$qD <- 0
    zero_df$Method <- "No identified genus"
    zero_df$q_lab <- factor(
      paste0("q = ", zero_df$q),
      levels = c("q = 0", "q = 1", "q = 2")
    )
    
    est_df <- bind_rows(est_df, zero_df)
  }
  
  est_df <- fix_duplicate_names(est_df)
  est_df$RegionPlot <- region_value
  
  sample_code_levels <- sample_codes[sample_codes %in% unique(est_df$SampleCode)]
  
  if (length(sample_code_levels) == 0) {
    sample_code_levels <- sample_codes
  }
  
  est_df$SampleCode <- factor(est_df$SampleCode, levels = sample_code_levels)
  
  pal_region <- setNames(
    scales::hue_pal()(length(sample_code_levels)),
    sample_code_levels
  )
  
  p_region <- ggplot(
    est_df,
    aes(
      x = SC,
      y = qD,
      color = SampleCode,
      group = SampleCode
    )
  ) +
    geom_line(linewidth = 0.75) +
    facet_wrap(
      ~ q_lab,
      nrow = 1,
      scales = "free_y"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    scale_color_manual(
      values = pal_region,
      breaks = sample_code_levels,
      name = "Sample code",
      drop = FALSE
    ) +
    labs(
      x = "Sample coverage",
      y = "Genus diversity"
    ) +
    theme_bw(base_size = 11) +
    theme(
      text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey85", color = "black"),
      strip.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(size = 7, color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  


  if (length(zero_codes) > 0) {
    zero_samples_df <- data.frame(
      Region = rep(region_value, length(zero_codes)),
      SampleCode = zero_codes,
      Marker = rep(marker_label, length(zero_codes)),
      stringsAsFactors = FALSE
    )
  } else {
    zero_samples_df <- make_empty_zero_table()
  }
  
  return(
    list(
      plot = p_region,
      estimates = est_df,
      zero_samples = zero_samples_df
    )
  )
}


run_rarefaction_marker <- function(marker_name,
                                   marker_label,
                                   base_candidates,
                                   ps_candidates,
                                   file_prefix) {
  
  cat("\n=====================================================\n")
  cat("PROCESANDO:", marker_label, "\n")
  cat("=====================================================\n")
  
  base <- find_existing_path(base_candidates)
  
  if (is.na(base)) {
    stop("Marker folder not found: ", marker_label)
  }
  
  ps_file <- find_existing_path(file.path(base, ps_candidates))
  
  if (is.na(ps_file)) {
    stop("Phyloseq object not found for: ", marker_label)
  }
  
  cat("\nCarpeta usada:\n")
  cat(base, "\n")
  
  cat("\nObjeto phyloseq usado:\n")
  cat(ps_file, "\n")
  
  ps <- readRDS(ps_file)
  
  cat("\nObjeto cargado:\n")
  cat("Muestras:", nsamples(ps), "\n")
  cat("ASVs:", ntaxa(ps), "\n")
  cat("Lecturas totales:", sum(sample_sums(ps)), "\n")
  
  meta_df <- data.frame(sample_data(ps), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(meta_df) <- sample_names(ps)
  meta_df <- fix_duplicate_names(meta_df)
  
  region_col <- find_col(
    colnames(meta_df),
    c("Region", "region", "Departamento", "departamento")
  )
  
  treat_col <- find_col(
    colnames(meta_df),
    c("Treatment", "treatment", "Tratamiento", "tratamiento", "Tratament")
  )
  
  year_col <- find_col(
    colnames(meta_df),
    c("Years_since_fire", "Years_since_fire.1", "YearsF", "Year", "year", "Año", "Anio")
  )
  
  if (is.na(region_col)) stop("Region column not found for: ", marker_label)
  if (is.na(treat_col))  stop("Treatment column not found for: ", marker_label)
  if (is.na(year_col))   stop("Years-since-fire column not found for: ", marker_label)
  
  meta_df$SampleID <- rownames(meta_df)
  meta_df$SampleLabel <- clean_sample_label(meta_df$SampleID)
  meta_df$RegionPlot <- normalize_region(meta_df[[region_col]])
  meta_df$TreatmentPlot <- normalize_treatment(meta_df[[treat_col]])
  meta_df$YearPlot <- normalize_year(meta_df[[year_col]])
  meta_df$YearNum <- suppressWarnings(as.numeric(meta_df$YearPlot))
  
  keep_samples <- meta_df$SampleID[
    meta_df$RegionPlot %in% c("Amazonas", "Cusco") &
      meta_df$TreatmentPlot %in% c("Control", "Wildfire")
  ]
  
  ps <- prune_samples(keep_samples, ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  meta_df <- meta_df[keep_samples, , drop = FALSE]
  meta_df <- fix_duplicate_names(meta_df)
  
  cat("\nSamples by region and treatment:\n")
  print(table(meta_df$RegionPlot, meta_df$TreatmentPlot))
  
  if (!"Genus" %in% rank_names(ps)) {
    stop("Taxonomy table lacks a Genus column for: ", marker_label)
  }
  
  tax_df <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  tax_df <- fix_duplicate_names(tax_df)
  tax_df$ASV_ID <- rownames(tax_df)
  tax_df$Genus_clean <- clean_genus(tax_df$Genus)
  
  valid_taxa <- tax_df$ASV_ID[!is.na(tax_df$Genus_clean)]
  
  if (length(valid_taxa) == 0) {
    stop("No identified genera for: ", marker_label)
  }
  
  ps_genus <- prune_taxa(valid_taxa, ps)
  
  tax_genus <- as.data.frame(tax_table(ps_genus), stringsAsFactors = FALSE)
  tax_genus <- fix_duplicate_names(tax_genus)
  tax_genus$ASV_ID <- rownames(tax_genus)
  tax_genus$Genus_clean <- clean_genus(tax_genus$Genus)
  
  cat("\n===== TAXONOMY SUMMARY =====\n")
  cat("ASVs totales:", ntaxa(ps), "\n")
  cat("ASVs assigned at genus level:", ntaxa(ps_genus), "\n")
  cat("Unique identified genera:", length(unique(na.omit(tax_genus$Genus_clean))), "\n")
  cat("Muestras conservadas:", nsamples(ps_genus), "\n")
  
  otu_mat <- as(otu_table(ps_genus), "matrix")
  
  if (!taxa_are_rows(ps_genus)) {
    otu_mat <- t(otu_mat)
  }
  
  tax_genus <- tax_genus[rownames(otu_mat), , drop = FALSE]
  
  genus_vector <- tax_genus$Genus_clean
  names(genus_vector) <- rownames(tax_genus)
  
  valid_rows <- !is.na(genus_vector)
  
  otu_valid <- otu_mat[valid_rows, , drop = FALSE]
  genus_valid <- genus_vector[valid_rows]
  
  genus_mat <- rowsum(
    otu_valid,
    group = genus_valid,
    reorder = TRUE
  )
  
  genus_mat <- genus_mat[, meta_df$SampleID, drop = FALSE]
  
  res_amz <- plot_region_rarefaction(
    genus_mat = genus_mat,
    meta_df = meta_df,
    region_value = "Amazonas",
    region_order = order_amazonas,
    marker_label = marker_label
  )
  
  res_cus <- plot_region_rarefaction(
    genus_mat = genus_mat,
    meta_df = meta_df,
    region_value = "Cusco",
    region_order = order_cusco,
    marker_label = marker_label
  )
  
  title_main <- ggdraw() +
    draw_label(
      paste0(marker_label, " - Genus-level rarefaction curves"),
      fontface = "bold",
      size = 18,
      x = 0.5,
      hjust = 0.5
    )
  
  title_amz <- ggdraw() +
    draw_label(
      "AMAZONAS",
      fontface = "bold",
      size = 15,
      x = 0,
      hjust = 0
    )
  
  title_cus <- ggdraw() +
    draw_label(
      "CUSCO",
      fontface = "bold",
      size = 15,
      x = 0,
      hjust = 0
    )
  
  block_amz <- plot_grid(
    title_amz,
    res_amz$plot,
    ncol = 1,
    rel_heights = c(0.07, 1)
  )
  
  block_cus <- plot_grid(
    title_cus,
    res_cus$plot,
    ncol = 1,
    rel_heights = c(0.07, 1)
  )
  
  fig_final <- plot_grid(
    title_main,
    block_amz,
    block_cus,
    ncol = 1,
    rel_heights = c(0.08, 1, 1)
  )
  
  print(fig_final)
  
  outdir <- file.path(
    base,
    "results/plots",
    paste0("rarefaction_curves_", file_prefix, "_genus_q012_coverage")
  )
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  estimates_all <- bind_rows(
    res_amz$estimates,
    res_cus$estimates
  )
  
  zero_samples_all <- bind_rows(
    res_amz$zero_samples,
    res_cus$zero_samples
  )
  
  write.csv(
    estimates_all,
    file.path(outdir, paste0(file_prefix, "_rarefaction_genus_q012_coverage_estimates.csv")),
    row.names = FALSE
  )
  
  write.csv(
    zero_samples_all,
    file.path(outdir, paste0(file_prefix, "_samples_without_identified_genus.csv")),
    row.names = FALSE
  )
  
  write.csv(
    genus_mat,
    file.path(outdir, paste0(file_prefix, "_genus_count_matrix_used_for_rarefaction.csv"))
  )
  
  write.csv(
    tax_genus[!is.na(tax_genus$Genus_clean), ],
    file.path(outdir, paste0(file_prefix, "_ASVs_with_identified_genus_used.csv")),
    row.names = FALSE
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_rarefaction_genus_q012_coverage.png")),
    fig_final,
    width = 17,
    height = 12,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_rarefaction_genus_q012_coverage.pdf")),
    fig_final,
    width = 17,
    height = 12,
    bg = "white"
  )
  
  cat("\nArchivos guardados en:\n")
  cat(outdir, "\n")
  
  cat("\nFigura PNG:\n")
  cat(file.path(outdir, paste0(file_prefix, "_rarefaction_genus_q012_coverage.png")), "\n")
  
  cat("\nFigura PDF:\n")
  cat(file.path(outdir, paste0(file_prefix, "_rarefaction_genus_q012_coverage.pdf")), "\n")
  
  return(fig_final)
}


p_bac_rarefaction <- run_rarefaction_marker(
  marker_name = "Bacteria",
  marker_label = "Bacterias 16S",
  base_candidates = c(
    "~/Rosmery/Tesis Rosmery-Metabarcoding/RESULT BACTERIA",
    "~/Tesis Rosmery-Metabarcoding/RESULT BACTERIA"
  ),
  ps_candidates = c(
    "data/processed_data/phyloseq/ps_sin_prevalencia.rds",
    "data/processed_data/phyloseq/ps_filt.rds",
    "data/processed_data/phyloseq/ps_noplast.rds",
    "data/processed_data/phyloseq/ps_raw.rds"
  ),
  file_prefix = "BAC"
)

p_arq_rarefaction <- run_rarefaction_marker(
  marker_name = "Archaea",
  marker_label = "Arqueas 16S",
  base_candidates = c(
    "~/Rosmery/Tesis Rosmery-Metabarcoding/RESULT ARQUEA",
    "~/Tesis Rosmery-Metabarcoding/RESULT ARQUEA"
  ),
  ps_candidates = c(
    "data/processed_data/phyloseq/ps_sin_prevalencia_ARQ.rds",
    "data/processed_data/phyloseq/ps_only_archaea_ARQ.rds",
    "data/processed_data/phyloseq/ps_raw_ARQ.rds"
  ),
  file_prefix = "ARQ"
)

p_its_rarefaction <- run_rarefaction_marker(
  marker_name = "Fungi_ITS",
  marker_label = "Hongos ITS",
  base_candidates = c(
    "~/Rosmery/Tesis Rosmery-Metabarcoding/RESULT FUNGAL",
    "~/Tesis Rosmery-Metabarcoding/RESULT FUNGAL"
  ),
  ps_candidates = c(
    "data/processed_data/phyloseq/ps_raw_ITS.rds",
    "data/processed_data/phyloseq/ps_filt_ITS.rds"
  ),
  file_prefix = "ITS"
)

cat("\n=====================================================\n")
cat("Rarefaction curves generated\n")
cat("Marcadores: Bacterias 16S, Arqueas 16S, Hongos ITS\n")
cat("Taxonomic rank: Genus\n")
cat("q = 0, q = 1, q = 2\n")
cat("Eje X: Sample coverage\n")
cat("Hora fin:", as.character(Sys.time()), "\n")
cat("=====================================================\n")


rm(list = ls())


packages_cran <- c(
  "ggplot2", "dplyr", "tidyr", "RColorBrewer",
  "scales", "iNEXT", "cowplot"
)

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("phyloseq", quietly = TRUE)) {
  stop("The phyloseq package is required to run this script.")
}

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(scales)
library(iNEXT)
library(cowplot)

cat("\n=====================================================\n")
cat("RAREFACTION CURVES BY MARKER - ASV LEVEL\n")
cat("q = 0, q = 1, q = 2 | Eje X: Sample coverage\n")
cat("=====================================================\n")
cat("Hora inicio:", as.character(Sys.time()), "\n")


MAKE_GROUPED_PAPER <- TRUE
MAKE_SAMPLE_SUPPLEMENT <- TRUE
KNOTS_INEXT <- 80


find_existing_path <- function(paths) {
  paths <- path.expand(paths)
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) return(NA_character_)
  return(hit)
}

fix_duplicate_names <- function(df) {
  rn <- rownames(df)
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- make.unique(colnames(df))
  rownames(df) <- rn
  return(df)
}

find_col <- function(names_vec, candidates) {
  hit <- candidates[candidates %in% names_vec][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  return(hit)
}

normalize_region <- function(x) {
  x <- trimws(as.character(x))
  ifelse(
    grepl("amazon", x, ignore.case = TRUE),
    "Amazonas",
    ifelse(
      grepl("cusco|cuzco", x, ignore.case = TRUE),
      "Cusco",
      x
    )
  )
}

normalize_treatment <- function(x) {
  x <- trimws(as.character(x))
  xl <- tolower(x)
  
  ifelse(
    xl %in% c("control", "ctrl", "unburned", "unburnt", "no burned", "no burn"),
    "Control",
    ifelse(
      xl %in% c("burned", "burnt", "wildfire", "fire", "incendio", "incendiado", "quemado"),
      "Wildfire",
      x
    )
  )
}

normalize_year <- function(x) {
  x <- trimws(as.character(x))
  xn <- suppressWarnings(as.numeric(x))
  ifelse(is.na(xn), x, as.character(as.integer(round(xn))))
}

clean_sample_label <- function(x) {
  x <- as.character(x)
  x <- sub("\\.raw$", "", x)
  x <- sub("_(BAC|ARQ|FUN)$", "", x)
  x <- sub("_(BAC|ARQ|FUN)\\.raw$", "", x)
  x <- gsub("_", "", x)
  return(x)
}

remove_contaminant_taxa <- function(ps) {
  
  tax_df <- tryCatch(
    as.data.frame(tax_table(ps), stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(tax_df)) {
    return(ps)
  }
  
  tax_df <- fix_duplicate_names(tax_df)
  
  bad_taxa <- apply(
    tax_df,
    1,
    function(z) {
      any(
        grepl(
          "chloroplast|mitochondria",
          paste(z, collapse = " "),
          ignore.case = TRUE
        )
      )
    }
  )
  
  bad_taxa[is.na(bad_taxa)] <- FALSE
  
  if (any(bad_taxa)) {
    cat("\nTaxones removidos por cloroplasto/mitocondria:", sum(bad_taxa), "\n")
    ps <- prune_taxa(!bad_taxa, ps)
  }
  
  return(ps)
}

make_empty_estimates <- function() {
  data.frame(
    Assemblage = character(),
    q = numeric(),
    SC = numeric(),
    qD = numeric(),
    Method = character(),
    q_lab = factor(
      character(),
      levels = c("q = 0", "q = 1", "q = 2")
    ),
    stringsAsFactors = FALSE
  )
}

make_empty_zero_table <- function() {
  data.frame(
    Marker = character(),
    Type = character(),
    Assemblage = character(),
    stringsAsFactors = FALSE
  )
}


order_amazonas <- c(
  "ACYU1", "ACYU2", "ACYU3",
  "AYU1", "AYU2", "AYU3",
  "ACKU1", "ACKU2", "ACKU3",
  "AKU1", "AKU2", "AKU3",
  "ACLE1", "ACLE2", "ACLE3",
  "ALE1", "ALE2", "ALE3",
  "ACTI1", "ACTI2", "ACTI3",
  "ATI1", "ATI2", "ATI3"
)

order_cusco <- c(
  "CCPAT1", "CCPAT2", "CCPAT3",
  "CPAT1", "CPAT2", "CPAT3",
  "CCPUC1", "CCPUC2", "CCPUC3",
  "CPUC1", "CPUC2", "CPUC3",
  "CCMAT1", "CCMAT2", "CCMAT3",
  "CMAT1", "CMAT2", "CMAT3",
  "CCHUC1", "CCHUC2", "CCHUC3",
  "CHUC1", "CHUC2", "CHUC3",
  "CCHAC1", "CCHAC2", "CCHAC3",
  "CHAC1", "CHAC2", "CHAC3"
)

sample_order_fixed <- c(order_amazonas, order_cusco)

sample_palette_fixed <- setNames(
  colorRampPalette(
    c(
      brewer.pal(12, "Paired"),
      brewer.pal(8, "Dark2"),
      brewer.pal(12, "Set3"),
      brewer.pal(8, "Set2"),
      brewer.pal(9, "Set1")
    )
  )(length(sample_order_fixed)),
  sample_order_fixed
)

group_palette <- c(
  "Amazonas - Control"  = "#1B9E77",
  "Amazonas - Wildfire" = "#D95F02",
  "Cusco - Control"     = "#7570B3",
  "Cusco - Wildfire"    = "#E7298A"
)


extract_inext_coverage <- function(inext_obj) {
  
  df <- NULL
  
  if (is.list(inext_obj$iNextEst) &&
      "coverage_based" %in% names(inext_obj$iNextEst)) {
    df <- inext_obj$iNextEst$coverage_based
  }
  
  if (is.null(df) &&
      is.list(inext_obj$iNextEst) &&
      "size_based" %in% names(inext_obj$iNextEst)) {
    df <- inext_obj$iNextEst$size_based
  }
  
  if (is.null(df) && is.data.frame(inext_obj$iNextEst)) {
    df <- inext_obj$iNextEst
  }
  
  if (is.null(df)) {
    stop("No se pudo extraer la tabla de estimaciones de iNEXT.")
  }
  
  df <- fix_duplicate_names(df)
  names_df <- colnames(df)
  
  assemblage_col <- if ("Assemblage" %in% names_df) {
    "Assemblage"
  } else if ("site" %in% names_df) {
    "site"
  } else {
    stop("Assemblage/site column not found in iNEXT.")
  }
  
  q_col <- if ("Order.q" %in% names_df) {
    "Order.q"
  } else if ("order" %in% names_df) {
    "order"
  } else if ("q" %in% names_df) {
    "q"
  } else {
    stop("Diversity order column not found in iNEXT.")
  }
  
  sc_col <- if ("SC" %in% names_df) {
    "SC"
  } else if ("SC.hat" %in% names_df) {
    "SC.hat"
  } else {
    stop("Sample coverage column not found in iNEXT.")
  }
  
  qd_col <- if ("qD" %in% names_df) {
    "qD"
  } else if ("qD.hat" %in% names_df) {
    "qD.hat"
  } else {
    stop("Diversity estimate column not found in iNEXT.")
  }
  
  method_col <- if ("Method" %in% names_df) {
    "Method"
  } else if ("method" %in% names_df) {
    "method"
  } else {
    NA_character_
  }
  
  out <- data.frame(
    Assemblage = as.character(df[[assemblage_col]]),
    q = as.numeric(as.character(df[[q_col]])),
    SC = as.numeric(df[[sc_col]]),
    qD = as.numeric(df[[qd_col]]),
    stringsAsFactors = FALSE
  )
  
  if (!is.na(method_col)) {
    out$Method <- as.character(df[[method_col]])
  } else {
    out$Method <- "iNEXT"
  }
  
  out <- out %>%
    filter(
      q %in% c(0, 1, 2),
      is.finite(SC),
      is.finite(qD)
    ) %>%
    mutate(
      SC = pmin(pmax(SC, 0), 1),
      q_lab = factor(
        paste0("q = ", q),
        levels = c("q = 0", "q = 1", "q = 2")
      )
    )
  
  return(out)
}

make_inext_estimates <- function(assemblage_list,
                                 marker_label,
                                 assemblage_type) {
  
  total_reads <- vapply(
    assemblage_list,
    function(x) sum(x, na.rm = TRUE),
    numeric(1)
  )
  
  zero_names <- names(total_reads)[total_reads <= 0 | is.na(total_reads)]
  valid_list <- assemblage_list[total_reads > 0]
  
  est_df <- make_empty_estimates()
  
  if (length(valid_list) > 0) {
    
    inext_obj <- iNEXT(
      valid_list,
      q = c(0, 1, 2),
      datatype = "abundance",
      knots = KNOTS_INEXT,
      se = FALSE,
      nboot = 0
    )
    
    est_df <- extract_inext_coverage(inext_obj)
  }
  
  zero_table <- make_empty_zero_table()
  
  if (length(zero_names) > 0) {
    
    cat("\n", marker_label, " - ", assemblage_type,
        ": groups/samples without valid reads; shown as a zero line:\n", sep = "")
    print(zero_names)
    
    zero_df <- expand.grid(
      Assemblage = zero_names,
      q = c(0, 1, 2),
      SC = seq(0, 1, length.out = KNOTS_INEXT),
      stringsAsFactors = FALSE
    )
    
    zero_df$qD <- 0
    zero_df$Method <- "No reads"
    zero_df$q_lab <- factor(
      paste0("q = ", zero_df$q),
      levels = c("q = 0", "q = 1", "q = 2")
    )
    
    est_df <- bind_rows(est_df, zero_df)
    
    zero_table <- data.frame(
      Marker = rep(marker_label, length(zero_names)),
      Type = rep(assemblage_type, length(zero_names)),
      Assemblage = zero_names,
      stringsAsFactors = FALSE
    )
  }
  
  return(
    list(
      estimates = est_df,
      zero_table = zero_table
    )
  )
}


prepare_marker_asv <- function(marker_name,
                               marker_label,
                               base_candidates,
                               ps_candidates,
                               file_prefix) {
  
  cat("\n=====================================================\n")
  cat("PREPARANDO:", marker_label, "\n")
  cat("=====================================================\n")
  
  base <- find_existing_path(base_candidates)
  
  if (is.na(base)) {
    stop("Marker folder not found: ", marker_label)
  }
  
  ps_file <- find_existing_path(file.path(base, ps_candidates))
  
  if (is.na(ps_file)) {
    stop("Phyloseq object not found for: ", marker_label)
  }
  
  cat("\nCarpeta usada:\n")
  cat(base, "\n")
  
  cat("\nObjeto phyloseq usado:\n")
  cat(ps_file, "\n")
  
  ps <- readRDS(ps_file)
  
  cat("\nObjeto cargado:\n")
  cat("Muestras:", nsamples(ps), "\n")
  cat("ASVs iniciales:", ntaxa(ps), "\n")
  cat("Lecturas iniciales:", sum(sample_sums(ps)), "\n")
  

  meta_df <- data.frame(sample_data(ps), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(meta_df) <- sample_names(ps)
  meta_df <- fix_duplicate_names(meta_df)
  
  region_col <- find_col(
    colnames(meta_df),
    c("Region", "region", "Departamento", "departamento")
  )
  
  treat_col <- find_col(
    colnames(meta_df),
    c("Treatment", "treatment", "Tratamiento", "tratamiento", "Tratament")
  )
  
  year_col <- find_col(
    colnames(meta_df),
    c("Years_since_fire", "Years_since_fire.1", "YearsF", "Year", "year", "Año", "Anio")
  )
  
  if (is.na(region_col)) stop("Region column not found for: ", marker_label)
  if (is.na(treat_col))  stop("Treatment column not found for: ", marker_label)
  if (is.na(year_col))   stop("Years-since-fire column not found for: ", marker_label)
  
  meta_df$SampleID <- rownames(meta_df)
  meta_df$SampleLabel <- clean_sample_label(meta_df$SampleID)
  meta_df$RegionPlot <- normalize_region(meta_df[[region_col]])
  meta_df$TreatmentPlot <- normalize_treatment(meta_df[[treat_col]])
  meta_df$YearPlot <- normalize_year(meta_df[[year_col]])
  meta_df$YearNum <- suppressWarnings(as.numeric(meta_df$YearPlot))
  
  keep_samples <- meta_df$SampleID[
    meta_df$RegionPlot %in% c("Amazonas", "Cusco") &
      meta_df$TreatmentPlot %in% c("Control", "Wildfire")
  ]
  
  ps <- prune_samples(keep_samples, ps)
  ps <- remove_contaminant_taxa(ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  meta_df <- meta_df[keep_samples, , drop = FALSE]
  meta_df <- fix_duplicate_names(meta_df)
  
  cat("\nSamples by region and treatment:\n")
  print(table(meta_df$RegionPlot, meta_df$TreatmentPlot))
  
  otu_mat <- as(otu_table(ps), "matrix")
  
  if (!taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }
  
  rownames(otu_mat) <- make.unique(rownames(otu_mat))
  
  missing_cols <- setdiff(meta_df$SampleID, colnames(otu_mat))
  
  if (length(missing_cols) > 0) {
    zero_add <- matrix(
      0,
      nrow = nrow(otu_mat),
      ncol = length(missing_cols),
      dimnames = list(rownames(otu_mat), missing_cols)
    )
    otu_mat <- cbind(otu_mat, zero_add)
  }
  
  otu_mat <- otu_mat[, meta_df$SampleID, drop = FALSE]
  
  cat("\nRarefaction dataset:\n")
  cat("Muestras:", ncol(otu_mat), "\n")
  cat("ASVs:", nrow(otu_mat), "\n")
  cat("Lecturas:", sum(otu_mat), "\n")
  
  return(
    list(
      marker_name = marker_name,
      marker_label = marker_label,
      file_prefix = file_prefix,
      base = base,
      ps_file = ps_file,
      meta_df = meta_df,
      otu_mat = otu_mat
    )
  )
}


make_sample_assemblages <- function(otu_mat, meta_df) {
  
  meta_ordered <- meta_df %>%
    mutate(
      SampleLabel_chr = as.character(SampleLabel),
      OrderFixed = match(SampleLabel_chr, sample_order_fixed),
      OrderFixed2 = ifelse(
        is.na(OrderFixed),
        length(sample_order_fixed) + row_number(),
        OrderFixed
      )
    ) %>%
    arrange(OrderFixed2)
  
  assemblages <- list()
  
  for (i in seq_len(nrow(meta_ordered))) {
    
    sid <- meta_ordered$SampleID[i]
    scode <- meta_ordered$SampleLabel_chr[i]
    
    v <- otu_mat[, sid, drop = TRUE]
    names(v) <- rownames(otu_mat)
    v <- round(as.numeric(v))
    names(v) <- rownames(otu_mat)
    v <- v[is.finite(v) & v > 0]
    
    assemblages[[scode]] <- v
  }
  
  return(
    list(
      assemblages = assemblages,
      meta_ordered = meta_ordered
    )
  )
}

make_group_assemblages <- function(otu_mat, meta_df) {
  
  group_info <- data.frame(
    RegionPlot = c("Amazonas", "Amazonas", "Cusco", "Cusco"),
    TreatmentPlot = c("Control", "Wildfire", "Control", "Wildfire"),
    Assemblage = c(
      "Amazonas - Control",
      "Amazonas - Wildfire",
      "Cusco - Control",
      "Cusco - Wildfire"
    ),
    stringsAsFactors = FALSE
  )
  
  assemblages <- list()
  
  for (i in seq_len(nrow(group_info))) {
    
    region_use <- group_info$RegionPlot[i]
    treatment_use <- group_info$TreatmentPlot[i]
    assemblage_name <- group_info$Assemblage[i]
    
    sample_ids <- meta_df$SampleID[
      meta_df$RegionPlot == region_use &
        meta_df$TreatmentPlot == treatment_use
    ]
    
    sample_ids <- intersect(sample_ids, colnames(otu_mat))
    
    if (length(sample_ids) == 0) {
      v <- numeric(0)
    } else {
      v <- rowSums(otu_mat[, sample_ids, drop = FALSE], na.rm = TRUE)
      names(v) <- rownames(otu_mat)
      v <- round(as.numeric(v))
      names(v) <- rownames(otu_mat)
      v <- v[is.finite(v) & v > 0]
    }
    
    assemblages[[assemblage_name]] <- v
  }
  
  return(
    list(
      assemblages = assemblages,
      group_info = group_info
    )
  )
}


plot_grouped_paper <- function(marker_data) {
  
  marker_label <- marker_data$marker_label
  file_prefix <- marker_data$file_prefix
  base <- marker_data$base
  otu_mat <- marker_data$otu_mat
  meta_df <- marker_data$meta_df
  
  group_obj <- make_group_assemblages(otu_mat, meta_df)
  
  inext_res <- make_inext_estimates(
    assemblage_list = group_obj$assemblages,
    marker_label = marker_label,
    assemblage_type = "Region_Treatment"
  )
  
  est_df <- inext_res$estimates %>%
    left_join(group_obj$group_info, by = "Assemblage")
  
  est_df$Assemblage <- factor(
    est_df$Assemblage,
    levels = group_obj$group_info$Assemblage
  )
  
  p <- ggplot(
    est_df,
    aes(
      x = SC,
      y = qD,
      color = Assemblage,
      group = Assemblage
    )
  ) +
    geom_line(linewidth = 1.1) +
    facet_wrap(
      ~ q_lab,
      nrow = 1,
      scales = "free_y"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    scale_color_manual(
      values = group_palette,
      name = "Region × treatment",
      drop = FALSE
    ) +
    labs(
      title = paste0(marker_label, " - Coverage-based rarefaction/extrapolation"),
      subtitle = "ASV-level Hill diversity grouped by region and treatment",
      x = "Sample coverage",
      y = "ASV diversity"
    ) +
    theme_bw(base_size = 12) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      strip.background = element_rect(fill = "grey85", color = "black"),
      strip.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(size = 9, color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(p)
  
  outdir <- file.path(
    base,
    "results/plots",
    paste0("rarefaction_", file_prefix, "_ASV_q012_coverage_PAPER_grouped")
  )
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    est_df,
    file.path(outdir, paste0(file_prefix, "_ASV_rarefaction_grouped_estimates.csv")),
    row.names = FALSE
  )
  
  write.csv(
    inext_res$zero_table,
    file.path(outdir, paste0(file_prefix, "_zero_groups.csv")),
    row.names = FALSE
  )
  
  write.csv(
    otu_mat,
    file.path(outdir, paste0(file_prefix, "_ASV_count_matrix_used.csv"))
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_ASV_rarefaction_q012_coverage_PAPER_grouped.png")),
    p,
    width = 13.5,
    height = 5.5,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_ASV_rarefaction_q012_coverage_PAPER_grouped.pdf")),
    p,
    width = 13.5,
    height = 5.5,
    bg = "white"
  )
  
  cat("\nFigura recomendada paper guardada en:\n")
  cat(outdir, "\n")
  
  return(p)
}


plot_sample_supplement <- function(marker_data) {
  
  marker_label <- marker_data$marker_label
  file_prefix <- marker_data$file_prefix
  base <- marker_data$base
  otu_mat <- marker_data$otu_mat
  meta_df <- marker_data$meta_df
  
  sample_obj <- make_sample_assemblages(otu_mat, meta_df)
  
  inext_res <- make_inext_estimates(
    assemblage_list = sample_obj$assemblages,
    marker_label = marker_label,
    assemblage_type = "Sample"
  )
  
  sample_meta <- sample_obj$meta_ordered %>%
    select(
      SampleID,
      SampleLabel,
      SampleLabel_chr,
      RegionPlot,
      TreatmentPlot,
      YearPlot,
      YearNum
    )
  
  est_df <- inext_res$estimates %>%
    left_join(
      sample_meta,
      by = c("Assemblage" = "SampleLabel_chr")
    )
  
  sample_levels <- sample_meta$SampleLabel_chr[
    sample_meta$SampleLabel_chr %in% unique(est_df$Assemblage)
  ]
  
  est_df$Assemblage <- factor(est_df$Assemblage, levels = sample_levels)
  est_df$RegionPlot <- factor(est_df$RegionPlot, levels = c("Amazonas", "Cusco"))
  est_df$TreatmentPlot <- factor(est_df$TreatmentPlot, levels = c("Control", "Wildfire"))
  
  sample_colors_use <- sample_palette_fixed[sample_levels]
  
  missing_colors <- sample_levels[is.na(sample_colors_use)]
  
  if (length(missing_colors) > 0) {
    extra_cols <- setNames(
      colorRampPalette(brewer.pal(8, "Dark2"))(length(missing_colors)),
      missing_colors
    )
    sample_colors_use[missing_colors] <- extra_cols
  }
  
  p <- ggplot(
    est_df,
    aes(
      x = SC,
      y = qD,
      color = Assemblage,
      group = Assemblage
    )
  ) +
    geom_line(linewidth = 0.65, alpha = 0.95) +
    facet_grid(
      RegionPlot ~ q_lab,
      scales = "free_y"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    scale_color_manual(
      values = sample_colors_use,
      breaks = sample_levels,
      name = "Sample code",
      drop = FALSE
    ) +
    labs(
      title = paste0(marker_label, " - Coverage-based rarefaction/extrapolation by sample"),
      subtitle = "ASV-level Hill diversity; fixed sample order and fixed sample colors",
      x = "Sample coverage",
      y = "ASV diversity"
    ) +
    guides(
      color = guide_legend(
        ncol = 2,
        byrow = TRUE,
        override.aes = list(linewidth = 1.2, alpha = 1)
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      strip.background = element_rect(fill = "grey85", color = "black"),
      strip.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(size = 6.5, color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(p)
  
  outdir <- file.path(
    base,
    "results/plots",
    paste0("rarefaction_", file_prefix, "_ASV_q012_coverage_SUPPLEMENT_sample_colors")
  )
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    est_df,
    file.path(outdir, paste0(file_prefix, "_ASV_rarefaction_sample_estimates.csv")),
    row.names = FALSE
  )
  
  write.csv(
    inext_res$zero_table,
    file.path(outdir, paste0(file_prefix, "_zero_samples.csv")),
    row.names = FALSE
  )
  
  write.csv(
    sample_meta,
    file.path(outdir, paste0(file_prefix, "_sample_order_used.csv")),
    row.names = FALSE
  )
  
  write.csv(
    otu_mat,
    file.path(outdir, paste0(file_prefix, "_ASV_count_matrix_used.csv"))
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_ASV_rarefaction_q012_coverage_SUPPLEMENT_sample_colors.png")),
    p,
    width = 17.5,
    height = 9.5,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_ASV_rarefaction_q012_coverage_SUPPLEMENT_sample_colors.pdf")),
    p,
    width = 17.5,
    height = 9.5,
    bg = "white"
  )
  
  cat("\nFigura suplementaria por muestra guardada en:\n")
  cat(outdir, "\n")
  
  return(p)
}


run_marker_rarefaction <- function(marker_name,
                                   marker_label,
                                   base_candidates,
                                   ps_candidates,
                                   file_prefix) {
  
  marker_data <- prepare_marker_asv(
    marker_name = marker_name,
    marker_label = marker_label,
    base_candidates = base_candidates,
    ps_candidates = ps_candidates,
    file_prefix = file_prefix
  )
  
  p_grouped <- NULL
  p_sample <- NULL
  
  if (MAKE_GROUPED_PAPER) {
    p_grouped <- plot_grouped_paper(marker_data)
  }
  
  if (MAKE_SAMPLE_SUPPLEMENT) {
    p_sample <- plot_sample_supplement(marker_data)
  }
  
  return(
    list(
      grouped_paper = p_grouped,
      sample_supplement = p_sample
    )
  )
}


plots_bac <- run_marker_rarefaction(
  marker_name = "Bacteria",
  marker_label = "Bacterias 16S",
  base_candidates = c(
    "~/Rosmery/Tesis Rosmery-Metabarcoding/RESULT BACTERIA",
    "~/Tesis Rosmery-Metabarcoding/RESULT BACTERIA"
  ),
  ps_candidates = c(
    "data/processed_data/phyloseq/ps_sin_prevalencia.rds",
    "data/processed_data/phyloseq/ps_filt.rds",
    "data/processed_data/phyloseq/ps_noplast.rds",
    "data/processed_data/phyloseq/ps_raw.rds"
  ),
  file_prefix = "BAC"
)

plots_arq <- run_marker_rarefaction(
  marker_name = "Archaea",
  marker_label = "Arqueas 16S",
  base_candidates = c(
    "~/Rosmery/Tesis Rosmery-Metabarcoding/RESULT ARQUEA",
    "~/Tesis Rosmery-Metabarcoding/RESULT ARQUEA"
  ),
  ps_candidates = c(
    "data/processed_data/phyloseq/ps_sin_prevalencia_ARQ.rds",
    "data/processed_data/phyloseq/ps_only_archaea_ARQ.rds",
    "data/processed_data/phyloseq/ps_raw_ARQ.rds"
  ),
  file_prefix = "ARQ"
)

plots_its <- run_marker_rarefaction(
  marker_name = "Fungi_ITS",
  marker_label = "Hongos ITS",
  base_candidates = c(
    "~/Rosmery/Tesis Rosmery-Metabarcoding/RESULT FUNGAL",
    "~/Tesis Rosmery-Metabarcoding/RESULT FUNGAL"
  ),
  ps_candidates = c(
    "data/processed_data/phyloseq/ps_raw_ITS.rds",
    "data/processed_data/phyloseq/ps_filt_ITS.rds"
  ),
  file_prefix = "ITS"
)

cat("\n=====================================================\n")
cat("Rarefaction curves generated\n")
cat("Nivel: ASV\n")
cat("Marcadores: Bacterias 16S, Arqueas 16S, Hongos ITS\n")
cat("Outputs: grouped curves and sample-level curves\n")
cat("Hora fin:", as.character(Sys.time()), "\n")
cat("=====================================================\n")
