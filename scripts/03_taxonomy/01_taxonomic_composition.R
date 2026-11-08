# Taxonomic composition
# Bacterial 16S, archaeal 16S and fungal ITS2 community composition across sampling groups.
# Input: project-specific phyloseq objects, site metadata, or soil data tables.
# Output: analysis tables and figures as defined by the project directory settings.


rm(list = ls())


packages_cran <- c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "scales")

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

cat("\n=====================================================\n")
cat("TAXONOMIC COMPOSITION - IDENTIFIED GENERA ONLY\n")
cat("SPECIFIED SAMPLE ORDER\n")
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

make_palette <- function(n) {
  
  base_cols <- c(
    brewer.pal(12, "Set3"),
    brewer.pal(8, "Set2"),
    brewer.pal(8, "Dark2"),
    brewer.pal(12, "Paired"),
    brewer.pal(9, "Set1")
  )
  
  base_cols <- unique(base_cols)
  
  if (n <= length(base_cols)) {
    return(base_cols[seq_len(n)])
  } else {
    return(colorRampPalette(base_cols)(n))
  }
}


plot_genus_composition_marker <- function(marker_name,
                                          marker_label,
                                          base_candidates,
                                          ps_candidates,
                                          out_folder_name,
                                          file_prefix,
                                          top_n = 30) {
  
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
  
  cat("\nSamples by region and treatment:\n")
  print(table(meta_df$RegionPlot, meta_df$TreatmentPlot))
  


  
  if (!"Genus" %in% rank_names(ps)) {
    stop("Taxonomy table lacks a Genus column for: ", marker_label)
  }
  
  tax_df <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  tax_df$ASV_ID <- rownames(tax_df)
  tax_df$Genus_clean <- clean_genus(tax_df$Genus)
  
  valid_taxa <- tax_df$ASV_ID[!is.na(tax_df$Genus_clean)]
  
  if (length(valid_taxa) == 0) {
    stop("No identified genera found for: ", marker_label)
  }
  
  ps_genus <- prune_taxa(valid_taxa, ps)
  ps_genus <- prune_samples(sample_sums(ps_genus) > 0, ps_genus)
  
  tax_genus <- as.data.frame(tax_table(ps_genus), stringsAsFactors = FALSE)
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
  
  sample_totals <- colSums(genus_mat, na.rm = TRUE)
  
  genus_rel <- sweep(genus_mat, 2, sample_totals, "/")
  genus_rel[!is.finite(genus_rel)] <- 0
  


  
  df_long <- as.data.frame(as.table(genus_rel), stringsAsFactors = FALSE)
  colnames(df_long) <- c("Genus", "SampleID", "Abundance")
  
  df_long$Abundance <- as.numeric(df_long$Abundance)
  
  meta_plot <- meta_df[colnames(genus_rel), , drop = FALSE]
  meta_plot$SampleID <- rownames(meta_plot)
  
  df_long <- df_long %>%
    left_join(
      meta_plot %>%
        select(SampleID, SampleLabel, RegionPlot, TreatmentPlot, YearPlot, YearNum),
      by = "SampleID"
    )
  


  
  top_genera <- df_long %>%
    group_by(Genus) %>%
    summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalAbundance)) %>%
    slice_head(n = top_n) %>%
    pull(Genus)
  
  df_long$GenusPlot <- ifelse(df_long$Genus %in% top_genera, df_long$Genus, "Other")
  
  df_plot <- df_long %>%
    group_by(SampleID, SampleLabel, RegionPlot, TreatmentPlot, YearPlot, YearNum, GenusPlot) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  

  df_plot <- df_plot %>%
    group_by(SampleID) %>%
    mutate(Abundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
    ungroup()
  


  
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
  
  sample_order_df <- df_plot %>%
    distinct(SampleID, SampleLabel, RegionPlot, TreatmentPlot, YearPlot, YearNum) %>%
    mutate(
      SampleLabel_chr = as.character(SampleLabel),
      OrderFixed = match(SampleLabel_chr, sample_order_fixed)
    )
  
  missing_order <- sample_order_df %>%
    filter(is.na(OrderFixed)) %>%
    pull(SampleLabel_chr) %>%
    unique()
  
  if (length(missing_order) > 0) {
    cat("\nNotice: these samples are not in the specified order and will be appended:\n")
    print(missing_order)
  }
  
  sample_order_df <- sample_order_df %>%
    mutate(
      OrderFixed2 = ifelse(
        is.na(OrderFixed),
        max(sample_order_fixed %in% sample_order_fixed) + row_number() + 1000,
        OrderFixed
      )
    ) %>%
    arrange(OrderFixed2)
  
  sample_levels <- sample_order_df$SampleID
  
  df_plot <- df_plot %>%
    filter(SampleID %in% sample_levels)
  
  df_plot$SampleID <- factor(df_plot$SampleID, levels = sample_levels)
  df_plot$SampleLabelFixed <- factor(
    as.character(df_plot$SampleLabel),
    levels = sample_order_df$SampleLabel_chr
  )
  
  df_plot$RegionPlot <- factor(df_plot$RegionPlot, levels = c("Amazonas", "Cusco"))
  df_plot$TreatmentPlot <- factor(df_plot$TreatmentPlot, levels = c("Control", "Wildfire"))
  


  
  genus_levels <- df_plot %>%
    group_by(GenusPlot) %>%
    summarise(Total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(Total)) %>%
    pull(GenusPlot)
  
  genus_levels <- c(setdiff(genus_levels, "Other"), "Other")
  
  df_plot$GenusPlot <- factor(df_plot$GenusPlot, levels = genus_levels)
  


  
  pal <- make_palette(length(genus_levels))
  names(pal) <- genus_levels
  
  if ("Other" %in% names(pal)) {
    pal["Other"] <- "grey70"
  }
  


  
  p <- ggplot(
    df_plot,
    aes(x = SampleLabelFixed, y = Abundance, fill = GenusPlot)
  ) +
    geom_col(
      width = 0.95,
      color = "black",
      linewidth = 0.05
    ) +
    facet_grid(
      ~ RegionPlot,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = pal,
      name = "Genus"
    ) +
    labs(
      title = paste0(marker_label, " - Taxonomic composition at genus level"),
      subtitle = "Identified genera only; 'Other' groups low-abundance identified genera",
      x = NULL,
      y = "Relative abundance"
    ) +
    theme_bw(base_size = 12) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 17,
        color = "black"
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 11,
        color = "black"
      ),
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 6,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 9,
        color = "black"
      ),
      axis.title.y = element_text(
        face = "bold",
        color = "black"
      ),
      strip.background = element_rect(
        fill = "grey90",
        color = "black"
      ),
      strip.text = element_text(
        face = "bold",
        size = 11,
        color = "black"
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.25),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(size = 8, color = "black"),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(p)
  


  
  outdir <- file.path(base, "results/plots", out_folder_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    df_plot,
    file.path(outdir, paste0(file_prefix, "_genus_relative_abundance_table.csv")),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(Genus = top_genera),
    file.path(outdir, paste0(file_prefix, "_top_", top_n, "_identified_genera.csv")),
    row.names = FALSE
  )
  
  write.csv(
    tax_genus[!is.na(tax_genus$Genus_clean), ],
    file.path(outdir, paste0(file_prefix, "_ASVs_with_identified_genus_used.csv")),
    row.names = FALSE
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_taxonomic_composition_Genus_identified_only.png")),
    p,
    width = 18,
    height = 8,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_taxonomic_composition_Genus_identified_only.pdf")),
    p,
    width = 18,
    height = 8,
    bg = "white"
  )
  
  cat("\nArchivos guardados en:\n")
  cat(outdir, "\n")
  
  cat("\nFigura PNG:\n")
  cat(file.path(outdir, paste0(file_prefix, "_taxonomic_composition_Genus_identified_only.png")), "\n")
  
  cat("\nFigura PDF:\n")
  cat(file.path(outdir, paste0(file_prefix, "_taxonomic_composition_Genus_identified_only.pdf")), "\n")
  
  return(p)
}


p_bacteria_genus <- plot_genus_composition_marker(
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
  out_folder_name = "taxonomy_composition_BAC_genus_identified_only",
  file_prefix = "BAC",
  top_n = 30
)


p_archaea_genus <- plot_genus_composition_marker(
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
  out_folder_name = "taxonomy_composition_ARQ_genus_identified_only",
  file_prefix = "ARQ",
  top_n = 30
)


p_its_genus <- plot_genus_composition_marker(
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
  out_folder_name = "taxonomy_composition_ITS_genus_identified_only",
  file_prefix = "ITS",
  top_n = 30
)

cat("\n=====================================================\n")
cat("Genus-level taxonomic plots generated\n")
cat("Criterion: identified genera only\n")
cat("Orden usado: orden fijo Amazonas + Cusco\n")
cat("Hora fin:", as.character(Sys.time()), "\n")
cat("=====================================================\n")


library(phyloseq)
library(dplyr)

base <- "/home/usuario1/Rosmery/Tesis Rosmery-Metabarcoding/RESULT ARQUEA"

ps <- readRDS(
  file.path(base, "data/processed_data/phyloseq/ps_sin_prevalencia_ARQ.rds")
)

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

tax_df <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
tax_df$ASV_ID <- rownames(tax_df)
tax_df$Genus_clean <- clean_genus(tax_df$Genus)

valid_taxa <- tax_df$ASV_ID[!is.na(tax_df$Genus_clean)]

ps_genus_tmp <- prune_taxa(valid_taxa, ps)

revision_cpat2 <- data.frame(
  SampleID = sample_names(ps),
  SampleLabel = clean_sample_label(sample_names(ps)),
  Reads_total_ARQ = sample_sums(ps),
  Reads_con_genero_identificado = sample_sums(ps_genus_tmp)[sample_names(ps)],
  stringsAsFactors = FALSE
)

revision_cpat2$Reads_con_genero_identificado[
  is.na(revision_cpat2$Reads_con_genero_identificado)
] <- 0

revision_cpat2 %>%
  filter(SampleLabel == "CPAT2")


rm(list = ls())


packages_cran <- c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "scales")

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

cat("\n=====================================================\n")
cat("ARCHAEAL TAXONOMIC COMPOSITION - GENUS LEVEL\n")
cat("SAMPLES WITHOUT GENUS ASSIGNMENTS = No identified genus\n")
cat("=====================================================\n")
cat("Hora inicio:", as.character(Sys.time()), "\n")


base_candidates <- path.expand(c(
  "~/Rosmery/Tesis Rosmery-Metabarcoding/RESULT ARQUEA",
  "~/Tesis Rosmery-Metabarcoding/RESULT ARQUEA"
))

base <- base_candidates[dir.exists(base_candidates)][1]

if (is.na(base)) {
  stop("Archaeal data directory not found.")
}

setwd(base)

cat("\nCarpeta usada:\n")
cat(base, "\n")


ps_file_candidates <- file.path(
  base,
  "data/processed_data/phyloseq",
  c(
    "ps_sin_prevalencia_ARQ.rds",
    "ps_only_archaea_ARQ.rds",
    "ps_raw_ARQ.rds"
  )
)

ps_file <- ps_file_candidates[file.exists(ps_file_candidates)][1]

if (is.na(ps_file)) {
  stop("Archaeal phyloseq object not found.")
}

cat("\nObjeto phyloseq usado:\n")
cat(ps_file, "\n")

ps <- readRDS(ps_file)

cat("\nObjeto cargado:\n")
cat("Muestras:", nsamples(ps), "\n")
cat("ASVs:", ntaxa(ps), "\n")
cat("Lecturas totales:", sum(sample_sums(ps)), "\n")


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

make_palette <- function(n) {
  
  base_cols <- c(
    brewer.pal(12, "Set3"),
    brewer.pal(8, "Set2"),
    brewer.pal(8, "Dark2"),
    brewer.pal(12, "Paired"),
    brewer.pal(9, "Set1")
  )
  
  base_cols <- unique(base_cols)
  
  if (n <= length(base_cols)) {
    return(base_cols[seq_len(n)])
  } else {
    return(colorRampPalette(base_cols)(n))
  }
}


meta_df <- data.frame(sample_data(ps), stringsAsFactors = FALSE, check.names = FALSE)
rownames(meta_df) <- sample_names(ps)

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

if (is.na(region_col)) stop("Region column not found.")
if (is.na(treat_col))  stop("Treatment column not found.")
if (is.na(year_col))   stop("Years-since-fire column not found.")

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

cat("\nSamples by region and treatment:\n")
print(table(meta_df$RegionPlot, meta_df$TreatmentPlot))


if (!"Genus" %in% rank_names(ps)) {
  stop("Taxonomy table lacks a Genus column.")
}

tax_df <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
tax_df$ASV_ID <- rownames(tax_df)
tax_df$Genus_clean <- clean_genus(tax_df$Genus)

valid_taxa <- tax_df$ASV_ID[!is.na(tax_df$Genus_clean)]

if (length(valid_taxa) == 0) {
  stop("No identified archaeal genera were found.")
}


ps_genus <- prune_taxa(valid_taxa, ps)

tax_genus <- as.data.frame(tax_table(ps_genus), stringsAsFactors = FALSE)
tax_genus$ASV_ID <- rownames(tax_genus)
tax_genus$Genus_clean <- clean_genus(tax_genus$Genus)

cat("\n===== ARCHAEAL TAXONOMY SUMMARY =====\n")
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

sample_totals <- colSums(genus_mat, na.rm = TRUE)

genus_rel <- sweep(genus_mat, 2, sample_totals, "/")
genus_rel[!is.finite(genus_rel)] <- 0


df_long <- as.data.frame(as.table(genus_rel), stringsAsFactors = FALSE)
colnames(df_long) <- c("Genus", "SampleID", "Abundance")

df_long$Abundance <- as.numeric(df_long$Abundance)

meta_plot <- meta_df[colnames(genus_rel), , drop = FALSE]
meta_plot$SampleID <- rownames(meta_plot)

df_long <- df_long %>%
  left_join(
    meta_plot %>%
      select(SampleID, SampleLabel, RegionPlot, TreatmentPlot, YearPlot, YearNum),
    by = "SampleID"
  )


top_n <- 30

top_genera <- df_long %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = top_n) %>%
  pull(Genus)

df_long$GenusPlot <- ifelse(df_long$Genus %in% top_genera, df_long$Genus, "Other")

df_plot <- df_long %>%
  group_by(SampleID, SampleLabel, RegionPlot, TreatmentPlot, YearPlot, YearNum, GenusPlot) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")


sample_sums_plot <- df_plot %>%
  group_by(SampleID) %>%
  summarise(SumAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

zero_samples <- sample_sums_plot$SampleID[
  sample_sums_plot$SumAbundance == 0 | is.na(sample_sums_plot$SumAbundance)
]

df_plot <- df_plot %>%
  left_join(sample_sums_plot, by = "SampleID") %>%
  mutate(
    Abundance = ifelse(
      SumAbundance > 0,
      Abundance / SumAbundance,
      0
    )
  ) %>%
  select(-SumAbundance)

if (length(zero_samples) > 0) {
  
  cat("\nSamples without an identified genus displayed as an empty bar:\n")
  print(clean_sample_label(zero_samples))
  
  zero_meta <- meta_plot[zero_samples, , drop = FALSE]
  
  zero_rows <- data.frame(
    SampleID = zero_meta$SampleID,
    SampleLabel = zero_meta$SampleLabel,
    RegionPlot = zero_meta$RegionPlot,
    TreatmentPlot = zero_meta$TreatmentPlot,
    YearPlot = zero_meta$YearPlot,
    YearNum = zero_meta$YearNum,
    GenusPlot = "No identified genus",
    Abundance = 1,
    stringsAsFactors = FALSE
  )
  
  df_plot <- df_plot %>%
    filter(!SampleID %in% zero_samples) %>%
    bind_rows(zero_rows)
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

sample_order_df <- df_plot %>%
  distinct(SampleID, SampleLabel, RegionPlot, TreatmentPlot, YearPlot, YearNum) %>%
  mutate(
    SampleLabel_chr = as.character(SampleLabel),
    OrderFixed = match(SampleLabel_chr, sample_order_fixed)
  )

missing_order <- sample_order_df %>%
  filter(is.na(OrderFixed)) %>%
  pull(SampleLabel_chr) %>%
  unique()

if (length(missing_order) > 0) {
  cat("\nNotice: these samples are not in the specified order and will be appended:\n")
  print(missing_order)
}

sample_order_df <- sample_order_df %>%
  mutate(
    OrderFixed2 = ifelse(
      is.na(OrderFixed),
      length(sample_order_fixed) + row_number(),
      OrderFixed
    )
  ) %>%
  arrange(OrderFixed2)

cat("\nORDEN FINAL DE MUESTRAS EN LA FIGURA:\n")
print(sample_order_df$SampleLabel_chr)

sample_levels <- sample_order_df$SampleID

df_plot <- df_plot %>%
  filter(SampleID %in% sample_levels)

df_plot$SampleID <- factor(df_plot$SampleID, levels = sample_levels)

df_plot$SampleLabelFixed <- factor(
  as.character(df_plot$SampleLabel),
  levels = sample_order_df$SampleLabel_chr
)

df_plot$RegionPlot <- factor(df_plot$RegionPlot, levels = c("Amazonas", "Cusco"))
df_plot$TreatmentPlot <- factor(df_plot$TreatmentPlot, levels = c("Control", "Wildfire"))


genus_levels <- df_plot %>%
  group_by(GenusPlot) %>%
  summarise(Total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(GenusPlot)

genus_levels <- c(
  setdiff(genus_levels, c("Other", "No identified genus")),
  "Other",
  "No identified genus"
)

genus_levels <- genus_levels[genus_levels %in% unique(df_plot$GenusPlot)]

df_plot$GenusPlot <- factor(df_plot$GenusPlot, levels = genus_levels)

pal <- make_palette(length(genus_levels))
names(pal) <- genus_levels

if ("Other" %in% names(pal)) {
  pal["Other"] <- "grey70"
}

if ("No identified genus" %in% names(pal)) {
  pal["No identified genus"] <- "white"
}


p_arq_genus <- ggplot(
  df_plot,
  aes(x = SampleLabelFixed, y = Abundance, fill = GenusPlot)
) +
  geom_col(
    width = 0.95,
    color = "black",
    linewidth = 0.12
  ) +
  facet_grid(
    ~ RegionPlot,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = pal,
    name = "Genus",
    drop = FALSE
  ) +
  labs(
    title = "Archaeal 16S taxonomic composition at genus level",
    subtitle = "Identified genera; samples without genus assignments are displayed as 'No identified genus'",
    x = NULL,
    y = "Relative abundance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 17,
      color = "black"
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11,
      color = "black"
    ),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 6,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    axis.title.y = element_text(
      face = "bold",
      color = "black"
    ),
    strip.background = element_rect(
      fill = "grey90",
      color = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11,
      color = "black"
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.25),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(p_arq_genus)


outdir <- file.path(base, "results/plots/taxonomy_composition_ARQ_genus_identified_only")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  df_plot,
  file.path(outdir, "ARQ_genus_relative_abundance_table.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(Genus = top_genera),
  file.path(outdir, paste0("ARQ_top_", top_n, "_identified_genera.csv")),
  row.names = FALSE
)

write.csv(
  tax_genus[!is.na(tax_genus$Genus_clean), ],
  file.path(outdir, "ARQ_ASVs_with_identified_genus_used.csv"),
  row.names = FALSE
)

ggsave(
  file.path(outdir, "ARQ_taxonomic_composition_Genus_identified_only.png"),
  p_arq_genus,
  width = 18,
  height = 8,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(outdir, "ARQ_taxonomic_composition_Genus_identified_only.pdf"),
  p_arq_genus,
  width = 18,
  height = 8,
  bg = "white"
)

cat("\nArchivos guardados en:\n")
cat(outdir, "\n")

cat("\nFigura PNG:\n")
cat(file.path(outdir, "ARQ_taxonomic_composition_Genus_identified_only.png"), "\n")

cat("\nFigura PDF:\n")
cat(file.path(outdir, "ARQ_taxonomic_composition_Genus_identified_only.pdf"), "\n")

cat("\n=====================================================\n")
cat("Archaeal composition plot generated\n")
cat("CPAT2 will be shown as 'No identified genus' if no genus has been assigned.\n")
cat("Hora fin:", as.character(Sys.time()), "\n")
cat("=====================================================\n")


rm(list = ls())


packages_cran <- c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "scales")

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

cat("\n=====================================================\n")
cat("TAXONOMIC COMPOSITION - TOP 30 TAXA\n")
cat("NIVELES: PHYLUM, CLASS, ORDER, FAMILY\n")
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

clean_taxon <- function(x) {
  x <- trimws(as.character(x))
  

  x <- gsub("^k__", "", x)
  x <- gsub("^p__", "", x)
  x <- gsub("^c__", "", x)
  x <- gsub("^o__", "", x)
  x <- gsub("^f__", "", x)
  x <- gsub("^g__", "", x)
  x <- gsub("^s__", "", x)
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

make_palette <- function(n) {
  
  base_cols <- c(
    brewer.pal(12, "Set3"),
    brewer.pal(8, "Set2"),
    brewer.pal(8, "Dark2"),
    brewer.pal(12, "Paired"),
    brewer.pal(9, "Set1"),
    brewer.pal(8, "Accent")
  )
  
  base_cols <- unique(base_cols)
  
  if (n <= length(base_cols)) {
    return(base_cols[seq_len(n)])
  } else {
    return(colorRampPalette(base_cols)(n))
  }
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


plot_taxonomic_level_marker <- function(marker_name,
                                        marker_label,
                                        base_candidates,
                                        ps_candidates,
                                        file_prefix,
                                        tax_level,
                                        top_n = 30) {
  
  cat("\n=====================================================\n")
  cat("PROCESANDO:", marker_label, "-", tax_level, "\n")
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
  
  cat("\nSamples by region and treatment:\n")
  print(table(meta_df$RegionPlot, meta_df$TreatmentPlot))
  


  
  if (!tax_level %in% rank_names(ps)) {
    stop("Taxonomic rank not available: ", tax_level, " en: ", marker_label)
  }
  
  tax_df <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  tax_df$ASV_ID <- rownames(tax_df)
  tax_df$Taxon_clean <- clean_taxon(tax_df[[tax_level]])
  
  valid_taxa <- tax_df$ASV_ID[!is.na(tax_df$Taxon_clean)]
  
  if (length(valid_taxa) == 0) {
    stop("No valid taxonomic assignments at ", tax_level, " para ", marker_label)
  }
  

  ps_tax <- prune_taxa(valid_taxa, ps)
  
  tax_used <- as.data.frame(tax_table(ps_tax), stringsAsFactors = FALSE)
  tax_used$ASV_ID <- rownames(tax_used)
  tax_used$Taxon_clean <- clean_taxon(tax_used[[tax_level]])
  
  cat("\n===== TAXONOMY SUMMARY =====\n")
  cat("Nivel:", tax_level, "\n")
  cat("ASVs totales:", ntaxa(ps), "\n")
  cat("ASVs con", tax_level, "identificado:", ntaxa(ps_tax), "\n")
  cat(tax_level, "unique identified taxa:", length(unique(na.omit(tax_used$Taxon_clean))), "\n")
  cat("Muestras conservadas:", nsamples(ps_tax), "\n")
  


  
  otu_mat <- as(otu_table(ps_tax), "matrix")
  
  if (!taxa_are_rows(ps_tax)) {
    otu_mat <- t(otu_mat)
  }
  
  tax_used <- tax_used[rownames(otu_mat), , drop = FALSE]
  
  taxon_vector <- tax_used$Taxon_clean
  names(taxon_vector) <- rownames(tax_used)
  
  valid_rows <- !is.na(taxon_vector)
  
  otu_valid <- otu_mat[valid_rows, , drop = FALSE]
  taxon_valid <- taxon_vector[valid_rows]
  
  taxon_mat <- rowsum(
    otu_valid,
    group = taxon_valid,
    reorder = TRUE
  )
  
  sample_totals <- colSums(taxon_mat, na.rm = TRUE)
  
  taxon_rel <- sweep(taxon_mat, 2, sample_totals, "/")
  taxon_rel[!is.finite(taxon_rel)] <- 0
  


  
  df_long <- as.data.frame(as.table(taxon_rel), stringsAsFactors = FALSE)
  colnames(df_long) <- c("Taxon", "SampleID", "Abundance")
  
  df_long$Abundance <- as.numeric(df_long$Abundance)
  
  meta_plot <- meta_df[colnames(taxon_rel), , drop = FALSE]
  meta_plot$SampleID <- rownames(meta_plot)
  
  df_long <- df_long %>%
    left_join(
      meta_plot %>%
        select(SampleID, SampleLabel, RegionPlot, TreatmentPlot, YearPlot, YearNum),
      by = "SampleID"
    )
  


  
  top_taxa <- df_long %>%
    group_by(Taxon) %>%
    summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalAbundance)) %>%
    slice_head(n = top_n) %>%
    pull(Taxon)
  
  df_long$TaxonPlot <- ifelse(df_long$Taxon %in% top_taxa, df_long$Taxon, "Other")
  
  df_plot <- df_long %>%
    group_by(SampleID, SampleLabel, RegionPlot, TreatmentPlot, YearPlot, YearNum, TaxonPlot) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  


  
  sample_sums_plot <- df_plot %>%
    group_by(SampleID) %>%
    summarise(SumAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  zero_samples <- sample_sums_plot$SampleID[
    sample_sums_plot$SumAbundance == 0 | is.na(sample_sums_plot$SumAbundance)
  ]
  
  df_plot <- df_plot %>%
    left_join(sample_sums_plot, by = "SampleID") %>%
    mutate(
      Abundance = ifelse(
        SumAbundance > 0,
        Abundance / SumAbundance,
        0
      )
    ) %>%
    select(-SumAbundance)
  
  no_id_label <- paste0("No identified ", tax_level)
  
  if (length(zero_samples) > 0) {
    
    cat("\nSamples unassigned at", tax_level, "shown as an empty bar:\n")
    print(clean_sample_label(zero_samples))
    
    zero_meta <- meta_plot[zero_samples, , drop = FALSE]
    
    zero_rows <- data.frame(
      SampleID = zero_meta$SampleID,
      SampleLabel = zero_meta$SampleLabel,
      RegionPlot = zero_meta$RegionPlot,
      TreatmentPlot = zero_meta$TreatmentPlot,
      YearPlot = zero_meta$YearPlot,
      YearNum = zero_meta$YearNum,
      TaxonPlot = no_id_label,
      Abundance = 1,
      stringsAsFactors = FALSE
    )
    
    df_plot <- df_plot %>%
      filter(!SampleID %in% zero_samples) %>%
      bind_rows(zero_rows)
  }
  


  
  sample_order_df <- df_plot %>%
    distinct(SampleID, SampleLabel, RegionPlot, TreatmentPlot, YearPlot, YearNum) %>%
    mutate(
      SampleLabel_chr = as.character(SampleLabel),
      OrderFixed = match(SampleLabel_chr, sample_order_fixed)
    )
  
  missing_order <- sample_order_df %>%
    filter(is.na(OrderFixed)) %>%
    pull(SampleLabel_chr) %>%
    unique()
  
  if (length(missing_order) > 0) {
    cat("\nNotice: these samples are not in the specified order and will be appended:\n")
    print(missing_order)
  }
  
  sample_order_df <- sample_order_df %>%
    mutate(
      OrderFixed2 = ifelse(
        is.na(OrderFixed),
        length(sample_order_fixed) + row_number(),
        OrderFixed
      )
    ) %>%
    arrange(OrderFixed2)
  
  cat("\nORDEN FINAL DE MUESTRAS EN LA FIGURA:\n")
  print(sample_order_df$SampleLabel_chr)
  
  sample_levels <- sample_order_df$SampleID
  
  df_plot <- df_plot %>%
    filter(SampleID %in% sample_levels)
  
  df_plot$SampleID <- factor(df_plot$SampleID, levels = sample_levels)
  
  df_plot$SampleLabelFixed <- factor(
    as.character(df_plot$SampleLabel),
    levels = sample_order_df$SampleLabel_chr
  )
  
  df_plot$RegionPlot <- factor(df_plot$RegionPlot, levels = c("Amazonas", "Cusco"))
  df_plot$TreatmentPlot <- factor(df_plot$TreatmentPlot, levels = c("Control", "Wildfire"))
  


  
  taxon_levels <- df_plot %>%
    group_by(TaxonPlot) %>%
    summarise(Total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(Total)) %>%
    pull(TaxonPlot)
  
  taxon_levels <- c(
    setdiff(taxon_levels, c("Other", no_id_label)),
    "Other",
    no_id_label
  )
  
  taxon_levels <- taxon_levels[taxon_levels %in% unique(df_plot$TaxonPlot)]
  
  df_plot$TaxonPlot <- factor(df_plot$TaxonPlot, levels = taxon_levels)
  
  pal <- make_palette(length(taxon_levels))
  names(pal) <- taxon_levels
  
  if ("Other" %in% names(pal)) {
    pal["Other"] <- "grey70"
  }
  
  if (no_id_label %in% names(pal)) {
    pal[no_id_label] <- "white"
  }
  


  
  p <- ggplot(
    df_plot,
    aes(x = SampleLabelFixed, y = Abundance, fill = TaxonPlot)
  ) +
    geom_col(
      width = 0.95,
      color = "black",
      linewidth = 0.12
    ) +
    facet_grid(
      ~ RegionPlot,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = pal,
      name = tax_level,
      drop = FALSE
    ) +
    labs(
      title = paste0(marker_label, " - Taxonomic composition at ", tax_level),
      subtitle = paste0("Top ", top_n, " ", tax_level, "; 'Other' agrupa taxones identificados de menor abundancia"),
      x = NULL,
      y = "Relative abundance"
    ) +
    theme_bw(base_size = 12) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 17,
        color = "black"
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 11,
        color = "black"
      ),
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 6,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 9,
        color = "black"
      ),
      axis.title.y = element_text(
        face = "bold",
        color = "black"
      ),
      strip.background = element_rect(
        fill = "grey90",
        color = "black"
      ),
      strip.text = element_text(
        face = "bold",
        size = 11,
        color = "black"
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.25),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(size = 8, color = "black"),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(p)
  


  
  out_folder_name <- paste0(
    "taxonomy_composition_",
    file_prefix,
    "_",
    tax_level,
    "_top30_identified_only"
  )
  
  outdir <- file.path(base, "results/plots", out_folder_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    df_plot,
    file.path(outdir, paste0(file_prefix, "_", tax_level, "_relative_abundance_table_top30.csv")),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(Taxon = top_taxa),
    file.path(outdir, paste0(file_prefix, "_top30_", tax_level, "_identified.csv")),
    row.names = FALSE
  )
  
  write.csv(
    tax_used[!is.na(tax_used$Taxon_clean), ],
    file.path(outdir, paste0(file_prefix, "_ASVs_with_identified_", tax_level, "_used.csv")),
    row.names = FALSE
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_taxonomic_composition_", tax_level, "_top30_identified_only.png")),
    p,
    width = 18,
    height = 8,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_taxonomic_composition_", tax_level, "_top30_identified_only.pdf")),
    p,
    width = 18,
    height = 8,
    bg = "white"
  )
  
  cat("\nArchivos guardados en:\n")
  cat(outdir, "\n")
  
  cat("\nFigura PNG:\n")
  cat(file.path(outdir, paste0(file_prefix, "_taxonomic_composition_", tax_level, "_top30_identified_only.png")), "\n")
  
  cat("\nFigura PDF:\n")
  cat(file.path(outdir, paste0(file_prefix, "_taxonomic_composition_", tax_level, "_top30_identified_only.pdf")), "\n")
  
  return(p)
}


tax_levels_to_plot <- c("Phylum", "Class", "Order", "Family")
top_n <- 30

markers <- list(
  
  Bacteria = list(
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
  ),
  
  Archaea = list(
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
  ),
  
  ITS = list(
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
)


plots <- list()

for (marker_id in names(markers)) {
  
  mk <- markers[[marker_id]]
  
  for (tax_level in tax_levels_to_plot) {
    
    plot_name <- paste(marker_id, tax_level, sep = "_")
    
    plots[[plot_name]] <- plot_taxonomic_level_marker(
      marker_name = mk$marker_name,
      marker_label = mk$marker_label,
      base_candidates = mk$base_candidates,
      ps_candidates = mk$ps_candidates,
      file_prefix = mk$file_prefix,
      tax_level = tax_level,
      top_n = top_n
    )
  }
}

cat("\n=====================================================\n")
cat("Top-30 taxonomic plots generated\n")
cat("Niveles: Phylum, Class, Order, Family\n")
cat("Marcadores: Bacterias, Arqueas, ITS\n")
cat("Orden usado: orden fijo Amazonas + Cusco\n")
cat("Hora fin:", as.character(Sys.time()), "\n")
cat("=====================================================\n")
