# Shared and exclusive genera
# Genus-level presence comparisons by region, treatment and years since fire.
# Input: project-specific phyloseq objects, site metadata, or soil data tables.
# Output: analysis tables and figures as defined by the project directory settings.


rm(list = ls())


if (!requireNamespace("phyloseq", quietly = TRUE)) {
  stop("The phyloseq package is required to run this script.")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork", repos = "https://cloud.r-project.org")
}

library(phyloseq)
library(ggplot2)
library(patchwork)

cat("\n=====================================================\n")
cat("VENN DIAGRAMS: CONTROL VS BURNED - IDENTIFIED GENERA\n")
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
  out <- ifelse(
    grepl("amazon", x, ignore.case = TRUE),
    "Amazonas",
    ifelse(
      grepl("cusco|cuzco", x, ignore.case = TRUE),
      "Cusco",
      x
    )
  )
  return(out)
}

normalize_treatment <- function(x) {
  x <- trimws(as.character(x))
  xl <- tolower(x)
  
  out <- ifelse(
    xl %in% c("control", "ctrl", "unburned", "unburnt", "no burned", "no burn"),
    "Control",
    ifelse(
      xl %in% c("burned", "burnt", "wildfire", "fire", "incendio", "incendiado", "quemado"),
      "Burned",
      x
    )
  )
  
  return(out)
}

clean_genus <- function(x) {
  x <- trimws(as.character(x))
  

  x <- gsub("^g__", "", x)
  

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

make_ellipse <- function(group_name, x0, y0, a = 2.55, b = 1.35, n = 600) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(
    x = x0 + a * cos(t),
    y = y0 + b * sin(t),
    Group = group_name
  )
}

make_venn_plot <- function(set_burned,
                           set_control,
                           region_title,
                           fill_burned,
                           fill_control) {
  
  only_burned <- setdiff(set_burned, set_control)
  shared <- intersect(set_burned, set_control)
  only_control <- setdiff(set_control, set_burned)
  
  total <- length(union(set_burned, set_control))
  
  if (total == 0) {
    return(
      ggplot() +
        theme_void() +
        annotate(
          "text",
          x = 0,
          y = 0,
          label = paste0(region_title, "\nNo identified genera"),
          fontface = "bold",
          size = 5
        )
    )
  }
  
  pct <- function(x) round(x / total * 100, 1)
  
  ellipse_df <- rbind(
    make_ellipse("Burned",  x0 = 0, y0 = 0.65),
    make_ellipse("Control", x0 = 0, y0 = -0.65)
  )
  
  label_burned <- paste0(
    length(only_burned), "\n(",
    pct(length(only_burned)), "%)"
  )
  
  label_shared <- paste0(
    length(shared), "\n(",
    pct(length(shared)), "%)"
  )
  
  label_control <- paste0(
    length(only_control), "\n(",
    pct(length(only_control)), "%)"
  )
  
  p <- ggplot() +
    geom_polygon(
      data = ellipse_df,
      aes(x = x, y = y, group = Group, fill = Group),
      alpha = 0.72,
      color = "black",
      linewidth = 0.75
    ) +
    annotate(
      "text",
      x = 0,
      y = 1.35,
      label = label_burned,
      fontface = "bold",
      size = 4.8,
      color = "black"
    ) +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = label_shared,
      fontface = "bold",
      size = 4.8,
      color = "black"
    ) +
    annotate(
      "text",
      x = 0,
      y = -1.35,
      label = label_control,
      fontface = "bold",
      size = 4.8,
      color = "black"
    ) +
    scale_fill_manual(
      name = NULL,
      values = c(
        "Burned" = fill_burned,
        "Control" = fill_control
      ),
      breaks = c("Burned", "Control")
    ) +
    coord_fixed(
      xlim = c(-2.85, 2.85),
      ylim = c(-2.05, 2.05),
      expand = FALSE
    ) +
    ggtitle(region_title) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16,
        color = "black"
      ),
      legend.position = "left",
      legend.text = element_text(
        face = "bold",
        size = 11,
        color = "black"
      ),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      text = element_text(color = "black")
    )
  
  return(p)
}


run_venn_marker <- function(marker_name,
                            marker_label,
                            base_candidates,
                            ps_candidates,
                            out_folder_name,
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
  
  if (!inherits(ps, "phyloseq")) {
    stop("Input file is not a valid phyloseq object.")
  }
  
  cat("\nObjeto cargado:\n")
  cat("Muestras:", nsamples(ps), "\n")
  cat("ASVs:", ntaxa(ps), "\n")
  cat("Lecturas totales:", sum(sample_sums(ps)), "\n")
  


  
  sd <- data.frame(sample_data(ps), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(sd) <- sample_names(ps)
  
  region_col <- find_col(
    colnames(sd),
    c("Region", "region", "Departamento", "departamento")
  )
  
  treat_col <- find_col(
    colnames(sd),
    c("Treatment", "treatment", "Tratamiento", "tratamiento", "Tratament")
  )
  
  if (is.na(region_col)) {
    stop("Region column not found for: ", marker_label)
  }
  
  if (is.na(treat_col)) {
    stop("Treatment column not found for: ", marker_label)
  }
  
  sd[[region_col]] <- normalize_region(sd[[region_col]])
  sd[[treat_col]]  <- normalize_treatment(sd[[treat_col]])
  
  sample_data(ps) <- sample_data(sd)
  
  cat("\nRegion column:", region_col, "\n")
  cat("Columna tratamiento:", treat_col, "\n")
  
  cat("\nSamples by region and treatment:\n")
  print(table(sd[[region_col]], sd[[treat_col]], useNA = "ifany"))
  

  keep_samples <- rownames(sd)[
    sd[[region_col]] %in% c("Amazonas", "Cusco") &
      sd[[treat_col]] %in% c("Control", "Burned")
  ]
  
  ps <- prune_samples(keep_samples, ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  


  
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
  
  sd_genus <- data.frame(sample_data(ps_genus), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(sd_genus) <- sample_names(ps_genus)
  sd_genus$SampleID <- rownames(sd_genus)
  
  get_genera_set <- function(region_value, treatment_value) {
    
    samples_use <- sd_genus$SampleID[
      sd_genus[[region_col]] == region_value &
        sd_genus[[treat_col]] == treatment_value
    ]
    
    samples_use <- intersect(samples_use, colnames(otu_mat))
    
    if (length(samples_use) == 0) {
      return(character(0))
    }
    
    taxa_present <- rownames(otu_mat)[
      rowSums(otu_mat[, samples_use, drop = FALSE]) > 0
    ]
    
    genera_present <- tax_genus$Genus_clean[
      match(taxa_present, tax_genus$ASV_ID)
    ]
    
    genera_present <- sort(unique(na.omit(genera_present)))
    
    return(genera_present)
  }
  


  
  cusco_burned   <- get_genera_set("Cusco", "Burned")
  cusco_control  <- get_genera_set("Cusco", "Control")
  
  amazon_burned  <- get_genera_set("Amazonas", "Burned")
  amazon_control <- get_genera_set("Amazonas", "Control")
  
  summarize_region <- function(region_name, burned_set, control_set) {
    only_burned <- setdiff(burned_set, control_set)
    shared <- intersect(burned_set, control_set)
    only_control <- setdiff(control_set, burned_set)
    total <- length(union(burned_set, control_set))
    
    data.frame(
      Marker = marker_name,
      Region = region_name,
      Only_Burned = length(only_burned),
      Shared = length(shared),
      Only_Control = length(only_control),
      Total_identified_genera = total,
      Percent_Only_Burned = ifelse(total > 0, round(length(only_burned) / total * 100, 1), 0),
      Percent_Shared = ifelse(total > 0, round(length(shared) / total * 100, 1), 0),
      Percent_Only_Control = ifelse(total > 0, round(length(only_control) / total * 100, 1), 0),
      stringsAsFactors = FALSE
    )
  }
  
  summary_table <- rbind(
    summarize_region("Cusco", cusco_burned, cusco_control),
    summarize_region("Amazonas", amazon_burned, amazon_control)
  )
  
  print(summary_table)
  


  
  make_presence_table <- function(region_name, burned_set, control_set) {
    all_genera <- sort(unique(c(burned_set, control_set)))
    
    data.frame(
      Marker = marker_name,
      Region = region_name,
      Genus = all_genera,
      Burned = all_genera %in% burned_set,
      Control = all_genera %in% control_set,
      Status = ifelse(
        all_genera %in% intersect(burned_set, control_set),
        "Shared",
        ifelse(all_genera %in% setdiff(burned_set, control_set),
               "Only_Burned",
               "Only_Control")
      ),
      stringsAsFactors = FALSE
    )
  }
  
  presence_table <- rbind(
    make_presence_table("Cusco", cusco_burned, cusco_control),
    make_presence_table("Amazonas", amazon_burned, amazon_control)
  )
  


  
  p_cusco <- make_venn_plot(
    set_burned = cusco_burned,
    set_control = cusco_control,
    region_title = "CUSCO",
    fill_burned = "#6FA7CC",
    fill_control = "#CBEFF2"
  )
  
  p_amazonas <- make_venn_plot(
    set_burned = amazon_burned,
    set_control = amazon_control,
    region_title = "AMAZONAS",
    fill_burned = "#8FC5B8",
    fill_control = "#D2F4D4"
  )
  
  plot_final <- p_cusco + p_amazonas +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = paste0(marker_label, ": identified genera shared between Control and Burned"),
      theme = theme(
        plot.title = element_text(
          size = 18,
          face = "bold",
          hjust = 0.5,
          color = "black"
        ),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
  
  print(plot_final)
  


  
  outdir <- file.path(base, "results/plots", out_folder_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    summary_table,
    file.path(outdir, paste0(file_prefix, "_Venn_Control_vs_Burned_identified_genera_summary.csv")),
    row.names = FALSE
  )
  
  write.csv(
    presence_table,
    file.path(outdir, paste0(file_prefix, "_Venn_Control_vs_Burned_identified_genera_presence.csv")),
    row.names = FALSE
  )
  
  write.csv(
    tax_genus[!is.na(tax_genus$Genus_clean), ],
    file.path(outdir, paste0(file_prefix, "_ASVs_with_identified_genus_used.csv")),
    row.names = FALSE
  )
  
  ggsave(
    filename = file.path(outdir, paste0(file_prefix, "_Venn_Control_vs_Burned_identified_genera_Cusco_Amazonas.png")),
    plot = plot_final,
    width = 13,
    height = 6.5,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(outdir, paste0(file_prefix, "_Venn_Control_vs_Burned_identified_genera_Cusco_Amazonas.pdf")),
    plot = plot_final,
    width = 13,
    height = 6.5,
    bg = "white"
  )
  
  cat("\nArchivos guardados en:\n")
  cat(outdir, "\n")
  
  cat("\nFigura PNG:\n")
  cat(file.path(outdir, paste0(file_prefix, "_Venn_Control_vs_Burned_identified_genera_Cusco_Amazonas.png")), "\n")
  
  cat("\nFigura PDF:\n")
  cat(file.path(outdir, paste0(file_prefix, "_Venn_Control_vs_Burned_identified_genera_Cusco_Amazonas.pdf")), "\n")
  
  return(plot_final)
}


p_bacteria <- run_venn_marker(
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
  out_folder_name = "venn_control_burned_BAC_identified_genera",
  file_prefix = "BAC"
)


p_archaea <- run_venn_marker(
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
  out_folder_name = "venn_control_burned_ARQ_identified_genera",
  file_prefix = "ARQ"
)


p_its <- run_venn_marker(
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
  out_folder_name = "venn_control_burned_ITS_identified_genera",
  file_prefix = "ITS"
)

cat("\n=====================================================\n")
cat("TODOS LOS VENN FINALIZADOS\n")
cat("Se generaron figuras para: Bacterias, Arqueas e ITS\n")
cat("Criterion: identified genera only\n")
cat("Hora fin:", as.character(Sys.time()), "\n")
cat("=====================================================\n")


rm(list = ls())


packages_cran <- c("ggplot2", "dplyr", "ggVennDiagram", "cowplot")

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
library(ggVennDiagram)
library(cowplot)

cat("\n=====================================================\n")
cat("VENN DIAGRAMS BY YEARS SINCE FIRE - IDENTIFIED GENERA\n")
cat("PERCENTAGE SUMMARIES\n")
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
      "Burned",
      x
    )
  )
}

normalize_year <- function(x) {
  x <- trimws(as.character(x))
  xn <- suppressWarnings(as.numeric(x))
  ifelse(is.na(xn), x, as.character(as.integer(round(xn))))
}

clean_genus <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^g__", "", x)
  
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

move_labels_outward <- function(df, amount = 0.08) {
  if (nrow(df) == 0) return(df)
  
  cx <- mean(range(df$X, na.rm = TRUE))
  cy <- mean(range(df$Y, na.rm = TRUE))
  
  df %>%
    rowwise() %>%
    mutate(
      dx0 = X - cx,
      dy0 = Y - cy,
      rr = sqrt(dx0^2 + dy0^2),
      dx = ifelse(is.na(rr) | rr == 0, 0, dx0 / rr),
      dy = ifelse(is.na(rr) | rr == 0, 0, dy0 / rr),
      X = X + dx * amount,
      Y = Y + dy * amount
    ) %>%
    ungroup()
}

make_legend_plot <- function(years, palette_named) {
  
  df <- data.frame(
    Year = factor(years, levels = years),
    x = 1,
    y = seq(length(years), 1),
    stringsAsFactors = FALSE
  )
  
  ggplot(df, aes(x = x, y = y, fill = Year)) +
    geom_tile(width = 0.45, height = 0.45, color = "black", linewidth = 0.35) +
    geom_text(aes(label = Year), x = 1.38, hjust = 0, size = 4.2, color = "black") +
    scale_fill_manual(values = palette_named[years]) +
    coord_cartesian(
      xlim = c(0.7, 2.0),
      ylim = c(0.5, length(years) + 0.6),
      clip = "off"
    ) +
    theme_void() +
    ggtitle("Year") +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0, color = "black"),
      legend.position = "none",
      plot.margin = margin(8, 8, 8, 8)
    )
}


make_multiyear_venn_plot <- function(sets,
                                     title_text,
                                     palette_named,
                                     region_value) {
  

  sets <- sets[lengths(sets) > 0]
  
  if (length(sets) < 2) {
    return(
      ggplot() +
        theme_void() +
        ggtitle(title_text) +
        annotate(
          "text",
          x = 0,
          y = 0,
          label = "Insufficient\nyear groups with identified genera",
          size = 5,
          fontface = "bold"
        ) +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black")
        )
    )
  }
  

  total_taxa <- length(unique(unlist(sets)))
  
  venn_obj <- Venn(sets)
  vd <- process_data(venn_obj)
  
  df_set <- vd$setEdge %>%
    left_join(vd$setData, by = "id")
  
  df_set$name <- as.character(df_set$name)
  df_set$name <- factor(df_set$name, levels = names(sets))
  


  
  region_df <- vd$regionLabel
  

  if ("item" %in% colnames(region_df) && is.list(region_df$item)) {
    region_df$count_real <- lengths(region_df$item)
  } else {
    region_df$count_real <- suppressWarnings(as.numeric(as.character(region_df$count)))
  }
  
  region_df <- region_df %>%
    filter(!is.na(count_real), count_real > 0) %>%
    mutate(
      percent_real = round(count_real / total_taxa * 100, 1),
      label = paste0(count_real, "\n(", percent_real, "%)")
    )
  
  region_df <- move_labels_outward(
    region_df,
    amount = ifelse(region_value == "Cusco", 0.035, 0.025)
  )
  
  set_sizes <- data.frame(
    name = names(sets),
    nset = vapply(sets, length, integer(1)),
    stringsAsFactors = FALSE
  )
  
  set_df <- vd$setLabel %>%
    mutate(name = as.character(name)) %>%
    left_join(set_sizes, by = "name")
  
  set_df <- move_labels_outward(
    set_df,
    amount = ifelse(region_value == "Cusco", 0.18, 0.13)
  )
  
  set_df <- set_df %>%
    mutate(label_ext = paste0(nset, " genus (", name, "y)"))
  
  ggplot() +
    geom_polygon(
      data = df_set,
      aes(x = X, y = Y, group = id, fill = name),
      alpha = 0.48,
      color = "black",
      linewidth = 0.75
    ) +
    geom_text(
      data = region_df,
      aes(x = X, y = Y, label = label),
      size = ifelse(region_value == "Cusco", 3.0, 3.2),
      fontface = "bold",
      lineheight = 0.90,
      color = "black"
    ) +
    geom_text(
      data = set_df,
      aes(x = X, y = Y, label = label_ext),
      size = 3.7,
      color = "black"
    ) +
    scale_fill_manual(values = palette_named[names(sets)], drop = FALSE) +
    coord_equal(clip = "off") +
    ggtitle(title_text) +
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
      plot.margin = margin(t = 20, r = 35, b = 20, l = 35),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "none",
      text = element_text(color = "black")
    )
}


run_multiyear_venn_marker <- function(marker_name,
                                      marker_label,
                                      base_candidates,
                                      ps_candidates,
                                      out_folder_name,
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
  
  sd <- data.frame(sample_data(ps), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(sd) <- sample_names(ps)
  
  region_col <- find_col(
    colnames(sd),
    c("Region", "region", "Departamento", "departamento")
  )
  
  treat_col <- find_col(
    colnames(sd),
    c("Treatment", "treatment", "Tratamiento", "tratamiento", "Tratament")
  )
  
  year_col <- find_col(
    colnames(sd),
    c("Years_since_fire", "Years_since_fire.1", "YearsF", "Year", "year", "Año", "Anio")
  )
  
  if (is.na(region_col)) stop("Region column not found for: ", marker_label)
  if (is.na(treat_col))  stop("Treatment column not found for: ", marker_label)
  if (is.na(year_col))   stop("Years-since-fire column not found for: ", marker_label)
  
  sd[[region_col]] <- normalize_region(sd[[region_col]])
  sd[[treat_col]]  <- normalize_treatment(sd[[treat_col]])
  sd[[year_col]]   <- normalize_year(sd[[year_col]])
  
  sample_data(ps) <- sample_data(sd)
  
  cat("\nRegion column:", region_col, "\n")
  cat("Columna tratamiento:", treat_col, "\n")
  cat("Years-since-fire column:", year_col, "\n")
  
  cat("\nSamples by region and treatment:\n")
  print(table(sd[[region_col]], sd[[treat_col]], useNA = "ifany"))
  
  cat("\nSamples by region, treatment and years since fire:\n")
  print(table(sd[[region_col]], sd[[treat_col]], sd[[year_col]], useNA = "ifany"))
  
  keep_samples <- rownames(sd)[
    sd[[region_col]] %in% c("Amazonas", "Cusco") &
      sd[[treat_col]] %in% c("Control", "Burned")
  ]
  
  ps <- prune_samples(keep_samples, ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
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
  
  sd_genus <- data.frame(sample_data(ps_genus), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(sd_genus) <- sample_names(ps_genus)
  sd_genus$SampleID <- rownames(sd_genus)
  
  get_genera_set <- function(region_value, treatment_value, year_value) {
    
    samples_use <- sd_genus$SampleID[
      sd_genus[[region_col]] == region_value &
        sd_genus[[treat_col]] == treatment_value &
        sd_genus[[year_col]] == year_value
    ]
    
    samples_use <- intersect(samples_use, colnames(otu_mat))
    
    if (length(samples_use) == 0) {
      return(character(0))
    }
    
    taxa_present <- rownames(otu_mat)[
      rowSums(otu_mat[, samples_use, drop = FALSE]) > 0
    ]
    
    genera_present <- tax_genus$Genus_clean[
      match(taxa_present, tax_genus$ASV_ID)
    ]
    
    genera_present <- sort(unique(na.omit(genera_present)))
    
    return(genera_present)
  }
  
  years_amazonas <- c("2", "4", "6", "10")
  years_cusco <- c("1", "5", "10", "15", "20")
  
  pal_amazonas <- c(
    "2"  = "#CBE8C8",
    "4"  = "#A8D9D1",
    "6"  = "#8CC7C2",
    "10" = "#2F9C83"
  )
  
  pal_cusco <- c(
    "1"  = "#CFE4EF",
    "5"  = "#9EC8E0",
    "10" = "#65ACD2",
    "15" = "#3F91C8",
    "20" = "#1F6FC2"
  )
  
  build_sets <- function(region_value, treatment_value, years_vec) {
    sets <- list()
    
    for (yr in years_vec) {
      sets[[yr]] <- get_genera_set(region_value, treatment_value, yr)
    }
    
    sets <- sets[years_vec]
    return(sets)
  }
  
  sets_amz_control <- build_sets("Amazonas", "Control", years_amazonas)
  sets_amz_burned  <- build_sets("Amazonas", "Burned",  years_amazonas)
  
  sets_cus_control <- build_sets("Cusco", "Control", years_cusco)
  sets_cus_burned  <- build_sets("Cusco", "Burned",  years_cusco)
  
  p_amz_control <- make_multiyear_venn_plot(
    sets = sets_amz_control,
    title_text = "Control",
    palette_named = pal_amazonas,
    region_value = "Amazonas"
  )
  
  p_amz_burned <- make_multiyear_venn_plot(
    sets = sets_amz_burned,
    title_text = "Burned",
    palette_named = pal_amazonas,
    region_value = "Amazonas"
  )
  
  p_cus_control <- make_multiyear_venn_plot(
    sets = sets_cus_control,
    title_text = "Control",
    palette_named = pal_cusco,
    region_value = "Cusco"
  )
  
  p_cus_burned <- make_multiyear_venn_plot(
    sets = sets_cus_burned,
    title_text = "Burned",
    palette_named = pal_cusco,
    region_value = "Cusco"
  )
  
  leg_amz <- make_legend_plot(years_amazonas, pal_amazonas)
  leg_cus <- make_legend_plot(years_cusco, pal_cusco)
  
  title_amz <- ggdraw() +
    draw_label("Amazonas", x = 0.02, hjust = 0, fontface = "bold", size = 22)
  
  title_cus <- ggdraw() +
    draw_label("Cusco", x = 0.02, hjust = 0, fontface = "bold", size = 22)
  
  row_amz <- plot_grid(
    p_amz_control,
    p_amz_burned,
    leg_amz,
    nrow = 1,
    rel_widths = c(1, 1, 0.18),
    align = "h"
  )
  
  row_cus <- plot_grid(
    p_cus_control,
    p_cus_burned,
    leg_cus,
    nrow = 1,
    rel_widths = c(1, 1, 0.18),
    align = "h"
  )
  
  fig_amz <- plot_grid(
    title_amz,
    row_amz,
    ncol = 1,
    rel_heights = c(0.08, 1)
  )
  
  fig_cus <- plot_grid(
    title_cus,
    row_cus,
    ncol = 1,
    rel_heights = c(0.08, 1)
  )
  
  fig_body <- plot_grid(
    fig_amz,
    fig_cus,
    ncol = 1,
    rel_heights = c(1, 1)
  )
  
  main_title <- ggdraw() +
    draw_label(
      paste0(marker_label, ": identified genera by treatment and years since fire"),
      fontface = "bold",
      size = 20
    )
  
  fig_final <- plot_grid(
    main_title,
    fig_body,
    ncol = 1,
    rel_heights = c(0.07, 1)
  )
  
  print(fig_final)
  
  make_set_table <- function(sets, region_value, treatment_value) {
    data.frame(
      Marker = marker_name,
      Region = region_value,
      Treatment = treatment_value,
      Year_after_wildfire = names(sets),
      N_identified_genera = vapply(sets, length, integer(1)),
      stringsAsFactors = FALSE
    )
  }
  
  set_sizes_table <- bind_rows(
    make_set_table(sets_amz_control, "Amazonas", "Control"),
    make_set_table(sets_amz_burned,  "Amazonas", "Burned"),
    make_set_table(sets_cus_control, "Cusco",    "Control"),
    make_set_table(sets_cus_burned,  "Cusco",    "Burned")
  )
  
  make_presence_table <- function(sets, region_value, treatment_value) {
    all_genera <- sort(unique(unlist(sets)))
    
    if (length(all_genera) == 0) {
      return(data.frame())
    }
    
    out <- data.frame(
      Marker = marker_name,
      Region = region_value,
      Treatment = treatment_value,
      Genus = all_genera,
      stringsAsFactors = FALSE
    )
    
    for (yr in names(sets)) {
      out[[paste0("Year_", yr)]] <- all_genera %in% sets[[yr]]
    }
    
    return(out)
  }
  
  presence_table <- bind_rows(
    make_presence_table(sets_amz_control, "Amazonas", "Control"),
    make_presence_table(sets_amz_burned,  "Amazonas", "Burned"),
    make_presence_table(sets_cus_control, "Cusco",    "Control"),
    make_presence_table(sets_cus_burned,  "Cusco",    "Burned")
  )
  
  print(set_sizes_table)
  
  outdir <- file.path(base, "results/plots", out_folder_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    set_sizes_table,
    file.path(outdir, paste0(file_prefix, "_multiyear_venn_identified_genera_set_sizes.csv")),
    row.names = FALSE
  )
  
  write.csv(
    presence_table,
    file.path(outdir, paste0(file_prefix, "_multiyear_venn_identified_genera_presence.csv")),
    row.names = FALSE
  )
  
  write.csv(
    tax_genus[!is.na(tax_genus$Genus_clean), ],
    file.path(outdir, paste0(file_prefix, "_ASVs_with_identified_genus_used.csv")),
    row.names = FALSE
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_Multiyear_Venn_identified_genera_Combined.png")),
    fig_final,
    width = 16,
    height = 11.5,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_Multiyear_Venn_identified_genera_Combined.pdf")),
    fig_final,
    width = 16,
    height = 11.5,
    bg = "white"
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_Multiyear_Venn_identified_genera_Amazonas.png")),
    fig_amz,
    width = 16,
    height = 5.8,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    file.path(outdir, paste0(file_prefix, "_Multiyear_Venn_identified_genera_Cusco.png")),
    fig_cus,
    width = 16,
    height = 5.8,
    dpi = 600,
    bg = "white"
  )
  
  cat("\nArchivos guardados en:\n")
  cat(outdir, "\n")
  
  return(fig_final)
}


p_bacteria <- run_multiyear_venn_marker(
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
  out_folder_name = "venn_multiyear_BAC_identified_genera_corrected_percent",
  file_prefix = "BAC"
)

p_archaea <- run_multiyear_venn_marker(
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
  out_folder_name = "venn_multiyear_ARQ_identified_genera_corrected_percent",
  file_prefix = "ARQ"
)

p_its <- run_multiyear_venn_marker(
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
  out_folder_name = "venn_multiyear_ITS_identified_genera_corrected_percent",
  file_prefix = "ITS"
)

cat("\n=====================================================\n")
cat("Venn diagrams by years since fire generated\n")
cat("Criterion: identified genera only\n")
cat("Percentages computed as region / total unique genera in the panel\n")
cat("Hora fin:", as.character(Sys.time()), "\n")
cat("=====================================================\n")


