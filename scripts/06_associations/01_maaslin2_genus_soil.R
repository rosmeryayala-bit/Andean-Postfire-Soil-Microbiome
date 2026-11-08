# Microbiome-soil associations
# Genus-level associations between microbial abundance and soil properties using MaAsLin2.
# Input: project-specific phyloseq objects, site metadata, or soil data tables.
# Output: analysis tables and figures as defined by the project directory settings.
options(stringsAsFactors = FALSE)

pkgs <- c("phyloseq","Maaslin2","ggplot2","ggrepel","patchwork","openxlsx")
miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if(length(miss) > 0){
  stop("Faltan paquetes: ", paste(miss, collapse = ", "))
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(Maaslin2)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(openxlsx)
})

root <- "/home/usuario1/Rosmery/Tesis Rosmery-Metabarcoding"

soil_candidates_files <- c(
  file.path(root, "RESULTADOS_SOIL_UNIQUE_UNITS_FINAL/TABLES/soil_metadata_unique_units.csv"),
  file.path(root, "RESULTADOS_SOIL_ALL_PARAMETERS_FINAL/TABLES/soil_metadata_clean.csv"),
  file.path(root, "RESULTADOS_SOIL_ONLY_FAST_FINAL/TABLES/soil_metadata_clean.csv")
)

soil_file <- soil_candidates_files[file.exists(soil_candidates_files)][1]

if(is.na(soil_file)){
  stop("Unique-unit soil data file not found.")
}

ps_files <- list(
  Bacteria = file.path(root, "RESULT BACTERIA/data/processed_data/phyloseq/ps_sin_prevalencia.rds"),
  Archaea  = file.path(root, "RESULT ARQUEA/data/processed_data/phyloseq/ps_sin_prevalencia_ARQ.rds"),
  Fungi    = file.path(root, "RESULT FUNGAL/data/processed_data/phyloseq/ps_raw_ITS.rds")
)

for(p in ps_files){
  if(!file.exists(p)) stop("No existe: ", p)
}

outdir <- file.path(root, "RESULTADOS_MAASLIN2_SOIL_Q1_FIXED_FINAL")

dirs <- c(
  outdir,
  file.path(outdir, "01_INPUT_TABLES"),
  file.path(outdir, "02_RAW_MAASLIN2"),
  file.path(outdir, "03_TABLES"),
  file.path(outdir, "04_FIGURES_PNG"),
  file.path(outdir, "05_FIGURES_PDF_EDITABLE"),
  file.path(outdir, "06_EXCEL")
)

invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

q_primary <- 0.10
q_exploratory <- 0.25
min_prev <- 0.20
min_mean_ab <- 0.0001

soil_vars_target <- c(
  "pH",
  "Electrical_conductivity",
  "Nitrogen",
  "Phosphorus",
  "Available_K",
  "Organic_carbon",
  "Organic_matter",
  "CEC",
  "Calcium",
  "Magnesium",
  "Exchangeable_K",
  "Sodium",
  "Aluminum",
  "Base_saturation"
)

pretty_names <- c(
  pH = "pH",
  Electrical_conductivity = "Electrical conductivity",
  Nitrogen = "Nitrogen",
  Phosphorus = "Phosphorus",
  Available_K = "Available K",
  Organic_carbon = "Organic carbon",
  Organic_matter = "Organic matter",
  CEC = "CEC",
  Calcium = "Calcium",
  Magnesium = "Magnesium",
  Exchangeable_K = "Exchangeable K",
  Sodium = "Sodium",
  Aluminum = "Aluminum",
  Base_saturation = "Base saturation"
)

pretty_var <- function(x){
  out <- x
  out[x %in% names(pretty_names)] <- pretty_names[x[x %in% names(pretty_names)]]
  out
}

norm_txt <- function(x){
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  gsub("[^a-z0-9]", "", x)
}

find_col <- function(df, candidates){
  nms <- names(df)
  idx <- match(norm_txt(candidates), norm_txt(nms))
  idx <- idx[!is.na(idx)]
  if(length(idx) == 0) return(NA_character_)
  nms[idx[1]]
}

parse_num <- function(x){
  if(is.numeric(x)) return(as.numeric(x))
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(gsub("[^0-9\\.\\-]", "", x)))
}

clean_sample <- function(x){
  x <- as.character(x)
  x <- gsub("\\.raw$", "", x)
  x <- gsub("_(BAC|ARQ|FUN|ITS)$", "", x, ignore.case = TRUE)
  x <- gsub("^A_", "A", x)
  x <- gsub("^C_", "C", x)
  x
}

std_region <- function(x){
  x0 <- tolower(as.character(x))
  ifelse(grepl("amazon", x0), "Amazonas",
         ifelse(grepl("cusco", x0), "Cusco", as.character(x)))
}

std_treatment <- function(x){
  x0 <- tolower(as.character(x))
  ifelse(grepl("wild|burn|incend|fire|quem", x0), "Burned",
         ifelse(grepl("control|ctrl", x0), "Control", as.character(x)))
}

make_key <- function(region, treatment, year){
  paste(std_region(region), std_treatment(treatment), parse_num(year), sep = "_")
}

valid_taxon <- function(x){
  x <- as.character(x)
  bad <- is.na(x) |
    x == "" |
    grepl("uncultured|unclassified|unknown|unidentified|metagenome|norank|ambiguous|incertae|chloroplast|mitochondria", x, ignore.case = TRUE)
  !bad
}

feature_clean <- function(x){
  x <- gsub("^Bacteria[\\._]+", "", x)
  x <- gsub("^Archaea[\\._]+", "", x)
  x <- gsub("^Fungi[\\._]+", "", x)
  x <- gsub("\\.", " ", x)
  x
}

theme_q1 <- function(){
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.30),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.70),
      strip.background = element_rect(fill = "grey93", color = "grey55"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
}


soil_u <- read.csv(soil_file, check.names = FALSE)
names(soil_u) <- make.unique(names(soil_u), sep = "_dup")

region_col_s <- find_col(soil_u, c("Region", "Región"))
treat_col_s  <- find_col(soil_u, c("Treatment", "Tratamiento", "Tratament", "Condition"))
year_col_s   <- find_col(soil_u, c("Year_num", "Year", "Years_since_fire", "Years since fire", "Año", "Anio"))

if(any(is.na(c(region_col_s, treat_col_s, year_col_s)))){
  stop("Required soil columns not found: Region/Treatment/Year.")
}

for(v in soil_vars_target){
  if(v %in% names(soil_u)){
    soil_u[[v]] <- parse_num(soil_u[[v]])
  }
}

soil_vars <- soil_vars_target[soil_vars_target %in% names(soil_u)]
soil_vars <- soil_vars[sapply(soil_u[soil_vars], function(x){
  sum(!is.na(x)) >= 6 && sd(x, na.rm = TRUE) > 0
})]

if(length(soil_vars) < 2){
  stop("Insufficient valid soil properties.")
}

soil_u$Region_std <- std_region(soil_u[[region_col_s]])
soil_u$Treatment_std <- std_treatment(soil_u[[treat_col_s]])
soil_u$Year_num <- parse_num(soil_u[[year_col_s]])
soil_u$SoilUnitKey <- make_key(soil_u$Region_std, soil_u$Treatment_std, soil_u$Year_num)

soil_u <- soil_u[!duplicated(soil_u$SoilUnitKey), ]

write.csv(
  soil_u,
  file.path(outdir, "01_INPUT_TABLES", "soil_unique_units_used_for_MaAsLin2.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(Soil_variable = soil_vars, Label = pretty_var(soil_vars)),
  file.path(outdir, "01_INPUT_TABLES", "soil_variables_used_for_MaAsLin2.csv"),
  row.names = FALSE
)

message("Soil file: ", soil_file)
message("Unique soil units: ", nrow(soil_u))
message("Soil variables:")
print(soil_vars)


ps_to_soilunit_genus <- function(ps, marker){

  sample_names(ps) <- clean_sample(sample_names(ps))

  meta <- as(sample_data(ps), "data.frame")
  meta <- data.frame(meta, check.names = FALSE, stringsAsFactors = FALSE)
  names(meta) <- make.unique(names(meta), sep = "_dup")
  meta$SampleID <- clean_sample(rownames(meta))

  region_col <- find_col(meta, c("Region", "Región"))
  treat_col  <- find_col(meta, c("Treatment", "Tratamiento", "Tratament", "Condition"))
  year_col   <- find_col(meta, c("Years_since_fire", "Years since fire", "Year", "Year_num", "Año", "Anio"))

  if(any(is.na(c(region_col, treat_col, year_col)))){
    stop(marker, ": Region/Treatment/Year not found in sample metadata.")
  }

  meta$Region_std <- std_region(meta[[region_col]])
  meta$Treatment_std <- std_treatment(meta[[treat_col]])
  meta$Year_num <- parse_num(meta[[year_col]])
  meta$SoilUnitKey <- make_key(meta$Region_std, meta$Treatment_std, meta$Year_num)

  common_samples <- intersect(sample_names(ps), meta$SampleID)
  ps <- prune_samples(common_samples, ps)
  meta <- meta[match(common_samples, meta$SampleID), , drop = FALSE]

  tax <- as(tax_table(ps), "matrix")
  tax_cols <- colnames(tax)
  genus_col <- tax_cols[match(norm_txt("Genus"), norm_txt(tax_cols))]

  if(is.na(genus_col)){
    stop(marker, ": Genus column not found.")
  }

  keep_tax <- valid_taxon(tax[, genus_col])
  ps <- prune_taxa(keep_tax, ps)

  ps_g <- tax_glom(ps, taxrank = genus_col, NArm = FALSE)

  taxg <- as(tax_table(ps_g), "matrix")
  genera <- as.character(taxg[, genus_col])
  genera <- make.unique(genera)

  otu <- as(otu_table(ps_g), "matrix")
  if(taxa_are_rows(ps_g)) otu <- t(otu)

  rownames(otu) <- clean_sample(rownames(otu))
  colnames(otu) <- make.names(paste(marker, genera, sep = "__"), unique = TRUE)

  otu_rel <- sweep(otu, 1, rowSums(otu, na.rm = TRUE), "/")
  otu_rel[!is.finite(otu_rel)] <- 0

  map <- meta[, c("SampleID", "SoilUnitKey"), drop = FALSE]
  map <- map[match(rownames(otu_rel), map$SampleID), , drop = FALSE]

  otu_df <- data.frame(SoilUnitKey = map$SoilUnitKey, otu_rel, check.names = FALSE)

  otu_unit <- aggregate(
    otu_df[, colnames(otu_rel), drop = FALSE],
    by = list(SoilUnitKey = otu_df$SoilUnitKey),
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  rownames(otu_unit) <- otu_unit$SoilUnitKey
  otu_unit$SoilUnitKey <- NULL

  otu_unit <- otu_unit[rownames(otu_unit) %in% soil_u$SoilUnitKey, , drop = FALSE]

  prev <- colMeans(otu_unit > 0, na.rm = TRUE)
  mean_ab <- colMeans(otu_unit, na.rm = TRUE)
  keep <- prev >= min_prev & mean_ab >= min_mean_ab

  otu_unit <- otu_unit[, keep, drop = FALSE]

  write.csv(
    data.frame(SoilUnitKey = rownames(otu_unit), otu_unit, check.names = FALSE),
    file.path(outdir, "01_INPUT_TABLES", paste0(marker, "_genus_relative_abundance_by_soil_unit.csv")),
    row.names = FALSE
  )

  message(marker, ": ", nrow(otu_unit), " soil units × ", ncol(otu_unit), " genera retained")

  otu_unit
}

features <- list()

for(marker in names(ps_files)){
  message("\nPreparing marker: ", marker)
  ps <- readRDS(ps_files[[marker]])
  features[[marker]] <- ps_to_soilunit_genus(ps, marker)
}


run_maaslin_one <- function(marker, ft, reg, soil_var){

  md <- soil_u[soil_u$Region_std == reg, , drop = FALSE]
  common <- intersect(rownames(ft), md$SoilUnitKey)

  md <- md[match(common, md$SoilUnitKey), , drop = FALSE]
  ft2 <- ft[common, , drop = FALSE]

  if(nrow(ft2) < 6 || ncol(ft2) < 2){
    return(NULL)
  }

  x <- md[[soil_var]]
  if(sd(x, na.rm = TRUE) == 0 || sum(!is.na(x)) < 6){
    return(NULL)
  }

  md_run <- data.frame(
    SoilPredictor = as.numeric(scale(x)),
    Treatment = factor(md$Treatment_std, levels = c("Control", "Burned")),
    Year = as.numeric(md$Year_num),
    row.names = md$SoilUnitKey
  )

  ft2 <- ft2[rownames(md_run), , drop = FALSE]

  keep_features <- apply(ft2, 2, function(z){
    sd(z, na.rm = TRUE) > 0 && mean(z > 0, na.rm = TRUE) >= min_prev
  })

  ft2 <- ft2[, keep_features, drop = FALSE]

  if(nrow(ft2) < 6 || ncol(ft2) < 2){
    return(NULL)
  }

  fixed_effects <- c("SoilPredictor")

  if(length(unique(md_run$Treatment)) > 1){
    fixed_effects <- c(fixed_effects, "Treatment")
  }

  if(sd(md_run$Year, na.rm = TRUE) > 0){
    fixed_effects <- c(fixed_effects, "Year")
  }

  refs <- NULL
  if("Treatment" %in% fixed_effects){
    refs <- c("Treatment,Control")
  }

  out_sub <- file.path(outdir, "02_RAW_MAASLIN2", marker, reg, soil_var)
  dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)

  Maaslin2(
    input_data = ft2,
    input_metadata = md_run,
    output = out_sub,
    fixed_effects = fixed_effects,
    reference = refs,
    normalization = "NONE",
    transform = "LOG",
    analysis_method = "LM",
    correction = "BH",
    standardize = FALSE,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    cores = 1
  )

  res_file <- file.path(out_sub, "all_results.tsv")

  if(!file.exists(res_file)){
    return(NULL)
  }

  res <- read.delim(res_file, check.names = FALSE)
  res$Marker <- marker
  res$Region <- reg
  res$Soil_variable <- soil_var
  res$Soil_variable_label <- pretty_var(soil_var)
  res$Model <- paste("Genus ~", soil_var, "+ Treatment + Year")
  res
}

all_results <- list()
settings_list <- list()

for(marker in names(features)){
  for(reg in c("Amazonas","Cusco")){
    for(sv in soil_vars){

      message("Running MaAsLin2: ", marker, " | ", reg, " | ", sv)

      tmp <- tryCatch(
        run_maaslin_one(marker, features[[marker]], reg, sv),
        error = function(e){
          message("ERROR: ", marker, " | ", reg, " | ", sv, " -> ", e$message)
          NULL
        }
      )

      if(!is.null(tmp)){
        all_results[[length(all_results) + 1]] <- tmp
      }

      settings_list[[length(settings_list) + 1]] <- data.frame(
        Marker = marker,
        Region = reg,
        Taxonomic_rank = "Genus",
        Soil_variable = sv,
        Model = paste("Genus ~", sv, "+ Treatment + Year"),
        Normalization = "NONE",
        Transformation = "LOG",
        Method = "LM",
        Prevalence_threshold = min_prev,
        Mean_abundance_threshold = min_mean_ab,
        Primary_q_threshold = q_primary,
        Exploratory_q_threshold = q_exploratory,
        stringsAsFactors = FALSE
      )
    }
  }
}

if(length(all_results) == 0){
  stop("MaAsLin2 no produjo resultados.")
}

res_all <- do.call(rbind, all_results)
settings_df <- do.call(rbind, settings_list)

soil_res <- res_all[res_all$metadata == "SoilPredictor", , drop = FALSE]

soil_res$coef <- as.numeric(soil_res$coef)
soil_res$pval <- as.numeric(soil_res$pval)
soil_res$qval <- as.numeric(soil_res$qval)

soil_res$Feature_clean <- feature_clean(soil_res$feature)
soil_res$Direction <- ifelse(soil_res$coef > 0, "Positive", "Negative")
soil_res$Significance <- ifelse(
  soil_res$qval <= q_primary,
  "Primary FDR q <= 0.10",
  ifelse(soil_res$qval <= q_exploratory, "Exploratory FDR q <= 0.25", "Not significant")
)

soil_res <- soil_res[order(soil_res$qval, soil_res$pval), ]

sig_primary <- soil_res[soil_res$qval <= q_primary, , drop = FALSE]
sig_exploratory <- soil_res[soil_res$qval <= q_exploratory, , drop = FALSE]

summary_list <- list()
idx <- 1

for(marker in unique(soil_res$Marker)){
  for(reg in unique(soil_res$Region)){
    for(sv in unique(soil_res$Soil_variable)){

      d <- soil_res[soil_res$Marker == marker & soil_res$Region == reg & soil_res$Soil_variable == sv, , drop = FALSE]

      if(nrow(d) == 0) next

      d_sig <- d[d$qval <= q_exploratory, , drop = FALSE]

      top_pos <- NA
      top_neg <- NA

      if(any(d_sig$coef > 0)){
        top_pos <- d_sig$Feature_clean[which.max(ifelse(d_sig$coef > 0, d_sig$coef, -Inf))]
      }

      if(any(d_sig$coef < 0)){
        top_neg <- d_sig$Feature_clean[which.min(ifelse(d_sig$coef < 0, d_sig$coef, Inf))]
      }

      summary_list[[idx]] <- data.frame(
        Marker = marker,
        Region = reg,
        Soil_variable = sv,
        Soil_variable_label = pretty_var(sv),
        N_tested = nrow(d),
        N_primary_q0.10 = sum(d$qval <= q_primary, na.rm = TRUE),
        N_exploratory_q0.25 = sum(d$qval <= q_exploratory, na.rm = TRUE),
        Top_positive_feature = top_pos,
        Top_negative_feature = top_neg,
        stringsAsFactors = FALSE
      )

      idx <- idx + 1
    }
  }
}

summary_assoc <- do.call(rbind, summary_list)


write.csv(settings_df, file.path(outdir, "03_TABLES", "Table_S16_MaAsLin2_model_settings.csv"), row.names = FALSE)
write.csv(sig_primary, file.path(outdir, "03_TABLES", "Table_S17_primary_significant_associations_q0.10.csv"), row.names = FALSE)
write.csv(sig_exploratory, file.path(outdir, "03_TABLES", "Table_S18_exploratory_associations_q0.25.csv"), row.names = FALSE)
write.csv(summary_assoc, file.path(outdir, "03_TABLES", "Table_S19_summary_soil_microbiome_associations.csv"), row.names = FALSE)
write.csv(soil_res, file.path(outdir, "03_TABLES", "Dataset_S4_full_MaAsLin2_soil_associations.csv"), row.names = FALSE)


dir.create(file.path(outdir, "04_FIGURES_PNG"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "05_FIGURES_PDF_EDITABLE"), recursive = TRUE, showWarnings = FALSE)

plot_base <- sig_primary
plot_level <- "Primary FDR q <= 0.10"

if(nrow(plot_base) == 0){
  plot_base <- sig_exploratory
  plot_level <- "Exploratory FDR q <= 0.25"
}

if(nrow(plot_base) == 0){
  plot_base <- soil_res[order(soil_res$qval, soil_res$pval), , drop = FALSE]
  plot_base <- plot_base[1:min(80, nrow(plot_base)), , drop = FALSE]
  plot_level <- "Top-ranked associations; no FDR-significant associations"
}

plot_base <- plot_base[order(plot_base$qval, -abs(plot_base$coef)), , drop = FALSE]
plot_base <- plot_base[1:min(90, nrow(plot_base)), , drop = FALSE]

plot_base$Feature_label <- paste(plot_base$Marker, plot_base$Feature_clean, sep = " | ")

p_heat <- ggplot(
  plot_base,
  aes(x = Soil_variable_label, y = reorder(Feature_label, coef), fill = coef)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(Marker ~ Region, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#7DA0CA",
    mid = "#F7F4ED",
    high = "#D97B66",
    midpoint = 0,
    name = "Coefficient"
  ) +
  labs(
    title = "Soil-microbiome associations identified by MaAsLin2",
    subtitle = plot_level,
    x = "Soil physicochemical variable",
    y = "Microbial genus"
  ) +
  theme_q1() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(outdir, "04_FIGURES_PNG", "Figure_S16_MaAsLin2_soil_microbiome_heatmap.png"),
  p_heat,
  width = 14,
  height = 10,
  dpi = 600
)

ggsave(
  file.path(outdir, "05_FIGURES_PDF_EDITABLE", "Figure_S16_MaAsLin2_soil_microbiome_heatmap.pdf"),
  p_heat,
  width = 14,
  height = 10
)

top_bar <- plot_base[order(-abs(plot_base$coef), plot_base$qval), , drop = FALSE]
top_bar <- top_bar[1:min(50, nrow(top_bar)), , drop = FALSE]
top_bar$Label <- paste(
  top_bar$Marker,
  top_bar$Region,
  top_bar$Soil_variable_label,
  top_bar$Feature_clean,
  sep = " | "
)

p_bar <- ggplot(
  top_bar,
  aes(x = coef, y = reorder(Label, coef), fill = Direction)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c("Positive" = "#D97B66", "Negative" = "#7DA0CA")) +
  labs(
    title = "Top soil-associated microbial genera",
    subtitle = plot_level,
    x = "MaAsLin2 coefficient",
    y = NULL,
    fill = "Direction"
  ) +
  theme_q1()

ggsave(
  file.path(outdir, "04_FIGURES_PNG", "Figure_S17_top_MaAsLin2_coefficients.png"),
  p_bar,
  width = 14,
  height = max(7, 0.24 * nrow(top_bar)),
  dpi = 600
)

ggsave(
  file.path(outdir, "05_FIGURES_PDF_EDITABLE", "Figure_S17_top_MaAsLin2_coefficients.pdf"),
  p_bar,
  width = 14,
  height = max(7, 0.24 * nrow(top_bar))
)

make_detail_panel <- function(marker_sel, reg_sel){

  d_all <- soil_res[
    soil_res$Marker == marker_sel & soil_res$Region == reg_sel,
    ,
    drop = FALSE
  ]

  if(nrow(d_all) == 0){
    return(NULL)
  }

  top_vars <- summary_assoc[
    summary_assoc$Marker == marker_sel & summary_assoc$Region == reg_sel,
    ,
    drop = FALSE
  ]

  if(nrow(top_vars) == 0){
    return(NULL)
  }

  top_vars <- top_vars[order(-top_vars$N_exploratory_q0.25, top_vars$Soil_variable_label), , drop = FALSE]
  top_vars <- top_vars$Soil_variable[1:min(4, nrow(top_vars))]

  plist <- list()

  for(sv in top_vars){

    dv <- d_all[d_all$Soil_variable == sv, , drop = FALSE]

    if(nrow(dv) == 0){
      next
    }

    dv$Status <- ifelse(
      dv$qval <= q_primary,
      "q <= 0.10",
      ifelse(dv$qval <= q_exploratory, "q <= 0.25", "ns")
    )

    lab_df <- dv[dv$qval <= q_exploratory, , drop = FALSE]

    if(nrow(lab_df) > 0){
      lab_df <- lab_df[order(lab_df$qval, -abs(lab_df$coef)), , drop = FALSE]
      lab_df <- lab_df[1:min(6, nrow(lab_df)), , drop = FALSE]
      dv$Label <- ifelse(dv$feature %in% lab_df$feature, dv$Feature_clean, NA)
    } else {
      dv$Label <- NA
    }

    p1 <- ggplot(dv, aes(x = coef, y = -log10(pmax(qval, 1e-300)))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
      geom_hline(yintercept = -log10(q_primary), linetype = "dotted", color = "grey55") +
      geom_point(aes(color = Status), alpha = 0.85, size = 1.9) +
      ggrepel::geom_text_repel(aes(label = Label), size = 2.6, max.overlaps = Inf) +
      scale_color_manual(
        values = c(
          "q <= 0.10" = "#08306B",
          "q <= 0.25" = "#6BAED6",
          "ns" = "grey82"
        )
      ) +
      labs(
        title = pretty_var(sv),
        x = "MaAsLin2 coefficient",
        y = expression(-log[10](q~value)),
        color = "FDR"
      ) +
      theme_q1() +
      theme(legend.position = "bottom")

    bar_df <- dv[dv$qval <= q_exploratory, , drop = FALSE]

    if(nrow(bar_df) > 0){

      bar_df <- bar_df[order(-abs(bar_df$coef), bar_df$qval), , drop = FALSE]
      bar_df <- bar_df[1:min(8, nrow(bar_df)), , drop = FALSE]

      p2 <- ggplot(bar_df, aes(x = coef, y = reorder(Feature_clean, coef), fill = Direction)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
        geom_col(width = 0.75) +
        scale_fill_manual(values = c("Positive" = "#D97B66", "Negative" = "#7DA0CA")) +
        labs(
          x = "MaAsLin2 coefficient",
          y = NULL,
          fill = "Direction"
        ) +
        theme_q1() +
        theme(legend.position = "bottom")

    } else {

      p2 <- ggplot() +
        annotate("text", x = 0, y = 0, label = "No q <= 0.25 associations", size = 4) +
        theme_void()
    }

    plist[[sv]] <- p1 / p2
  }

  if(length(plist) == 0){
    return(NULL)
  }

  panel <- patchwork::wrap_plots(plist, ncol = 2) +
    patchwork::plot_annotation(
      title = paste0(marker_sel, " | ", reg_sel, " | MaAsLin2 soil associations")
    )

  ggsave(
    file.path(outdir, "04_FIGURES_PNG", paste0("Figure_S18_", marker_sel, "_", reg_sel, "_detailed_panels.png")),
    panel,
    width = 16,
    height = 12,
    dpi = 600
  )

  ggsave(
    file.path(outdir, "05_FIGURES_PDF_EDITABLE", paste0("Figure_S18_", marker_sel, "_", reg_sel, "_detailed_panels.pdf")),
    panel,
    width = 16,
    height = 12
  )
}

for(marker in unique(soil_res$Marker)){
  for(reg in unique(soil_res$Region)){
    make_detail_panel(marker, reg)
  }
}


dir.create(file.path(outdir, "06_EXCEL"), recursive = TRUE, showWarnings = FALSE)

wb <- openxlsx::createWorkbook()

add_sheet <- function(name, df){
  name <- substr(name, 1, 31)
  openxlsx::addWorksheet(wb, name)

  if(is.null(df) || nrow(df) == 0){
    openxlsx::writeData(wb, name, data.frame(Message = "No results"))
  } else {
    openxlsx::writeData(wb, name, df)
    openxlsx::freezePane(wb, name, firstRow = TRUE)
    openxlsx::setColWidths(wb, name, cols = 1:ncol(df), widths = "auto")
  }
}

add_sheet("S16_model_settings", settings_df)
add_sheet("S17_primary_q0.10", sig_primary)
add_sheet("S18_exploratory_q0.25", sig_exploratory)
add_sheet("S19_summary", summary_assoc)
add_sheet("Dataset_S4_full", soil_res)

openxlsx::saveWorkbook(
  wb,
  file.path(outdir, "06_EXCEL", "MAASLIN2_SOIL_Q1_FIXED_FINAL_RESULTS.xlsx"),
  overwrite = TRUE
)

message("\nDONE.")
message("Output folder: ", outdir)
message("Excel: ", file.path(outdir, "06_EXCEL", "MAASLIN2_SOIL_Q1_FIXED_FINAL_RESULTS.xlsx"))


base_dir <- "/home/usuario1/Rosmery/Tesis Rosmery-Metabarcoding/RESULTADOS_MAASLIN2_SOIL_Q1_FIXED_FINAL"
dir.exists(base_dir)
