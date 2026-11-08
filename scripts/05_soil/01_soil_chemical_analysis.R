# Soil chemical analyses
# Region-specific edaphic characterisation and statistical models.
# Input: project-specific phyloseq objects, site metadata, or soil data tables.
# Output: analysis tables and figures as defined by the project directory settings.


suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(broom)
})


base <- "~/Tesis Rosmery-Metabarcoding/RESULT BACTERIA"
outdir <- file.path(base, "results/soil_univariate_by_region")
figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)


if (!exists("ps")) stop("The phyloseq object 'ps' has not been loaded.")

soil_df <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
soil_df$SampleID <- rownames(soil_df)

soil_df[] <- lapply(soil_df, function(x) {
  if (is.factor(x)) as.character(x) else x
})

soil_df <- as.data.frame(soil_df, check.names = FALSE, stringsAsFactors = FALSE)

cat("\n=====================================================\n")
cat("COLUMNAS DETECTADAS EN sample_data(ps)\n")
cat("=====================================================\n")
print(colnames(soil_df))
cat("\n")


find_col <- function(df_names, candidates) {
  hit <- candidates[candidates %in% df_names][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}


region_col <- find_col(colnames(soil_df), c("Region", "region"))
treat_col  <- find_col(colnames(soil_df), c("Treatment", "treatment", "Tratamiento", "tratamiento"))
year_col   <- find_col(colnames(soil_df), c("Year", "year", "Years_since_fire", "Years_since_fire.1"))

if (is.na(region_col)) stop("Region column not found.")
if (is.na(treat_col))  stop("Treatment column not found.")
if (is.na(year_col))   stop("Years-since-fire column not found.")


soil_aliases <- list(
  pH   = c("pH"),
  EC   = c("Electrical conductivity (mS/m)", "Electrical.conductivity..mS.m."),
  N    = c("N (mg/g)", "N..mg.g."),
  P    = c("P (ppm)", "P..ppm."),
  K    = c("K (ppm)", "K..ppm."),
  C    = c("C (%)", "C...."),
  OM   = c("Organic matter (%)", "Organic.matter...."),
  CEC  = c("CEC (meq/100g)", "CEC..meq.100g."),
  Ca   = c("Ca (meq/100g)", "Ca..meq.100g."),
  Mg   = c("Mg (meq/100g)", "Mg..meq.100g."),
  K_ex = c("K+  (meq/100g)", "K....meq.100g."),
  Na   = c("Na  (meq/100g)", "Na...meq.100g."),
  Al   = c("Al+3  (meq/100g)", "Al.3...meq.100g."),
  BS   = c("%  Base saturation", "X...Base.saturation")
)

soil_cols <- sapply(soil_aliases, function(x) find_col(colnames(soil_df), x))

missing_soil <- names(soil_cols)[is.na(soil_cols)]
if (length(missing_soil) > 0) {
  cat("ATENCION: No se encontraron estas variables de suelo:\n")
  print(missing_soil)
  cat("\n")
}

soil_cols <- soil_cols[!is.na(soil_cols)]

if (length(soil_cols) == 0) {
  stop("No soil property variables found.")
}

cat("Soil property columns selected for analysis:\n")
print(soil_cols)
cat("\n")


cols_to_take <- c("SampleID", region_col, treat_col, year_col, unname(soil_cols))
cols_to_take <- unique(cols_to_take)

missing_take <- setdiff(cols_to_take, colnames(soil_df))
if (length(missing_take) > 0) {
  stop("Estas columnas no existen en soil_df: ", paste(missing_take, collapse = ", "))
}

soil_use <- soil_df[, cols_to_take, drop = FALSE]
soil_use <- as.data.frame(soil_use, check.names = FALSE, stringsAsFactors = FALSE)


colnames(soil_use)[colnames(soil_use) == region_col] <- "Region"
colnames(soil_use)[colnames(soil_use) == treat_col]  <- "Treatment"
colnames(soil_use)[colnames(soil_use) == year_col]   <- "Years_since_fire"


for (nm in names(soil_cols)) {
  colnames(soil_use)[colnames(soil_use) == soil_cols[[nm]]] <- nm
}


soil_use$Region <- trimws(as.character(soil_use$Region))
soil_use$Treatment <- trimws(as.character(soil_use$Treatment))
soil_use$Years_since_fire <- trimws(as.character(soil_use$Years_since_fire))

soil_use$Treatment <- dplyr::case_when(
  tolower(soil_use$Treatment) %in% c("control") ~ "Control",
  tolower(soil_use$Treatment) %in% c("burned", "wildfire", "fire", "incendio", "quemado") ~ "Burned",
  TRUE ~ soil_use$Treatment
)

soil_use$Region <- factor(soil_use$Region, levels = c("Amazonas", "Cusco"))
soil_use$Treatment <- factor(soil_use$Treatment, levels = c("Control", "Burned"))
soil_use$RegionYear <- paste(soil_use$Region, soil_use$Years_since_fire, sep = "_")

soil_var_names <- names(soil_cols)


for (v in soil_var_names) {
  soil_use[[v]] <- as.numeric(as.character(soil_use[[v]]))
}

pretty_names <- c(
  pH   = "pH",
  EC   = "Electrical conductivity",
  N    = "Nitrogen",
  P    = "Phosphorus",
  K    = "Potassium",
  C    = "Carbon",
  OM   = "Organic matter",
  CEC  = "CEC",
  Ca   = "Calcium",
  Mg   = "Magnesium",
  K_ex = "Exchangeable K",
  Na   = "Sodium",
  Al   = "Aluminum",
  BS   = "Base saturation"
)

cat("=====================================================\n")
cat("REVISION DE DATOS LIMPIOS\n")
cat("=====================================================\n")
cat("Numero de muestras:", nrow(soil_use), "\n")
cat("Variables analizadas:\n")
print(soil_var_names)
cat("\n")

cat("Vista previa de soil_use:\n")
print(head(soil_use))
cat("\n")


soil_long <- tidyr::pivot_longer(
  data = soil_use,
  cols = all_of(soil_var_names),
  names_to = "Variable",
  values_to = "Value"
)

soil_summary <- soil_long %>%
  filter(!is.na(Value)) %>%
  group_by(Region, Years_since_fire, Treatment, Variable) %>%
  summarise(
    n = n(),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Variable_pretty = recode(Variable, !!!pretty_names)) %>%
  select(Region, Years_since_fire, Treatment, Variable, Variable_pretty, everything())

write_csv(
  soil_summary,
  file.path(tabdir, "soil_descriptive_summary_by_region_year_treatment.csv")
)

cat("Tabla descriptiva generada correctamente.\n\n")


run_soil_model <- function(df_region, variable_name) {
  
  dd <- df_region %>%
    dplyr::select(SampleID, Region, Treatment, Years_since_fire, all_of(variable_name)) %>%
    dplyr::rename(Value = all_of(variable_name)) %>%
    dplyr::filter(!is.na(Value))
  
  if (nrow(dd) < 6) {
    return(tibble(
      Variable = variable_name,
      term = NA_character_,
      df = NA_real_,
      sumsq = NA_real_,
      meansq = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_
    ))
  }
  
  dd$Years_since_fire <- factor(
    dd$Years_since_fire,
    levels = sort(unique(dd$Years_since_fire))
  )
  
  fit <- aov(Value ~ Years_since_fire * Treatment, data = dd)
  out <- broom::tidy(fit)
  
  out %>%
    mutate(Variable = variable_name) %>%
    select(Variable, everything())
}


model_results <- bind_rows(
  lapply(levels(soil_use$Region), function(reg) {
    df_reg <- soil_use %>% filter(Region == reg)
    
    bind_rows(
      lapply(soil_var_names, function(v) {
        run_soil_model(df_reg, v) %>%
          mutate(Region = reg)
      })
    )
  })
) %>%
  mutate(
    Variable_pretty = recode(Variable, !!!pretty_names),
    Significance = case_when(
      is.na(p.value) ~ "",
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ "ns"
    )
  )

write_csv(
  model_results,
  file.path(tabdir, "soil_anova_by_region_variable.csv")
)

cat("Resultados de modelos generados correctamente.\n\n")


soil_results_clean <- model_results %>%
  filter(!is.na(term)) %>%
  mutate(
    Term = case_when(
      term == "Years_since_fire" ~ "Year",
      term == "Treatment" ~ "Treatment",
      term == "Years_since_fire:Treatment" ~ "Year × Treatment",
      term == "Residuals" ~ "Residuals",
      TRUE ~ term
    )
  ) %>%
  select(
    Region,
    Variable,
    Variable_pretty,
    Term,
    df,
    statistic,
    p.value,
    Significance
  )

write_csv(
  soil_results_clean,
  file.path(tabdir, "soil_anova_clean_for_thesis.csv")
)

cat("Tabla limpia para tesis generada correctamente.\n\n")


plot_soil_variable <- function(df, variable_name) {
  
  dd <- df %>%
    dplyr::select(Region, Treatment, Years_since_fire, all_of(variable_name)) %>%
    dplyr::rename(Value = all_of(variable_name)) %>%
    dplyr::filter(!is.na(Value))
  
  if (nrow(dd) == 0) return(NULL)
  
  dd$Year <- as.character(dd$Years_since_fire)
  
  p <- ggplot(dd, aes(x = Year, y = Value, fill = Treatment, color = Treatment)) +
    geom_boxplot(
      position = position_dodge(width = 0.72),
      width = 0.58,
      outlier.shape = NA,
      alpha = 0.82
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.72),
      size = 2,
      alpha = 0.85
    ) +
    facet_wrap(~Region, scales = "free_x", nrow = 1) +
    labs(
      title = paste("Soil variable -", pretty_names[[variable_name]]),
      x = "Years since fire",
      y = pretty_names[[variable_name]]
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  safe_name <- gsub("[^A-Za-z0-9]+", "_", variable_name)
  safe_name <- gsub("_+", "_", safe_name)
  safe_name <- gsub("^_|_$", "", safe_name)
  
  ggsave(
    filename = file.path(figdir, paste0("soil_", safe_name, ".png")),
    plot = p,
    width = 11,
    height = 5.5,
    dpi = 500
  )
  
  ggsave(
    filename = file.path(figdir, paste0("soil_", safe_name, ".pdf")),
    plot = p,
    width = 11,
    height = 5.5
  )
  
  return(p)
}

plots_soil <- lapply(soil_var_names, function(v) plot_soil_variable(soil_use, v))


distinct_check <- soil_use %>%
  group_by(Region, Years_since_fire, Treatment) %>%
  summarise(across(all_of(soil_var_names), ~n_distinct(.x[!is.na(.x)])), .groups = "drop")

write_csv(
  distinct_check,
  file.path(tabdir, "soil_distinct_values_check.csv")
)

cat("Tabla de verificacion de valores distintos por grupo guardada.\n\n")


cat("\n=====================================================\n")
cat("ANALISIS 1 COMPLETADO CORRECTAMENTE\n")
cat("Archivos guardados en:\n")
cat(" - ", tabdir, "\n")
cat(" - ", figdir, "\n")
cat("=====================================================\n")


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})


if (!exists("soil_use")) stop("No existe 'soil_use' en el entorno.")

soil_var_names <- c("pH", "EC", "N", "P", "K", "C", "OM", "CEC", "Ca", "Mg", "K_ex", "Na", "Al", "BS")


missing_vars <- setdiff(soil_var_names, colnames(soil_use))
if (length(missing_vars) > 0) {
  stop("Faltan estas variables en soil_use: ", paste(missing_vars, collapse = ", "))
}


soil_collapsed <- soil_use %>%
  group_by(Region, Years_since_fire, Treatment, RegionYear) %>%
  summarise(
    across(all_of(soil_var_names), ~mean(.x, na.rm = TRUE)),
    n_rows_collapsed = n(),
    .groups = "drop"
  )


base <- "~/Tesis Rosmery-Metabarcoding/RESULT BACTERIA"
outdir <- file.path(base, "results/soil_collapsed")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write_csv(soil_collapsed, file.path(outdir, "soil_collapsed_by_region_year_treatment.csv"))

cat("\n=====================================================\n")
cat("TABLA DE SUELO COLAPSADA\n")
cat("=====================================================\n")
print(soil_collapsed)
cat("\n")
cat("Filas finales:", nrow(soil_collapsed), "\n")
cat("Archivo guardado en:\n")
cat(file.path(outdir, "soil_collapsed_by_region_year_treatment.csv"), "\n")
cat("=====================================================\n")


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
})


base <- "~/Tesis Rosmery-Metabarcoding/RESULT BACTERIA"
indir <- file.path(base, "results/soil_collapsed")
outdir <- file.path(base, "results/soil_pca")
figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

soil_file <- file.path(indir, "soil_collapsed_by_region_year_treatment.csv")
if (!file.exists(soil_file)) {
  stop("No existe: ", soil_file)
}

soil_collapsed <- read_csv(soil_file, show_col_types = FALSE)


soil_collapsed$Region <- factor(soil_collapsed$Region, levels = c("Amazonas", "Cusco"))
soil_collapsed$Treatment <- factor(soil_collapsed$Treatment, levels = c("Control", "Burned"))
soil_collapsed$Years_since_fire <- as.character(soil_collapsed$Years_since_fire)

soil_vars <- c("pH", "EC", "N", "P", "K", "C", "OM", "CEC", "Ca", "Mg", "K_ex", "Na", "Al", "BS")

missing_vars <- setdiff(soil_vars, colnames(soil_collapsed))
if (length(missing_vars) > 0) {
  stop("Faltan variables en soil_collapsed: ", paste(missing_vars, collapse = ", "))
}


na_count <- sapply(soil_collapsed[soil_vars], function(x) sum(is.na(x)))
cat("NA por variable:\n")
print(na_count)
cat("\n")


year_colors <- c(
  "1"  = "#B3E5FC",
  "2"  = "#A7E9AF",
  "4"  = "#B45A73",
  "5"  = "#FFC1D5",
  "6"  = "#FFD6A5",
  "10" = "#4CB5AB",
  "15" = "#C7D8A5",
  "20" = "#01579B"
)

pretty_names <- c(
  pH   = "pH",
  EC   = "Electrical conductivity",
  N    = "Nitrogen",
  P    = "Phosphorus",
  K    = "Potassium",
  C    = "Carbon",
  OM   = "Organic matter",
  CEC  = "CEC",
  Ca   = "Calcium",
  Mg   = "Magnesium",
  K_ex = "Exchangeable K",
  Na   = "Sodium",
  Al   = "Aluminum",
  BS   = "Base saturation"
)


run_pca_region <- function(df, region_name) {
  
  suppressPackageStartupMessages({
    library(ggrepel)
  })
  
  cat("\n=====================================================\n")
  cat("PCA REGION:", region_name, "\n")
  cat("=====================================================\n")
  
  region_dir_fig <- file.path(figdir, region_name)
  region_dir_tab <- file.path(tabdir, region_name)
  dir.create(region_dir_fig, recursive = TRUE, showWarnings = FALSE)
  dir.create(region_dir_tab, recursive = TRUE, showWarnings = FALSE)
  
  dd <- df %>% filter(Region == region_name)
  

  X <- dd %>% select(all_of(soil_vars))
  

  keep_vars <- sapply(X, function(x) sd(x, na.rm = TRUE) > 0)
  X <- X[, keep_vars, drop = FALSE]
  
  if (ncol(X) < 2) {
    stop("No hay suficientes variables con varianza en ", region_name)
  }
  

  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  

  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  pc1 <- round(var_exp[1] * 100, 2)
  pc2 <- round(var_exp[2] * 100, 2)
  
  cat("Varianza explicada:\n")
  cat("PC1:", pc1, "%\n")
  cat("PC2:", pc2, "%\n\n")
  

  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  colnames(scores)[1:2] <- c("PC1", "PC2")
  scores$Region <- dd$Region
  scores$Treatment <- factor(dd$Treatment, levels = c("Control", "Burned"))
  scores$Years_since_fire <- as.character(dd$Years_since_fire)
  scores$RegionYear <- dd$RegionYear
  

  set.seed(123)
  scores$PC1_plot <- jitter(scores$PC1, amount = 0.08)
  scores$PC2_plot <- jitter(scores$PC2, amount = 0.08)
  

  loadings <- as.data.frame(pca$rotation[, 1:2, drop = FALSE])
  colnames(loadings)[1:2] <- c("PC1", "PC2")
  loadings$Variable <- rownames(loadings)
  rownames(loadings) <- NULL
  loadings$Variable_pretty <- recode(loadings$Variable, !!!pretty_names)
  

  loadings$contrib_PC1 <- (loadings$PC1^2 / sum(loadings$PC1^2)) * 100
  loadings$contrib_PC2 <- (loadings$PC2^2 / sum(loadings$PC2^2)) * 100
  loadings$contrib_total <- loadings$contrib_PC1 + loadings$contrib_PC2
  

  write_csv(scores, file.path(region_dir_tab, paste0("PCA_scores_", region_name, ".csv")))
  write_csv(loadings, file.path(region_dir_tab, paste0("PCA_loadings_", region_name, ".csv")))
  
  variance_tbl <- tibble(
    Region = region_name,
    PC = paste0("PC", seq_along(var_exp)),
    Variance_explained = round(var_exp * 100, 4)
  )
  write_csv(variance_tbl, file.path(region_dir_tab, paste0("PCA_variance_", region_name, ".csv")))
  
  contrib_tbl <- loadings %>%
    arrange(desc(contrib_total))
  write_csv(contrib_tbl, file.path(region_dir_tab, paste0("PCA_variable_contributions_", region_name, ".csv")))
  

  n_arrows <- min(10, nrow(loadings))
  loadings_plot <- loadings %>%
    arrange(desc(contrib_total)) %>%
    slice_head(n = n_arrows)
  

  range_pc1_scores <- max(scores$PC1_plot) - min(scores$PC1_plot)
  range_pc2_scores <- max(scores$PC2_plot) - min(scores$PC2_plot)
  range_pc1_load   <- max(loadings_plot$PC1) - min(loadings_plot$PC1)
  range_pc2_load   <- max(loadings_plot$PC2) - min(loadings_plot$PC2)
  
  if (range_pc1_load == 0) range_pc1_load <- 1
  if (range_pc2_load == 0) range_pc2_load <- 1
  
  arrow_scale <- min(
    range_pc1_scores / range_pc1_load,
    range_pc2_scores / range_pc2_load
  ) * 0.42
  
  loadings_plot <- loadings_plot %>%
    mutate(
      PC1_arrow = PC1 * arrow_scale,
      PC2_arrow = PC2 * arrow_scale
    )
  

  year_levels_present <- unique(scores$Years_since_fire)
  year_levels_present <- year_levels_present[order(as.numeric(year_levels_present))]
  

  year_values <- year_colors[year_levels_present]
  

  p <- ggplot(scores, aes(PC1_plot, PC2_plot)) +
    geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey55") +
    

    stat_ellipse(
      aes(group = Treatment, linetype = Treatment),
      linewidth = 0.8,
      color = "grey55",
      show.legend = FALSE
    ) +
    

    geom_segment(
      data = loadings_plot,
      aes(
        x = 0, y = 0,
        xend = PC1_arrow, yend = PC2_arrow,
        color = contrib_total
      ),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.16, "cm")),
      linewidth = 0.9
    ) +
    
    geom_text_repel(
      data = loadings_plot,
      aes(
        x = PC1_arrow,
        y = PC2_arrow,
        label = Variable_pretty,
        color = contrib_total
      ),
      inherit.aes = FALSE,
      size = 4.3,
      fontface = "bold",
      box.padding = 0.32,
      point.padding = 0.20,
      segment.color = "grey60",
      max.overlaps = 100,
      show.legend = FALSE
    ) +
    

    geom_point(
      aes(fill = Years_since_fire, shape = Treatment),
      size = 5.2,
      color = "black",
      stroke = 0.5,
      alpha = 0.96
    ) +
    
    scale_fill_manual(
      values = year_values,
      breaks = year_levels_present,
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c("Control" = 21, "Burned" = 24)
    ) +
    scale_color_gradient(
      low = "#5E3C99",
      high = "#E66101"
    ) +
    
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(shape = 21, size = 5, color = "black")
      ),
      shape = guide_legend(
        order = 2,
        override.aes = list(fill = "white", size = 5, color = "black")
      ),
      color = guide_colorbar(order = 3)
    ) +
    
    labs(
      title = paste("PCA of soil properties -", region_name),
      x = paste0("PC1 (", pc1, "%)"),
      y = paste0("PC2 (", pc2, "%)"),
      fill = "Years since fire",
      shape = "Treatment",
      color = "Contribution"
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  ggsave(
    file.path(region_dir_fig, paste0("PCA_soil_", region_name, "_final.png")),
    p, width = 9.5, height = 7, dpi = 700, bg = "white"
  )
  
  ggsave(
    file.path(region_dir_fig, paste0("PCA_soil_", region_name, "_final.pdf")),
    p, width = 9.5, height = 7, bg = "white"
  )
  
  return(list(
    scores = scores,
    loadings = loadings,
    variance = variance_tbl,
    contributions = contrib_tbl,
    plot = p
  ))
}


res_amz_pca <- run_pca_region(soil_collapsed, "Amazonas")
res_cus_pca <- run_pca_region(soil_collapsed, "Cusco")


pca_summary <- bind_rows(
  res_amz_pca$variance %>% slice(1:2),
  res_cus_pca$variance %>% slice(1:2)
)

cat("=====================================================\n")
cat("RESUMEN PCA SUELO\n")
cat("=====================================================\n")
print(pca_summary)
cat("\n")

write_csv(pca_summary, file.path(outdir, "PCA_soil_summary_PC1_PC2.csv"))


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(vegan)
  library(tibble)
})


base <- "~/Tesis Rosmery-Metabarcoding/RESULT BACTERIA"
indir <- file.path(base, "results/soil_collapsed")
outdir <- file.path(base, "results/soil_permanova")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

soil_file <- file.path(indir, "soil_collapsed_by_region_year_treatment.csv")
if (!file.exists(soil_file)) {
  stop("No existe: ", soil_file)
}

soil_collapsed <- read_csv(soil_file, show_col_types = FALSE)


soil_vars <- c("pH", "EC", "N", "P", "K", "C", "OM", "CEC", "Ca", "Mg", "K_ex", "Na", "Al", "BS")

missing_vars <- setdiff(soil_vars, colnames(soil_collapsed))
if (length(missing_vars) > 0) {
  stop("Faltan variables en soil_collapsed: ", paste(missing_vars, collapse = ", "))
}

soil_collapsed$Region <- factor(soil_collapsed$Region, levels = c("Amazonas", "Cusco"))
soil_collapsed$Treatment <- factor(soil_collapsed$Treatment, levels = c("Control", "Burned"))
soil_collapsed$Years_since_fire <- as.character(soil_collapsed$Years_since_fire)


run_soil_permanova_region <- function(df, region_name) {
  
  cat("\n=====================================================\n")
  cat("PERMANOVA SUELO - REGION:", region_name, "\n")
  cat("=====================================================\n")
  
  region_dir <- file.path(outdir, region_name)
  dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)
  
  dd <- df %>% filter(Region == region_name)
  

  X <- dd %>% select(all_of(soil_vars))
  

  keep_vars <- sapply(X, function(x) sd(x, na.rm = TRUE) > 0)
  X <- X[, keep_vars, drop = FALSE]
  
  if (ncol(X) < 2) {
    stop("No hay suficientes variables con varianza en ", region_name)
  }
  
  vars_used <- colnames(X)
  

  X_scaled <- scale(X)
  

  dist_soil <- dist(X_scaled, method = "euclidean")
  
  dd$Years_since_fire <- factor(
    dd$Years_since_fire,
    levels = sort(unique(dd$Years_since_fire))
  )
  


  set.seed(123)
  perm_global <- adonis2(
    dist_soil ~ Years_since_fire * Treatment,
    data = dd,
    permutations = 999
  )
  
  cat("PERMANOVA GLOBAL:\n")
  print(perm_global)
  cat("\n")
  
  perm_global_df <- as.data.frame(perm_global)
  perm_global_df$Term <- rownames(perm_global_df)
  rownames(perm_global_df) <- NULL
  
  write_csv(
    perm_global_df,
    file.path(region_dir, paste0("soil_PERMANOVA_global_", region_name, ".csv"))
  )
  


  set.seed(123)
  perm_terms <- adonis2(
    dist_soil ~ Years_since_fire * Treatment,
    data = dd,
    permutations = 999,
    by = "terms"
  )
  
  cat("PERMANOVA POR TERMINOS (by = 'terms'):\n")
  print(perm_terms)
  cat("\n")
  
  perm_terms_df <- as.data.frame(perm_terms)
  perm_terms_df$Term <- rownames(perm_terms_df)
  rownames(perm_terms_df) <- NULL
  
  write_csv(
    perm_terms_df,
    file.path(region_dir, paste0("soil_PERMANOVA_terms_", region_name, ".csv"))
  )
  


  bd_treat <- betadisper(dist_soil, dd$Treatment)
  bd_year  <- betadisper(dist_soil, dd$Years_since_fire)
  
  cat("BETADISPER Treatment - ANOVA:\n")
  print(anova(bd_treat))
  cat("\n")
  
  cat("BETADISPER Treatment - PERMUTEST:\n")
  bd_treat_perm <- permutest(bd_treat, permutations = 999)
  print(bd_treat_perm)
  cat("\n")
  
  cat("BETADISPER Year - ANOVA:\n")
  print(anova(bd_year))
  cat("\n")
  
  cat("BETADISPER Year - PERMUTEST:\n")
  bd_year_perm <- permutest(bd_year, permutations = 999)
  print(bd_year_perm)
  cat("\n")
  
  write_csv(
    data.frame(Term = rownames(as.data.frame(anova(bd_treat))),
               as.data.frame(anova(bd_treat))),
    file.path(region_dir, paste0("soil_BETADISPER_treatment_anova_", region_name, ".csv"))
  )
  
  write_csv(
    data.frame(Term = rownames(as.data.frame(anova(bd_year))),
               as.data.frame(anova(bd_year))),
    file.path(region_dir, paste0("soil_BETADISPER_year_anova_", region_name, ".csv"))
  )
  
  capture.output(
    bd_treat_perm,
    file = file.path(region_dir, paste0("soil_BETADISPER_treatment_permutest_", region_name, ".txt"))
  )
  
  capture.output(
    bd_year_perm,
    file = file.path(region_dir, paste0("soil_BETADISPER_year_permutest_", region_name, ".txt"))
  )
  


  get_val <- function(df_tab, term_name, col_name) {
    x <- df_tab %>% filter(Term == term_name)
    if (nrow(x) == 0) return(NA_real_)
    x[[col_name]][1]
  }
  
  summary_tbl <- tibble(
    Region = region_name,
    n_samples = nrow(dd),
    n_variables = ncol(X),
    Year_R2 = get_val(perm_terms_df, "Years_since_fire", "R2"),
    Year_F = get_val(perm_terms_df, "Years_since_fire", "F"),
    Year_p = get_val(perm_terms_df, "Years_since_fire", "Pr(>F)"),
    Treatment_R2 = get_val(perm_terms_df, "Treatment", "R2"),
    Treatment_F = get_val(perm_terms_df, "Treatment", "F"),
    Treatment_p = get_val(perm_terms_df, "Treatment", "Pr(>F)"),
    Interaction_R2 = get_val(perm_terms_df, "Years_since_fire:Treatment", "R2"),
    Interaction_F = get_val(perm_terms_df, "Years_since_fire:Treatment", "F"),
    Interaction_p = get_val(perm_terms_df, "Years_since_fire:Treatment", "Pr(>F)"),
    Model_R2 = get_val(perm_global_df, "Model", "R2"),
    Model_F = get_val(perm_global_df, "Model", "F"),
    Model_p = get_val(perm_global_df, "Model", "Pr(>F)")
  )
  
  cat("RESUMEN LIMPIO:\n")
  print(summary_tbl)
  cat("\n")
  
  write_csv(
    summary_tbl,
    file.path(region_dir, paste0("soil_PERMANOVA_summary_", region_name, ".csv"))
  )
  
  return(list(
    global = perm_global_df,
    terms = perm_terms_df,
    summary = summary_tbl
  ))
}


res_amz_soil_perm <- run_soil_permanova_region(soil_collapsed, "Amazonas")
res_cus_soil_perm <- run_soil_permanova_region(soil_collapsed, "Cusco")


soil_perm_summary_all <- bind_rows(
  res_amz_soil_perm$summary,
  res_cus_soil_perm$summary
)

cat("=====================================================\n")
cat("RESUMEN GLOBAL PERMANOVA SUELO\n")
cat("=====================================================\n")
print(soil_perm_summary_all)
cat("\n")

write_csv(
  soil_perm_summary_all,
  file.path(outdir, "soil_PERMANOVA_summary_all_regions.csv")
)

cat("Archivos guardados en:\n")
cat(outdir, "\n")


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(broom)
})

base <- "~/Tesis Rosmery-Metabarcoding/RESULT BACTERIA"
pca_dir <- file.path(base, "results/soil_pca/tables")
outdir <- file.path(base, "results/soil_pca_stats")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

run_pca_manova <- function(region_name) {
  
  cat("\n=====================================================\n")
  cat("MANOVA PCA SUELO -", region_name, "\n")
  cat("=====================================================\n")
  
  score_file <- file.path(pca_dir, region_name, paste0("PCA_scores_", region_name, ".csv"))
  if (!file.exists(score_file)) stop("No existe: ", score_file)
  
  df <- read_csv(score_file, show_col_types = FALSE)
  
  df$Treatment <- factor(df$Treatment, levels = c("Control", "Burned"))
  df$Years_since_fire <- factor(df$Years_since_fire,
                                levels = sort(unique(df$Years_since_fire)))
  

  fit_manova <- manova(cbind(PC1, PC2) ~ Years_since_fire + Treatment, data = df)
  
  cat("MANOVA (Pillai):\n")
  manova_pillai <- summary(fit_manova, test = "Pillai")
  print(manova_pillai)
  cat("\n")
  

  cat("MANOVA (Wilks):\n")
  manova_wilks <- summary(fit_manova, test = "Wilks")
  print(manova_wilks)
  cat("\n")
  
  capture.output(
    manova_pillai,
    file = file.path(outdir, paste0("MANOVA_PCA_", region_name, "_Pillai.txt"))
  )
  
  capture.output(
    manova_wilks,
    file = file.path(outdir, paste0("MANOVA_PCA_", region_name, "_Wilks.txt"))
  )
  
  return(list(
    pillai = manova_pillai,
    wilks = manova_wilks
  ))
}

res_amz_manova <- run_pca_manova("Amazonas")
res_cus_manova <- run_pca_manova("Cusco")


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})


base <- "~/Tesis Rosmery-Metabarcoding/RESULT BACTERIA"
indir <- file.path(base, "results/soil_collapsed")
outdir <- file.path(base, "results/soil_univariate_final")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

soil_file <- file.path(indir, "soil_collapsed_by_region_year_treatment.csv")
if (!file.exists(soil_file)) stop("No existe: ", soil_file)

soil <- read_csv(soil_file, show_col_types = FALSE)


soil_vars <- c("pH", "EC", "N", "P", "K", "C", "OM", "CEC", "Ca", "Mg", "K_ex", "Na", "Al", "BS")

pretty_names <- c(
  pH   = "pH",
  EC   = "Electrical conductivity",
  N    = "Nitrogen",
  P    = "Phosphorus",
  K    = "Potassium",
  C    = "Carbon",
  OM   = "Organic matter",
  CEC  = "CEC",
  Ca   = "Calcium",
  Mg   = "Magnesium",
  K_ex = "Exchangeable K",
  Na   = "Sodium",
  Al   = "Aluminum",
  BS   = "Base saturation"
)

soil$Region <- factor(soil$Region, levels = c("Amazonas", "Cusco"))
soil$Treatment <- factor(soil$Treatment, levels = c("Control", "Burned"))
soil$Years_since_fire <- factor(soil$Years_since_fire,
                                levels = sort(unique(as.numeric(soil$Years_since_fire))))


soil_desc <- soil %>%
  pivot_longer(cols = all_of(soil_vars), names_to = "Variable", values_to = "Value") %>%
  group_by(Region, Years_since_fire, Treatment, Variable) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Variable_pretty = recode(Variable, !!!pretty_names))

write_csv(soil_desc, file.path(outdir, "soil_descriptive_collapsed.csv"))


dat <- alpha_df %>%
  left_join(soil_meta, by = "SampleID", suffix = c(".alpha", ".soil"))


if (!"Region" %in% colnames(dat)) {
  if ("Region.alpha" %in% colnames(dat) && "Region.soil" %in% colnames(dat)) {
    dat$Region <- ifelse(!is.na(dat$Region.alpha), dat$Region.alpha, dat$Region.soil)
  } else if ("Region.alpha" %in% colnames(dat)) {
    dat$Region <- dat$Region.alpha
  } else if ("Region.soil" %in% colnames(dat)) {
    dat$Region <- dat$Region.soil
  }
}

if (!"Treatment" %in% colnames(dat)) {
  if ("Treatment.alpha" %in% colnames(dat) && "Treatment.soil" %in% colnames(dat)) {
    dat$Treatment <- ifelse(!is.na(dat$Treatment.alpha), dat$Treatment.alpha, dat$Treatment.soil)
  } else if ("Treatment.alpha" %in% colnames(dat)) {
    dat$Treatment <- dat$Treatment.alpha
  } else if ("Treatment.soil" %in% colnames(dat)) {
    dat$Treatment <- dat$Treatment.soil
  }
}

if (!"Years_since_fire" %in% colnames(dat)) {
  if ("Years_since_fire.alpha" %in% colnames(dat) && "Years_since_fire.soil" %in% colnames(dat)) {
    dat$Years_since_fire <- ifelse(!is.na(dat$Years_since_fire.alpha), dat$Years_since_fire.alpha, dat$Years_since_fire.soil)
  } else if ("Years_since_fire.alpha" %in% colnames(dat)) {
    dat$Years_since_fire <- dat$Years_since_fire.alpha
  } else if ("Years_since_fire.soil" %in% colnames(dat)) {
    dat$Years_since_fire <- dat$Years_since_fire.soil
  }
}

dat$Region <- as.character(dat$Region)
dat$Treatment <- as.character(dat$Treatment)
dat$Years_since_fire <- as.character(dat$Years_since_fire)

write_csv(dat, file.path(tabdir, "alpha_soil_merged.csv"))

cat("Columnas finales de dat:\n")
print(colnames(dat))
cat("\n")


run_alpha_soil_cor <- function(df_region, region_name) {
  
  cat("\n=====================================================\n")
  cat("CORRELACIONES ALFA-SUELO -", region_name, "\n")
  cat("=====================================================\n")
  
  alpha_ok <- alpha_metrics[alpha_metrics %in% colnames(df_region)]
  soil_ok  <- soil_var_names[soil_var_names %in% colnames(df_region)]
  
  cat("Metricas alfa encontradas:\n")
  print(alpha_ok)
  cat("Variables de suelo encontradas:\n")
  print(soil_ok)
  cat("\n")
  
  res <- list()
  
  for (a in alpha_ok) {
    for (s in soil_ok) {
      
      dd <- df_region %>%
        dplyr::select(all_of(a), all_of(s)) %>%
        dplyr::filter(!is.na(.data[[a]]), !is.na(.data[[s]]))
      
      if (nrow(dd) < 4) next
      
      tt <- suppressWarnings(cor.test(dd[[a]], dd[[s]], method = "spearman", exact = FALSE))
      
      res[[paste(a, s, sep = "_")]] <- tibble(
        Region = region_name,
        Alpha_metric = a,
        Soil_variable = s,
        rho = unname(tt$estimate),
        p_value = tt$p.value,
        n = nrow(dd)
      )
    }
  }
  
  out <- bind_rows(res) %>%
    mutate(
      p_adjust_bh = p.adjust(p_value, method = "BH"),
      Significance = case_when(
        p_adjust_bh < 0.001 ~ "***",
        p_adjust_bh < 0.01  ~ "**",
        p_adjust_bh < 0.05  ~ "*",
        p_adjust_bh < 0.1   ~ ".",
        TRUE ~ "ns"
      ),
      Alpha_pretty = recode(Alpha_metric, !!!pretty_alpha),
      Soil_pretty = recode(Soil_variable, !!!pretty_soil)
    ) %>%
    arrange(p_adjust_bh, desc(abs(rho)))
  
  write_csv(out, file.path(tabdir, paste0("alpha_soil_correlations_", region_name, ".csv")))
  
  print(out)
  return(out)
}

res_amz_cor <- run_alpha_soil_cor(dat %>% dplyr::filter(Region == "Amazonas"), "Amazonas")
res_cus_cor <- run_alpha_soil_cor(dat %>% dplyr::filter(Region == "Cusco"), "Cusco")

all_cor <- bind_rows(res_amz_cor, res_cus_cor)
write_csv(all_cor, file.path(tabdir, "alpha_soil_correlations_all_regions.csv"))


library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)

cat("\n===== KRUSKAL-WALLIS / WILCOXON - ALPHA DIVERSITY BY REGION =====\n")
cat("Hora inicio:", as.character(Sys.time()), "\n")


base <- "~/Tesis Rosmery-Metabarcoding/RESULT BACTERIA"


outdir <- file.path(base, "results/alpha_diversity/nonparametric_by_region")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


alpha_use <- alpha_df %>%
  select(
    SampleID,
    Region,
    Treatment,
    Years_since_fire,
    Observed,
    Shannon,
    Simpson
  ) %>%
  pivot_longer(
    cols = c(Observed, Shannon, Simpson),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Region = as.factor(Region),
    Treatment = as.factor(Treatment),
    Years_since_fire = as.factor(Years_since_fire),
    Metric = as.factor(Metric)
  )


cat("\nResumen de datos usados:\n")
print(alpha_use %>% count(Region, Treatment, Years_since_fire, Metric))

write.csv(
  alpha_use,
  file.path(outdir, "alpha_long_observed_shannon_simpson.csv"),
  row.names = FALSE
)


wilcox_treatment_by_region <- alpha_use %>%
  group_by(Region, Metric) %>%
  wilcox_test(Value ~ Treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

print(wilcox_treatment_by_region)

write.csv(
  wilcox_treatment_by_region,
  file.path(outdir, "wilcoxon_treatment_by_region_alpha.csv"),
  row.names = FALSE
)


kruskal_year_by_region <- alpha_use %>%
  group_by(Region, Metric) %>%
  kruskal_test(Value ~ Years_since_fire)

print(kruskal_year_by_region)

write.csv(
  kruskal_year_by_region,
  file.path(outdir, "kruskal_years_by_region_alpha.csv"),
  row.names = FALSE
)


dunn_year_by_region <- alpha_use %>%
  group_by(Region, Metric) %>%
  dunn_test(Value ~ Years_since_fire, p.adjust.method = "BH") %>%
  add_significance("p.adj")

print(dunn_year_by_region)

write.csv(
  dunn_year_by_region,
  file.path(outdir, "dunn_posthoc_years_by_region_alpha.csv"),
  row.names = FALSE
)


wilcox_treatment_by_region_year <- alpha_use %>%
  group_by(Region, Years_since_fire, Metric) %>%
  wilcox_test(Value ~ Treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

print(wilcox_treatment_by_region_year)

write.csv(
  wilcox_treatment_by_region_year,
  file.path(outdir, "wilcoxon_treatment_by_region_year_alpha.csv"),
  row.names = FALSE
)


alpha_control <- alpha_use %>%
  filter(Treatment == "Control")

kruskal_year_control_by_region <- alpha_control %>%
  group_by(Region, Metric) %>%
  kruskal_test(Value ~ Years_since_fire)

print(kruskal_year_control_by_region)

write.csv(
  kruskal_year_control_by_region,
  file.path(outdir, "kruskal_years_control_only_by_region_alpha.csv"),
  row.names = FALSE
)


dunn_year_control_by_region <- alpha_control %>%
  group_by(Region, Metric) %>%
  dunn_test(Value ~ Years_since_fire, p.adjust.method = "BH") %>%
  add_significance("p.adj")

print(dunn_year_control_by_region)

write.csv(
  dunn_year_control_by_region,
  file.path(outdir, "dunn_years_control_only_by_region_alpha.csv"),
  row.names = FALSE
)


p_alpha_region_year <- ggplot(
  alpha_use,
  aes(x = Years_since_fire, y = Value, fill = Treatment)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.85,
    color = "black"
  ) +
  geom_jitter(
    aes(color = Treatment),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 2,
    alpha = 0.8
  ) +
  facet_grid(Metric ~ Region, scales = "free_y") +
  theme_bw(base_size = 13) +
  labs(
    x = "Years since fire",
    y = "Alpha diversity",
    title = "Alpha diversity across post-fire chronosequences by region",
    fill = "Treatment",
    color = "Treatment"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    legend.position = "bottom"
  )

print(p_alpha_region_year)

ggsave(
  file.path(outdir, "alpha_diversity_region_year_treatment.png"),
  p_alpha_region_year,
  width = 12,
  height = 9,
  dpi = 300
)

ggsave(
  file.path(outdir, "alpha_diversity_region_year_treatment.pdf"),
  p_alpha_region_year,
  width = 12,
  height = 9
)


cat("\n===== Non-parametric analyses completed =====\n")
cat("Archivos guardados en:\n")
cat(outdir, "\n")
cat("Hora fin:", as.character(Sys.time()), "\n")
