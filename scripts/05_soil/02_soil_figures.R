# Soil data visualisation
# PCA, soil profiles and correlation visualisations.
# Input: project-specific phyloseq objects, site metadata, or soil data tables.
# Output: analysis tables and figures as defined by the project directory settings.
options(stringsAsFactors = FALSE)

pkgs <- c("ggplot2","ggrepel","patchwork","pheatmap","openxlsx")
miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if(length(miss) > 0){
  stop("Faltan paquetes: ", paste(miss, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(openxlsx)
  library(grid)
})

root <- "/home/usuario1/Rosmery/Tesis Rosmery-Metabarcoding"

soil_file <- file.path(
  root,
  "RESULTADOS_SOIL_UNIQUE_UNITS_FINAL/TABLES/soil_metadata_unique_units.csv"
)

if(!file.exists(soil_file)){
  soil_file <- file.path(
    root,
    "RESULTADOS_SOIL_UNIQUE_UNITS_FINAL/TABLES/soil_metadata_clean.csv"
  )
}

if(!file.exists(soil_file)){
  stop("Unique-unit soil metadata not found.")
}

outdir <- file.path(root, "RESULTADOS_SOIL_FIGURES_Q1_PUBLICATION_FINAL")

dirs <- c(
  outdir,
  file.path(outdir, "01_PDF_EDITABLE"),
  file.path(outdir, "02_PNG_600DPI"),
  file.path(outdir, "03_TABLES"),
  file.path(outdir, "04_EXCEL")
)

invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

soil <- read.csv(soil_file, check.names = FALSE)

if(!"SampleID" %in% names(soil)){
  soil$SampleID <- paste(soil$Region, soil$Treatment, soil$Year, sep = "_")
}

if(!"Year_num" %in% names(soil)){
  soil$Year_num <- suppressWarnings(as.numeric(as.character(soil$Year)))
}

soil$Region <- as.character(soil$Region)
soil$Treatment <- as.character(soil$Treatment)

soil$Treatment_plot <- ifelse(
  grepl("wild|burn|incend|fire|quem", soil$Treatment, ignore.case = TRUE),
  "Burned",
  "Control"
)

soil$Region <- factor(soil$Region, levels = c("Amazonas","Cusco"))
soil$Treatment_plot <- factor(soil$Treatment_plot, levels = c("Control","Burned"))
soil$Year_num <- suppressWarnings(as.numeric(soil$Year_num))

soil_candidates <- c(
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

for(v in soil_candidates){
  if(v %in% names(soil)){
    soil[[v]] <- suppressWarnings(as.numeric(gsub(",", ".", as.character(soil[[v]]))))
  }
}

soil_vars <- soil_candidates[soil_candidates %in% names(soil)]
soil_vars <- soil_vars[sapply(soil[soil_vars], function(x){
  sum(!is.na(x)) >= 6 && sd(x, na.rm = TRUE) > 0
})]

message("Unique soil units: ", nrow(soil))
message("Variables de suelo usadas:")
print(soil_vars)

if(length(soil_vars) < 2){
  stop("Fewer than two valid soil properties.")
}

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

soil <- soil[order(soil$Region, soil$Year_num, soil$Treatment_plot), ]

soil$Sample_label <- paste0(
  ifelse(soil$Region == "Amazonas", "AMZ", "CUS"),
  "-",
  soil$Year_num,
  "-",
  ifelse(soil$Treatment_plot == "Control", "C", "B")
)

write.csv(
  soil,
  file.path(outdir, "03_TABLES", "DATA_soil_unique_units_ordered_used.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(Soil_variable = soil_vars, Label = pretty_var(soil_vars)),
  file.path(outdir, "03_TABLES", "S10_soil_variables_used.csv"),
  row.names = FALSE
)


treat_cols <- c("Control" = "#8ecae6", "Burned" = "#f4a261")
region_cols <- c("Amazonas" = "#FF9E9A", "Cusco" = "#A2D2FF")

year_levels_global <- as.character(sort(unique(soil$Year_num)))
year_shapes <- c(16, 17, 15, 18, 3, 4, 8, 7, 9)
year_shapes <- setNames(year_shapes[seq_along(year_levels_global)], year_levels_global)

theme_paper <- function(){
  theme_bw(base_size = 13) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.45),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black", size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      plot.subtitle = element_text(face = "bold", hjust = 0.5, size = 12),
      strip.text = element_text(face = "bold", size = 11),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10)
    )
}

prepare_matrix <- function(df, vars){
  mat <- as.matrix(df[, vars, drop = FALSE])
  rownames(mat) <- df$SampleID

  for(j in seq_len(ncol(mat))){
    med <- median(mat[,j], na.rm = TRUE)
    if(!is.finite(med)) med <- 0
    mat[is.na(mat[,j]), j] <- med
  }

  matz <- scale(mat)
  matz[!is.finite(matz)] <- 0
  matz
}


run_pca_region <- function(df, region_name){

  df_region <- df[df$Region == region_name, , drop = FALSE]
  matz <- prepare_matrix(df_region, soil_vars)

  pca <- prcomp(matz, center = FALSE, scale. = FALSE)
  expl <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)

  scores <- data.frame(
    SampleID = rownames(pca$x),
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    Region = df_region$Region,
    Treatment = df_region$Treatment_plot,
    Year = factor(as.character(df_region$Year_num), levels = year_levels_global),
    stringsAsFactors = FALSE
  )

  loadings <- data.frame(
    Variable = rownames(pca$rotation),
    PC1 = pca$rotation[,1],
    PC2 = pca$rotation[,2],
    stringsAsFactors = FALSE
  )

  loadings$Label <- pretty_var(loadings$Variable)

  loadings$Contribution <- 100 * rowSums(pca$rotation[,1:2, drop = FALSE]^2) /
    sum(rowSums(pca$rotation[,1:2, drop = FALSE]^2))

  mult <- min(
    diff(range(scores$PC1)) / diff(range(loadings$PC1)),
    diff(range(scores$PC2)) / diff(range(loadings$PC2))
  ) * 0.50

  if(!is.finite(mult)) mult <- 1

  loadings$PC1_plot <- loadings$PC1 * mult
  loadings$PC2_plot <- loadings$PC2 * mult

  write.csv(
    scores,
    file.path(outdir, "03_TABLES", paste0("S12_PCA_scores_", region_name, ".csv")),
    row.names = FALSE
  )

  write.csv(
    loadings[order(-loadings$Contribution), ],
    file.path(outdir, "03_TABLES", paste0("S12_PCA_loadings_contribution_", region_name, ".csv")),
    row.names = FALSE
  )

  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey75", linewidth = 0.45) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey75", linewidth = 0.45) +

    stat_ellipse(
      data = scores,
      aes(x = PC1, y = PC2, linetype = Treatment),
      type = "norm",
      level = 0.68,
      color = "grey55",
      linewidth = 0.8,
      show.legend = FALSE
    ) +

    geom_segment(
      data = loadings,
      aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot, color = Contribution),
      arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"),
      linewidth = 0.75,
      inherit.aes = FALSE
    ) +

    ggrepel::geom_text_repel(
      data = loadings,
      aes(x = PC1_plot, y = PC2_plot, label = Label, color = Contribution),
      size = 3.5,
      fontface = "bold",
      box.padding = 0.35,
      point.padding = 0.15,
      segment.color = "grey60",
      max.overlaps = Inf,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +

    geom_point(
      data = scores,
      aes(x = PC1, y = PC2, fill = Treatment, shape = Year),
      size = 3.8,
      color = "black",
      stroke = 0.55
    ) +

    scale_fill_manual(values = treat_cols, name = "Treatment") +
    scale_shape_manual(values = year_shapes, drop = FALSE, name = "Years since fire") +
    scale_linetype_manual(values = c("Control" = "solid", "Burned" = "dashed")) +
    scale_color_gradientn(
      colors = c("#6A4C93", "#C06C84", "#F4A261"),
      limits = c(0, 12),
      name = "Variable contribution\n(PC1–PC2, %)"
    ) +

    labs(
      title = toupper(region_name),
      subtitle = "(PCA)",
      x = paste0("PC1 (", expl[1], "%)"),
      y = paste0("PC2 (", expl[2], "%)")
    ) +

    guides(
      fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 4, color = "black")),
      shape = guide_legend(order = 2, override.aes = list(size = 4, fill = "grey80", color = "black")),
      color = guide_colorbar(order = 3)
    ) +

    theme_paper() +
    theme(legend.position = "right")

  list(plot = p, scores = scores, loadings = loadings, explained = expl)
}

res_amz <- run_pca_region(soil, "Amazonas")
res_cus <- run_pca_region(soil, "Cusco")

p_pca <- (res_amz$plot + res_cus$plot) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

ggsave(
  file.path(outdir, "02_PNG_600DPI", "S12_Soil_PCA_biplot_combined.png"),
  p_pca, width = 15.5, height = 7.4, dpi = 600
)

ggsave(
  file.path(outdir, "01_PDF_EDITABLE", "S12_Soil_PCA_biplot_combined.pdf"),
  p_pca, width = 15.5, height = 7.4
)


matz <- prepare_matrix(soil, soil_vars)
heat_mat <- t(matz)
rownames(heat_mat) <- pretty_var(rownames(heat_mat))
colnames(heat_mat) <- soil$Sample_label

ann_col <- data.frame(
  Region = soil$Region,
  Treatment = soil$Treatment_plot,
  Year = as.character(soil$Year_num)
)
rownames(ann_col) <- soil$Sample_label

year_unique <- sort(unique(ann_col$Year))
year_cols <- setNames(
  c("#A8DADC", "#B8E0D2", "#F7D6A3", "#F4A261", "#E9AFA3", "#CDB4DB", "#B5E48C", "#90BE6D")[seq_along(year_unique)],
  year_unique
)

ann_colors <- list(
  Region = region_cols,
  Treatment = treat_cols,
  Year = year_cols
)

group_key <- paste(soil$Region, soil$Year_num, sep = "_")
gaps_col <- which(group_key[-length(group_key)] != group_key[-1])
gaps_col <- unique(c(gaps_col, which(soil$Region[-length(soil$Region)] != soil$Region[-1])))

soil_hm_cols <- colorRampPalette(c("#B8B8FF", "#FFF8E7", "#FFB4A2"))(100)

png(
  file.path(outdir, "02_PNG_600DPI", "S13_Soil_scaled_heatmap_ordered.png"),
  width = 3900, height = 2600, res = 300
)
pheatmap::pheatmap(
  heat_mat,
  color = soil_hm_cols,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  border_color = "grey82",
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 9,
  angle_col = 45,
  main = "Scaled soil physicochemical variables"
)
dev.off()

pdf(
  file.path(outdir, "01_PDF_EDITABLE", "S13_Soil_scaled_heatmap_ordered.pdf"),
  width = 13, height = 8.5
)
pheatmap::pheatmap(
  heat_mat,
  color = soil_hm_cols,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  border_color = "grey82",
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 9,
  angle_col = 45,
  main = "Scaled soil physicochemical variables"
)
dev.off()

write.csv(
  data.frame(Sample_label = soil$Sample_label, soil[, c("SampleID","Region","Treatment_plot","Year_num")]),
  file.path(outdir, "03_TABLES", "S13_heatmap_sample_order.csv"),
  row.names = FALSE
)


soil_cor <- cor(
  as.matrix(soil[, soil_vars, drop = FALSE]),
  method = "spearman",
  use = "pairwise.complete.obs"
)

rownames(soil_cor) <- pretty_var(rownames(soil_cor))
colnames(soil_cor) <- pretty_var(colnames(soil_cor))

cor_cols <- colorRampPalette(c("#A8DADC", "#FFF8E7", "#E9C46A"))(100)

png(
  file.path(outdir, "02_PNG_600DPI", "S14_Soil_Spearman_correlation_heatmap.png"),
  width = 2800, height = 2500, res = 300
)
pheatmap::pheatmap(
  soil_cor,
  color = cor_cols,
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize = 10,
  fontsize_number = 7,
  main = "Spearman correlation among soil physicochemical variables"
)
dev.off()

pdf(
  file.path(outdir, "01_PDF_EDITABLE", "S14_Soil_Spearman_correlation_heatmap.pdf"),
  width = 9.5, height = 9
)
pheatmap::pheatmap(
  soil_cor,
  color = cor_cols,
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize = 10,
  fontsize_number = 7,
  main = "Spearman correlation among soil physicochemical variables"
)
dev.off()

write.csv(
  soil_cor,
  file.path(outdir, "03_TABLES", "S14_Soil_Spearman_correlation_matrix.csv")
)


all_loadings <- rbind(
  data.frame(Region = "Amazonas", res_amz$loadings),
  data.frame(Region = "Cusco", res_cus$loadings)
)

avg_contrib <- aggregate(
  Contribution ~ Variable + Label,
  data = all_loadings,
  FUN = mean
)

avg_contrib <- avg_contrib[order(-avg_contrib$Contribution), ]
top_vars <- head(avg_contrib$Variable, 8)

write.csv(
  avg_contrib,
  file.path(outdir, "03_TABLES", "S12_average_PCA_variable_contribution.csv"),
  row.names = FALSE
)

make_profile_plot <- function(vars_use, filename_prefix, title_text, fig_height){

  soil_long <- reshape(
    soil[, c("SampleID","Region","Treatment_plot","Year_num", vars_use)],
    varying = vars_use,
    v.names = "Value",
    timevar = "Soil_variable",
    times = vars_use,
    direction = "long"
  )

  rownames(soil_long) <- NULL
  soil_long$Soil_variable_label <- pretty_var(soil_long$Soil_variable)
  soil_long$Soil_variable_label <- factor(soil_long$Soil_variable_label, levels = pretty_var(vars_use))
  soil_long$YearF_region <- factor(soil_long$Year_num)

  pd <- position_dodge(width = 0.35)

  p <- ggplot(
    soil_long,
    aes(x = YearF_region, y = Value, color = Treatment_plot, group = Treatment_plot)
  ) +
    geom_line(position = pd, linewidth = 0.5, alpha = 0.8) +
    geom_point(
      aes(shape = Treatment_plot),
      position = pd,
      size = 2.5,
      stroke = 0.8
    ) +
    facet_grid(
      Soil_variable_label ~ Region,
      scales = "free",
      space = "free_x"
    ) +
    scale_color_manual(values = treat_cols, name = "Treatment") +
    scale_shape_manual(values = c("Control" = 16, "Burned" = 17), name = "Treatment") +
    labs(
      title = title_text,
      x = "Years since fire",
      y = "Soil variable value"
    ) +
    theme_paper() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      strip.text.y = element_text(angle = 0, hjust = 0, size = 9),
      legend.position = "right"
    )

  ggsave(
    file.path(outdir, "02_PNG_600DPI", paste0(filename_prefix, ".png")),
    p,
    width = 14,
    height = fig_height,
    dpi = 600
  )

  ggsave(
    file.path(outdir, "01_PDF_EDITABLE", paste0(filename_prefix, ".pdf")),
    p,
    width = 14,
    height = fig_height
  )

  p
}

p_profiles_key <- make_profile_plot(
  top_vars,
  "S15A_Key_soil_physicochemical_profiles",
  "Key soil physicochemical profiles across post-fire chronosequences",
  max(7, length(top_vars) * 0.85)
)

p_profiles_all <- make_profile_plot(
  soil_vars,
  "S15B_Extended_all_soil_physicochemical_profiles",
  "Extended soil physicochemical profiles across post-fire chronosequences",
  max(10, length(soil_vars) * 0.72)
)


wb <- openxlsx::createWorkbook()

add_sheet <- function(name, df){
  name <- substr(name, 1, 31)
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeData(wb, name, df)
  openxlsx::freezePane(wb, name, firstRow = TRUE)
  openxlsx::setColWidths(wb, name, cols = 1:ncol(df), widths = "auto")
}

add_sheet("Soil_units_ordered", soil)
add_sheet("Soil_variables", data.frame(Soil_variable = soil_vars, Label = pretty_var(soil_vars)))
add_sheet("PCA_scores_Amazonas", res_amz$scores)
add_sheet("PCA_loadings_Amazonas", res_amz$loadings)
add_sheet("PCA_scores_Cusco", res_cus$scores)
add_sheet("PCA_loadings_Cusco", res_cus$loadings)
add_sheet("Average_PCA_contribution", avg_contrib)
add_sheet("Spearman_correlation", as.data.frame(soil_cor))

openxlsx::saveWorkbook(
  wb,
  file.path(outdir, "04_EXCEL", "SOIL_FIGURES_Q1_PUBLICATION_FINAL_RESULTS.xlsx"),
  overwrite = TRUE
)

message("\nDONE.")
message("Output folder: ", outdir)
message("Main PCA PDF: ", file.path(outdir, "01_PDF_EDITABLE", "S12_Soil_PCA_biplot_combined.pdf"))
message("Ordered heatmap PDF: ", file.path(outdir, "01_PDF_EDITABLE", "S13_Soil_scaled_heatmap_ordered.pdf"))
message("Correlation PDF: ", file.path(outdir, "01_PDF_EDITABLE", "S14_Soil_Spearman_correlation_heatmap.pdf"))
message("Key profiles PDF: ", file.path(outdir, "01_PDF_EDITABLE", "S15A_Key_soil_physicochemical_profiles.pdf"))
message("Extended profiles PDF: ", file.path(outdir, "01_PDF_EDITABLE", "S15B_Extended_all_soil_physicochemical_profiles.pdf"))
