# Alpha diversity
# Genus richness and Hill diversity inference with model assumption checks.
# Input: project-specific phyloseq objects, site metadata, or soil data tables.
# Output: analysis tables and figures as defined by the project directory settings.


rm(list = ls())


packages_cran <- c(
  "dplyr", "tidyr", "ggplot2", "patchwork",
  "openxlsx", "car", "emmeans", "ARTool",
  "broom", "zip", "svglite"
)

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(openxlsx)
library(car)
library(emmeans)
library(ARTool)
library(broom)
library(zip)
library(svglite)

cat("\n=====================================================\n")
cat("ALPHA DIVERSITY AT GENUS LEVEL\n")
cat("D0 q=0 | D1 q=1 | D2 q=2\n")
cat("=====================================================\n")


project_root <- path.expand("~/Rosmery/Tesis Rosmery-Metabarcoding")

genus_supp_dir <- file.path(
  project_root,
  "RESULTADOS_ALPHA_DIVERSITY_GENUS_SUPUESTOS"
)

alpha_long_file <- file.path(
  genus_supp_dir,
  "genus_alpha_values",
  "alpha_diversity_GENUS_D0_D1_D2_values_long.csv"
)

alpha_wide_file <- file.path(
  genus_supp_dir,
  "genus_alpha_values",
  "alpha_diversity_GENUS_D0_D1_D2_values_wide.csv"
)

decision_file <- file.path(
  genus_supp_dir,
  "tables",
  "08_DECISION_TABLE_which_test_to_apply_GENUS.csv"
)

genus_qc_marker_file <- file.path(
  genus_supp_dir,
  "genus_alpha_values",
  "genus_level_QC_by_marker.csv"
)

genus_qc_sample_file <- file.path(
  genus_supp_dir,
  "genus_alpha_values",
  "genus_level_QC_by_sample.csv"
)

if (!file.exists(alpha_long_file)) {
  stop("Required alpha-diversity values file not found: alpha_diversity_GENUS_D0_D1_D2_values_long.csv")
}

if (!file.exists(decision_file)) {
  stop("Required statistical decision table not found: 08_DECISION_TABLE_which_test_to_apply_GENUS.csv")
}

outdir <- file.path(
  project_root,
  "RESULTADOS_ALPHA_DIVERSITY_GENUS_FINAL_PAPER"
)

fig_dir <- file.path(outdir, "01_FIGURES")
tab_dir <- file.path(outdir, "02_TABLES_CSV")
excel_dir <- file.path(outdir, "03_EXCEL_TABLES")
supp_dir <- file.path(outdir, "04_SUPPLEMENTARY")
text_dir <- file.path(outdir, "05_TEXT_FOR_MANUSCRIPT")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(excel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)


alpha_long <- read.csv(alpha_long_file, stringsAsFactors = FALSE)
decision_table <- read.csv(decision_file, stringsAsFactors = FALSE)

alpha_long <- alpha_long %>%
  mutate(
    Marker = as.character(Marker),
    Prefix = as.character(Prefix),
    SampleID = as.character(SampleID),
    SampleLabel = as.character(SampleLabel),
    Region = as.character(Region),
    Treatment = as.character(Treatment),
    Year = as.character(Year),
    Year_num = as.numeric(Year_num),
    Metric = as.character(Metric),
    Metric_label = as.character(Metric_label),
    Value = as.numeric(Value)
  ) %>%
  filter(
    Prefix %in% c("BAC", "ARQ", "ITS"),
    Region %in% c("Amazonas", "Cusco"),
    Treatment %in% c("Control", "Wildfire"),
    Metric %in% c(
      "D0_q0_Genus_Richness",
      "D1_q1_Genus_Shannon_Hill",
      "D2_q2_Genus_Simpson_Hill"
    ),
    is.finite(Value)
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Control", "Wildfire")),
    Region = factor(Region, levels = c("Amazonas", "Cusco")),
    Prefix = factor(Prefix, levels = c("BAC", "ARQ", "ITS")),
    Metric_label = factor(
      Metric_label,
      levels = c(
        "Genus richness D0 (q=0)",
        "Genus Shannon diversity D1 (q=1)",
        "Genus Simpson diversity D2 (q=2)"
      )
    )
  )

decision_table <- decision_table %>%
  mutate(
    Marker = as.character(Marker),
    Prefix = as.character(Prefix),
    Region = as.character(Region),
    Metric = as.character(Metric),
    Metric_label = as.character(Metric_label),
    Suggested_global_test = as.character(Suggested_global_test),
    Suggested_pairwise_test = as.character(Suggested_pairwise_test),
    Final_test_family = case_when(
      grepl("ANOVA factorial", Suggested_global_test, ignore.case = TRUE) ~ "Parametric_ANOVA",
      grepl("Welch", Suggested_global_test, ignore.case = TRUE) ~ "Welch_or_robust",
      grepl("ART|Kruskal|no param", Suggested_global_test, ignore.case = TRUE) ~ "Nonparametric_ART",
      TRUE ~ "Review"
    )
  )


marker_names <- c(
  "BAC" = "Bacteria 16S",
  "ARQ" = "Archaea 16S",
  "ITS" = "Fungi ITS"
)

metric_titles <- c(
  "Genus richness D0 (q=0)" = "Genus richness D0 (q=0)",
  "Genus Shannon diversity D1 (q=1)" = "Genus Shannon diversity D1 (q=1)",
  "Genus Simpson diversity D2 (q=2)" = "Genus Simpson diversity D2 (q=2)"
)

treatment_colors <- c(
  "Control" = "#9BD7E8",
  "Wildfire" = "#F4A261"
)

p_format <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.4f", p)
  )
}

p_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

round2 <- function(x) {
  ifelse(is.na(x), NA, round(x, 2))
}

safe_percent_change <- function(wild, control) {
  ifelse(
    is.na(wild) | is.na(control) | control == 0,
    NA_real_,
    ((wild - control) / control) * 100
  )
}


descriptive_full <- alpha_long %>%
  group_by(
    Marker,
    Prefix,
    Region,
    Year,
    Year_num,
    Treatment,
    Metric,
    Metric_label
  ) %>%
  summarise(
    n = n(),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    se = sd / sqrt(n),
    median = median(Value, na.rm = TRUE),
    q1 = quantile(Value, 0.25, na.rm = TRUE),
    q3 = quantile(Value, 0.75, na.rm = TRUE),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    mean_sd = paste0(sprintf("%.2f", mean), " ± ", sprintf("%.2f", sd)),
    median_iqr = paste0(
      sprintf("%.2f", median),
      " [",
      sprintf("%.2f", q1),
      "–",
      sprintf("%.2f", q3),
      "]"
    ),
    .groups = "drop"
  ) %>%
  arrange(Prefix, Region, Year_num, Treatment, Metric_label)

descriptive_paper <- descriptive_full %>%
  transmute(
    Marker,
    Region,
    Year,
    Treatment,
    Metric = Metric_label,
    n,
    `Mean ± SD` = mean_sd,
    `Median [IQR]` = median_iqr,
    Min = round2(min),
    Max = round2(max)
  )


run_parametric_anova <- function(df) {
  
  df <- df %>%
    filter(is.finite(Value)) %>%
    mutate(
      Treatment = factor(Treatment, levels = c("Control", "Wildfire")),
      Year = factor(Year, levels = as.character(sort(unique(as.numeric(as.character(Year))))))
    )
  
  marker <- unique(df$Marker)[1]
  prefix <- unique(as.character(df$Prefix))[1]
  region <- unique(as.character(df$Region))[1]
  metric <- unique(df$Metric)[1]
  metric_label <- unique(as.character(df$Metric_label))[1]
  
  old_contrasts <- options("contrasts")
  options(contrasts = c("contr.sum", "contr.poly"))
  
  fit <- tryCatch(
    lm(Value ~ Treatment * Year, data = df),
    error = function(e) NULL
  )
  
  options(contrasts = old_contrasts$contrasts)
  
  if (is.null(fit)) {
    return(
      data.frame(
        Marker = marker,
        Prefix = prefix,
        Region = region,
        Metric = metric,
        Metric_label = metric_label,
        Test_family = "Parametric_ANOVA",
        Global_test = "Type III factorial ANOVA",
        Term = NA_character_,
        df1 = NA_real_,
        df2 = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        partial_eta2 = NA_real_,
        status = "Error fitting model",
        stringsAsFactors = FALSE
      )
    )
  }
  
  a3 <- tryCatch(
    car::Anova(fit, type = 3),
    error = function(e) NULL
  )
  
  if (is.null(a3)) {
    return(
      data.frame(
        Marker = marker,
        Prefix = prefix,
        Region = region,
        Metric = metric,
        Metric_label = metric_label,
        Test_family = "Parametric_ANOVA",
        Global_test = "Type III factorial ANOVA",
        Term = NA_character_,
        df1 = NA_real_,
        df2 = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        partial_eta2 = NA_real_,
        status = "Error in ANOVA",
        stringsAsFactors = FALSE
      )
    )
  }
  
  a3_df <- as.data.frame(a3)
  a3_df$Term <- rownames(a3_df)
  rownames(a3_df) <- NULL
  
  ss_col <- grep("Sum Sq|Sum.Sq", colnames(a3_df), value = TRUE)[1]
  df_col <- grep("^Df$", colnames(a3_df), value = TRUE)[1]
  f_col <- grep("F value|F.value|F$", colnames(a3_df), value = TRUE)[1]
  p_col <- grep("Pr\\(>F\\)|Pr", colnames(a3_df), value = TRUE)[1]
  
  ss_resid <- a3_df[[ss_col]][a3_df$Term == "Residuals"]
  df_resid <- a3_df[[df_col]][a3_df$Term == "Residuals"]
  
  out <- a3_df %>%
    filter(!Term %in% c("(Intercept)", "Residuals")) %>%
    transmute(
      Marker = marker,
      Prefix = prefix,
      Region = region,
      Metric = metric,
      Metric_label = metric_label,
      Test_family = "Parametric_ANOVA",
      Global_test = "Type III factorial ANOVA",
      Term = Term,
      df1 = .data[[df_col]],
      df2 = df_resid,
      statistic = .data[[f_col]],
      p_value = .data[[p_col]],
      partial_eta2 = ifelse(
        is.na(ss_resid) | ss_resid == 0,
        NA_real_,
        .data[[ss_col]] / (.data[[ss_col]] + ss_resid)
      ),
      status = "OK"
    )
  
  return(out)
}

run_nonparametric_art <- function(df) {
  
  df <- df %>%
    filter(is.finite(Value)) %>%
    mutate(
      Treatment = factor(Treatment, levels = c("Control", "Wildfire")),
      Year = factor(Year, levels = as.character(sort(unique(as.numeric(as.character(Year))))))
    )
  
  marker <- unique(df$Marker)[1]
  prefix <- unique(as.character(df$Prefix))[1]
  region <- unique(as.character(df$Region))[1]
  metric <- unique(df$Metric)[1]
  metric_label <- unique(as.character(df$Metric_label))[1]
  
  art_fit <- tryCatch(
    ARTool::art(Value ~ Treatment * Year, data = df),
    error = function(e) NULL
  )
  
  if (is.null(art_fit)) {
    return(
      data.frame(
        Marker = marker,
        Prefix = prefix,
        Region = region,
        Metric = metric,
        Metric_label = metric_label,
        Test_family = "Nonparametric_ART",
        Global_test = "Aligned rank transform ANOVA",
        Term = NA_character_,
        df1 = NA_real_,
        df2 = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        partial_eta2 = NA_real_,
        status = "Error fitting ART model",
        stringsAsFactors = FALSE
      )
    )
  }
  
  art_aov <- tryCatch(
    anova(art_fit),
    error = function(e) NULL
  )
  
  if (is.null(art_aov)) {
    return(
      data.frame(
        Marker = marker,
        Prefix = prefix,
        Region = region,
        Metric = metric,
        Metric_label = metric_label,
        Test_family = "Nonparametric_ART",
        Global_test = "Aligned rank transform ANOVA",
        Term = NA_character_,
        df1 = NA_real_,
        df2 = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        partial_eta2 = NA_real_,
        status = "Error in ART ANOVA",
        stringsAsFactors = FALSE
      )
    )
  }
  
  art_df <- as.data.frame(art_aov)
  art_df$Term <- rownames(art_df)
  rownames(art_df) <- NULL
  
  df1_col <- grep("^Df$|Df$", colnames(art_df), value = TRUE)[1]
  df2_col <- grep("Df.res|Df\\.res|Den Df|DenDF", colnames(art_df), value = TRUE)[1]
  f_col <- grep("^F$|F value|F.value", colnames(art_df), value = TRUE)[1]
  p_col <- grep("Pr\\(>F\\)|Pr", colnames(art_df), value = TRUE)[1]
  
  if (is.na(df2_col)) {
    art_df$df2_tmp <- NA_real_
    df2_col <- "df2_tmp"
  }
  
  out <- art_df %>%
    filter(Term %in% c("Treatment", "Year", "Treatment:Year")) %>%
    transmute(
      Marker = marker,
      Prefix = prefix,
      Region = region,
      Metric = metric,
      Metric_label = metric_label,
      Test_family = "Nonparametric_ART",
      Global_test = "Aligned rank transform ANOVA",
      Term = Term,
      df1 = .data[[df1_col]],
      df2 = .data[[df2_col]],
      statistic = .data[[f_col]],
      p_value = .data[[p_col]],
      partial_eta2 = NA_real_,
      status = "OK"
    )
  
  return(out)
}

run_welch_global <- function(df) {
  
  df <- df %>%
    filter(is.finite(Value)) %>%
    mutate(
      Treatment = factor(Treatment, levels = c("Control", "Wildfire")),
      Year = factor(Year, levels = as.character(sort(unique(as.numeric(as.character(Year)))))),
      Group = interaction(Treatment, Year, drop = TRUE)
    )
  
  marker <- unique(df$Marker)[1]
  prefix <- unique(as.character(df$Prefix))[1]
  region <- unique(as.character(df$Region))[1]
  metric <- unique(df$Metric)[1]
  metric_label <- unique(as.character(df$Metric_label))[1]
  
  wt <- tryCatch(
    oneway.test(Value ~ Group, data = df, var.equal = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(wt)) {
    return(
      data.frame(
        Marker = marker,
        Prefix = prefix,
        Region = region,
        Metric = metric,
        Metric_label = metric_label,
        Test_family = "Welch_or_robust",
        Global_test = "Welch one-way test on Treatment-Year groups",
        Term = "Treatment-Year groups",
        df1 = NA_real_,
        df2 = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        partial_eta2 = NA_real_,
        status = "Error in Welch test",
        stringsAsFactors = FALSE
      )
    )
  }
  
  data.frame(
    Marker = marker,
    Prefix = prefix,
    Region = region,
    Metric = metric,
    Metric_label = metric_label,
    Test_family = "Welch_or_robust",
    Global_test = "Welch one-way test on Treatment-Year groups",
    Term = "Treatment-Year groups",
    df1 = unname(wt$parameter[1]),
    df2 = unname(wt$parameter[2]),
    statistic = unname(wt$statistic),
    p_value = wt$p.value,
    partial_eta2 = NA_real_,
    status = "OK",
    stringsAsFactors = FALSE
  )
}


run_pairwise_emmeans <- function(df) {
  
  df <- df %>%
    filter(is.finite(Value)) %>%
    mutate(
      Treatment = factor(Treatment, levels = c("Control", "Wildfire")),
      Year = factor(Year, levels = as.character(sort(unique(as.numeric(as.character(Year))))))
    )
  
  marker <- unique(df$Marker)[1]
  prefix <- unique(as.character(df$Prefix))[1]
  region <- unique(as.character(df$Region))[1]
  metric <- unique(df$Metric)[1]
  metric_label <- unique(as.character(df$Metric_label))[1]
  
  old_contrasts <- options("contrasts")
  options(contrasts = c("contr.sum", "contr.poly"))
  
  fit <- tryCatch(
    lm(Value ~ Treatment * Year, data = df),
    error = function(e) NULL
  )
  
  options(contrasts = old_contrasts$contrasts)
  
  if (is.null(fit)) return(data.frame())
  
  emm <- tryCatch(
    emmeans::emmeans(fit, ~ Treatment | Year),
    error = function(e) NULL
  )
  
  if (is.null(emm)) return(data.frame())
  
  cont <- as.data.frame(
    summary(
      pairs(emm, reverse = FALSE),
      infer = TRUE,
      adjust = "none"
    )
  )
  
  means_tab <- df %>%
    group_by(Year, Treatment) %>%
    summarise(
      n = n(),
      mean = mean(Value, na.rm = TRUE),
      sd = sd(Value, na.rm = TRUE),
      median = median(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Treatment,
      values_from = c(n, mean, sd, median),
      names_sep = "_"
    )
  
  out <- cont %>%
    mutate(
      Year = as.character(Year),
      Marker = marker,
      Prefix = prefix,
      Region = region,
      Metric = metric,
      Metric_label = metric_label,
      Pairwise_test = "EMMeans contrast after ANOVA",
      comparison = as.character(contrast),
      statistic = t.ratio,
      p_value = p.value
    ) %>%
    left_join(means_tab, by = "Year") %>%
    transmute(
      Marker,
      Prefix,
      Region,
      Metric,
      Metric_label,
      Year,
      Pairwise_test,
      comparison,
      n_control = n_Control,
      n_wildfire = n_Wildfire,
      mean_control = mean_Control,
      mean_wildfire = mean_Wildfire,
      sd_control = sd_Control,
      sd_wildfire = sd_Wildfire,
      median_control = median_Control,
      median_wildfire = median_Wildfire,
      estimate_control_minus_wildfire = estimate,
      percent_change_wildfire_vs_control = safe_percent_change(mean_Wildfire, mean_Control),
      SE = SE,
      df = .data$df,
      statistic = statistic,
      p_value = p_value,
      status = "OK"
    )
  
  return(out)
}

run_pairwise_welch <- function(df) {
  
  df <- df %>%
    filter(is.finite(Value)) %>%
    mutate(
      Treatment = factor(Treatment, levels = c("Control", "Wildfire")),
      Year = factor(Year, levels = as.character(sort(unique(as.numeric(as.character(Year))))))
    )
  
  marker <- unique(df$Marker)[1]
  prefix <- unique(as.character(df$Prefix))[1]
  region <- unique(as.character(df$Region))[1]
  metric <- unique(df$Metric)[1]
  metric_label <- unique(as.character(df$Metric_label))[1]
  
  out <- df %>%
    group_by(Year) %>%
    group_modify(
      function(dat, key) {
        
        n_control <- sum(dat$Treatment == "Control")
        n_wildfire <- sum(dat$Treatment == "Wildfire")
        
        mean_control <- mean(dat$Value[dat$Treatment == "Control"], na.rm = TRUE)
        mean_wildfire <- mean(dat$Value[dat$Treatment == "Wildfire"], na.rm = TRUE)
        sd_control <- sd(dat$Value[dat$Treatment == "Control"], na.rm = TRUE)
        sd_wildfire <- sd(dat$Value[dat$Treatment == "Wildfire"], na.rm = TRUE)
        median_control <- median(dat$Value[dat$Treatment == "Control"], na.rm = TRUE)
        median_wildfire <- median(dat$Value[dat$Treatment == "Wildfire"], na.rm = TRUE)
        
        tt <- tryCatch(
          t.test(Value ~ Treatment, data = dat, var.equal = FALSE),
          error = function(e) NULL
        )
        
        if (is.null(tt)) {
          return(
            data.frame(
              Pairwise_test = "Welch t-test",
              comparison = "Control vs Wildfire",
              n_control = n_control,
              n_wildfire = n_wildfire,
              mean_control = mean_control,
              mean_wildfire = mean_wildfire,
              sd_control = sd_control,
              sd_wildfire = sd_wildfire,
              median_control = median_control,
              median_wildfire = median_wildfire,
              estimate_control_minus_wildfire = NA_real_,
              percent_change_wildfire_vs_control = safe_percent_change(mean_wildfire, mean_control),
              SE = NA_real_,
              df = NA_real_,
              statistic = NA_real_,
              p_value = NA_real_,
              status = "Error"
            )
          )
        }
        
        data.frame(
          Pairwise_test = "Welch t-test",
          comparison = "Control vs Wildfire",
          n_control = n_control,
          n_wildfire = n_wildfire,
          mean_control = mean_control,
          mean_wildfire = mean_wildfire,
          sd_control = sd_control,
          sd_wildfire = sd_wildfire,
          median_control = median_control,
          median_wildfire = median_wildfire,
          estimate_control_minus_wildfire = mean_control - mean_wildfire,
          percent_change_wildfire_vs_control = safe_percent_change(mean_wildfire, mean_control),
          SE = NA_real_,
          df = unname(tt$parameter),
          statistic = unname(tt$statistic),
          p_value = tt$p.value,
          status = "OK"
        )
      }
    ) %>%
    ungroup() %>%
    mutate(
      Marker = marker,
      Prefix = prefix,
      Region = region,
      Metric = metric,
      Metric_label = metric_label,
      Year = as.character(Year)
    ) %>%
    select(
      Marker,
      Prefix,
      Region,
      Metric,
      Metric_label,
      Year,
      Pairwise_test,
      comparison,
      n_control,
      n_wildfire,
      mean_control,
      mean_wildfire,
      sd_control,
      sd_wildfire,
      median_control,
      median_wildfire,
      estimate_control_minus_wildfire,
      percent_change_wildfire_vs_control,
      SE,
      df,
      statistic,
      p_value,
      status
    )
  
  return(out)
}

run_pairwise_wilcox <- function(df) {
  
  df <- df %>%
    filter(is.finite(Value)) %>%
    mutate(
      Treatment = factor(Treatment, levels = c("Control", "Wildfire")),
      Year = factor(Year, levels = as.character(sort(unique(as.numeric(as.character(Year))))))
    )
  
  marker <- unique(df$Marker)[1]
  prefix <- unique(as.character(df$Prefix))[1]
  region <- unique(as.character(df$Region))[1]
  metric <- unique(df$Metric)[1]
  metric_label <- unique(as.character(df$Metric_label))[1]
  
  out <- df %>%
    group_by(Year) %>%
    group_modify(
      function(dat, key) {
        
        n_control <- sum(dat$Treatment == "Control")
        n_wildfire <- sum(dat$Treatment == "Wildfire")
        
        mean_control <- mean(dat$Value[dat$Treatment == "Control"], na.rm = TRUE)
        mean_wildfire <- mean(dat$Value[dat$Treatment == "Wildfire"], na.rm = TRUE)
        sd_control <- sd(dat$Value[dat$Treatment == "Control"], na.rm = TRUE)
        sd_wildfire <- sd(dat$Value[dat$Treatment == "Wildfire"], na.rm = TRUE)
        median_control <- median(dat$Value[dat$Treatment == "Control"], na.rm = TRUE)
        median_wildfire <- median(dat$Value[dat$Treatment == "Wildfire"], na.rm = TRUE)
        
        wt <- tryCatch(
          wilcox.test(Value ~ Treatment, data = dat, exact = FALSE),
          error = function(e) NULL
        )
        
        if (is.null(wt)) {
          return(
            data.frame(
              Pairwise_test = "Wilcoxon rank-sum",
              comparison = "Control vs Wildfire",
              n_control = n_control,
              n_wildfire = n_wildfire,
              mean_control = mean_control,
              mean_wildfire = mean_wildfire,
              sd_control = sd_control,
              sd_wildfire = sd_wildfire,
              median_control = median_control,
              median_wildfire = median_wildfire,
              estimate_control_minus_wildfire = NA_real_,
              percent_change_wildfire_vs_control = safe_percent_change(mean_wildfire, mean_control),
              SE = NA_real_,
              df = NA_real_,
              statistic = NA_real_,
              p_value = NA_real_,
              status = "Error"
            )
          )
        }
        
        data.frame(
          Pairwise_test = "Wilcoxon rank-sum",
          comparison = "Control vs Wildfire",
          n_control = n_control,
          n_wildfire = n_wildfire,
          mean_control = mean_control,
          mean_wildfire = mean_wildfire,
          sd_control = sd_control,
          sd_wildfire = sd_wildfire,
          median_control = median_control,
          median_wildfire = median_wildfire,
          estimate_control_minus_wildfire = mean_control - mean_wildfire,
          percent_change_wildfire_vs_control = safe_percent_change(mean_wildfire, mean_control),
          SE = NA_real_,
          df = NA_real_,
          statistic = unname(wt$statistic),
          p_value = wt$p.value,
          status = "OK"
        )
      }
    ) %>%
    ungroup() %>%
    mutate(
      Marker = marker,
      Prefix = prefix,
      Region = region,
      Metric = metric,
      Metric_label = metric_label,
      Year = as.character(Year)
    ) %>%
    select(
      Marker,
      Prefix,
      Region,
      Metric,
      Metric_label,
      Year,
      Pairwise_test,
      comparison,
      n_control,
      n_wildfire,
      mean_control,
      mean_wildfire,
      sd_control,
      sd_wildfire,
      median_control,
      median_wildfire,
      estimate_control_minus_wildfire,
      percent_change_wildfire_vs_control,
      SE,
      df,
      statistic,
      p_value,
      status
    )
  
  return(out)
}


global_list <- list()
pairwise_list <- list()

for (i in seq_len(nrow(decision_table))) {
  
  dec_i <- decision_table[i, ]
  
  df_i <- alpha_long %>%
    filter(
      as.character(Prefix) == dec_i$Prefix,
      as.character(Region) == dec_i$Region,
      Metric == dec_i$Metric
    )
  
  if (nrow(df_i) == 0) next
  
  if (dec_i$Final_test_family == "Parametric_ANOVA") {
    
    global_list[[length(global_list) + 1]] <- run_parametric_anova(df_i)
    pairwise_list[[length(pairwise_list) + 1]] <- run_pairwise_emmeans(df_i)
    
  } else if (dec_i$Final_test_family == "Welch_or_robust") {
    
    global_list[[length(global_list) + 1]] <- run_welch_global(df_i)
    pairwise_list[[length(pairwise_list) + 1]] <- run_pairwise_welch(df_i)
    
  } else {
    
    global_list[[length(global_list) + 1]] <- run_nonparametric_art(df_i)
    pairwise_list[[length(pairwise_list) + 1]] <- run_pairwise_wilcox(df_i)
  }
}

global_tests_full <- bind_rows(global_list) %>%
  mutate(
    p_value_formatted = p_format(p_value),
    significance = p_stars(p_value)
  ) %>%
  arrange(Prefix, Region, Metric_label, Global_test, Term)

pairwise_full <- bind_rows(pairwise_list) %>%
  group_by(Marker, Prefix, Region, Metric, Metric_label) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    p_value_formatted = p_format(p_value),
    p_adj_BH_formatted = p_format(p_adj_BH),
    significance = p_stars(p_value),
    significance_BH = p_stars(p_adj_BH),
    Direction_mean = case_when(
      mean_wildfire > mean_control ~ "Wildfire > Control",
      mean_wildfire < mean_control ~ "Wildfire < Control",
      TRUE ~ "Wildfire = Control"
    ),
    Direction_median = case_when(
      median_wildfire > median_control ~ "Wildfire > Control",
      median_wildfire < median_control ~ "Wildfire < Control",
      TRUE ~ "Wildfire = Control"
    )
  ) %>%
  arrange(Prefix, Region, as.numeric(Year), Metric_label)


global_tests_paper <- global_tests_full %>%
  transmute(
    Marker,
    Region,
    Metric = Metric_label,
    `Test applied` = Global_test,
    Term,
    df1 = round(df1, 3),
    df2 = round(df2, 3),
    Statistic = round(statistic, 4),
    `Partial eta squared` = round(partial_eta2, 4),
    `p value` = p_value_formatted,
    Significance = significance,
    Status = status
  )

pairwise_paper <- pairwise_full %>%
  transmute(
    Marker,
    Region,
    Year,
    Metric = Metric_label,
    `Pairwise test` = Pairwise_test,
    Comparison = comparison,
    `n Control` = n_control,
    `n Wildfire` = n_wildfire,
    `Mean Control` = round2(mean_control),
    `Mean Wildfire` = round2(mean_wildfire),
    `Median Control` = round2(median_control),
    `Median Wildfire` = round2(median_wildfire),
    `Wildfire change (%)` = round(percent_change_wildfire_vs_control, 2),
    Direction = Direction_mean,
    Statistic = round(statistic, 4),
    `p value` = p_value_formatted,
    `p adjusted BH` = p_adj_BH_formatted,
    Significance = significance_BH,
    Status = status
  )

pairwise_significant_paper <- pairwise_paper %>%
  filter(Significance %in% c("*", "**", "***"))

if (nrow(pairwise_significant_paper) == 0) {
  pairwise_significant_paper <- data.frame(
    Message = "No significant pairwise Control vs Wildfire comparisons after BH correction."
  )
}

assumptions_decision_paper <- decision_table %>%
  transmute(
    Marker,
    Region,
    Metric = Metric_label,
    `Residual normality p` = p_format(Residual_normality_p),
    `Residual normality` = Residual_normality_result,
    `Levene p` = p_format(Levene_p),
    `Variance homogeneity` = Variance_result,
    `Minimum n per group` = min_n_evaluable_group,
    `Maximum n per group` = max_n_evaluable_group,
    `Missing genus-level samples` = total_missing_genus_level,
    `Balanced after genus filtering` = balanced_design_after_genus_filter,
    `Selected global test` = case_when(
      Final_test_family == "Parametric_ANOVA" ~ "Type III factorial ANOVA",
      Final_test_family == "Welch_or_robust" ~ "Welch robust alternative",
      Final_test_family == "Nonparametric_ART" ~ "Aligned rank transform ANOVA",
      TRUE ~ "Manual review"
    ),
    `Selected pairwise test` = case_when(
      Final_test_family == "Parametric_ANOVA" ~ "EMMeans contrasts by year",
      Final_test_family == "Welch_or_robust" ~ "Welch t-test by year",
      Final_test_family == "Nonparametric_ART" ~ "Wilcoxon rank-sum by year",
      TRUE ~ "Manual review"
    )
  ) %>%
  arrange(Marker, Region, Metric)


annotation_df <- pairwise_full %>%
  mutate(
    Year_num = as.numeric(Year),
    label = case_when(
      !is.na(p_adj_BH) & p_adj_BH < 0.001 ~ "***",
      !is.na(p_adj_BH) & p_adj_BH < 0.01 ~ "**",
      !is.na(p_adj_BH) & p_adj_BH < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  filter(label != "") %>%
  select(Prefix, Region, Metric_label, Year, Year_num, label, p_adj_BH)

y_positions <- alpha_long %>%
  group_by(Prefix, Region, Metric_label, Year) %>%
  summarise(
    y_max = max(Value, na.rm = TRUE),
    y_min = min(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Prefix, Region, Metric_label) %>%
  mutate(
    range_panel = max(y_max, na.rm = TRUE) - min(y_min, na.rm = TRUE),
    range_panel = ifelse(range_panel == 0, max(y_max, na.rm = TRUE) * 0.15, range_panel),
    y_position = y_max + range_panel * 0.12
  ) %>%
  ungroup()

annotation_df <- annotation_df %>%
  left_join(
    y_positions,
    by = c("Prefix", "Region", "Metric_label", "Year")
  )


make_single_panel <- function(marker_prefix,
                              region_i,
                              metric_i,
                              show_stats = FALSE,
                              show_title = TRUE,
                              show_x_title = TRUE,
                              show_y_title = TRUE) {
  
  df <- alpha_long %>%
    filter(
      as.character(Prefix) == marker_prefix,
      as.character(Region) == region_i,
      as.character(Metric_label) == metric_i
    ) %>%
    mutate(
      Year_factor = factor(
        Year,
        levels = as.character(sort(unique(as.numeric(Year))))
      )
    )
  
  ann <- annotation_df %>%
    filter(
      as.character(Prefix) == marker_prefix,
      as.character(Region) == region_i,
      as.character(Metric_label) == metric_i
    ) %>%
    mutate(
      Year_factor = factor(
        Year,
        levels = as.character(sort(unique(as.numeric(df$Year))))
      )
    )
  
  title_txt <- ifelse(show_title, metric_titles[[metric_i]], "")
  x_txt <- ifelse(show_x_title, "Years since fire", "")
  y_txt <- ifelse(show_y_title, "Effective diversity of genera (Hill number)", "")
  
  p <- ggplot(
    df,
    aes(
      x = Year_factor,
      y = Value,
      fill = Treatment
    )
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.55,
      outlier.shape = NA,
      alpha = 0.80,
      linewidth = 0.45,
      color = "gray25"
    ) +
    geom_jitter(
      aes(color = Treatment),
      position = position_jitterdodge(
        jitter.width = 0.08,
        dodge.width = 0.75
      ),
      size = 1.8,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      aes(group = Treatment),
      position = position_dodge(width = 0.75),
      shape = 18,
      size = 2.8,
      color = "black",
      show.legend = FALSE
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = region_i,
      hjust = -0.10,
      vjust = 1.45,
      fontface = "bold",
      size = 4.2
    ) +
    scale_fill_manual(values = treatment_colors, drop = FALSE) +
    scale_color_manual(values = treatment_colors, drop = FALSE) +
    labs(
      title = title_txt,
      x = x_txt,
      y = y_txt,
      fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10, color = "black"),
      legend.position = "top",
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.55),
      plot.margin = margin(5, 6, 5, 6)
    )
  
  if (show_stats && nrow(ann) > 0) {
    p <- p +
      geom_text(
        data = ann,
        aes(
          x = Year_factor,
          y = y_position,
          label = label
        ),
        inherit.aes = FALSE,
        size = 5,
        fontface = "bold"
      )
  }
  
  return(p)
}

make_marker_figure <- function(marker_prefix, show_stats = FALSE) {
  
  fig_type <- ifelse(show_stats, "WITH_SIGNIFICANCE", "CLEAN")
  
  p_A_D0 <- make_single_panel(marker_prefix, "Amazonas", "Genus richness D0 (q=0)", show_stats, TRUE, FALSE, TRUE)
  p_A_D1 <- make_single_panel(marker_prefix, "Amazonas", "Genus Shannon diversity D1 (q=1)", show_stats, TRUE, FALSE, TRUE)
  p_A_D2 <- make_single_panel(marker_prefix, "Amazonas", "Genus Simpson diversity D2 (q=2)", show_stats, TRUE, FALSE, TRUE)
  
  p_C_D0 <- make_single_panel(marker_prefix, "Cusco", "Genus richness D0 (q=0)", show_stats, FALSE, TRUE, TRUE)
  p_C_D1 <- make_single_panel(marker_prefix, "Cusco", "Genus Shannon diversity D1 (q=1)", show_stats, FALSE, TRUE, TRUE)
  p_C_D2 <- make_single_panel(marker_prefix, "Cusco", "Genus Simpson diversity D2 (q=2)", show_stats, FALSE, TRUE, TRUE)
  
  combined <- (p_A_D0 + p_A_D1 + p_A_D2) /
    (p_C_D0 + p_C_D1 + p_C_D2) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "top",
      legend.justification = "right"
    )
  
  combined <- combined +
    plot_annotation(
      title = paste0(marker_names[[marker_prefix]], " alpha diversity at genus level"),
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      )
    )
  
  png_file <- file.path(fig_dir, paste0(marker_prefix, "_alpha_diversity_GENUS_Hill_D0_D1_D2_", fig_type, ".png"))
  tiff_file <- file.path(fig_dir, paste0(marker_prefix, "_alpha_diversity_GENUS_Hill_D0_D1_D2_", fig_type, ".tiff"))
  pdf_file <- file.path(fig_dir, paste0(marker_prefix, "_alpha_diversity_GENUS_Hill_D0_D1_D2_", fig_type, ".pdf"))
  svg_file <- file.path(fig_dir, paste0(marker_prefix, "_alpha_diversity_GENUS_Hill_D0_D1_D2_", fig_type, ".svg"))
  
  ggsave(png_file, plot = combined, width = 17, height = 7.2, dpi = 400)
  ggsave(tiff_file, plot = combined, width = 17, height = 7.2, dpi = 400, compression = "lzw")
  ggsave(pdf_file, plot = combined, width = 17, height = 7.2)
  ggsave(svg_file, plot = combined, width = 17, height = 7.2)
  
  cat("\nFigura generada:", marker_prefix, fig_type, "\n")
  cat(png_file, "\n")
  
  return(combined)
}

p_bac_clean <- make_marker_figure("BAC", show_stats = FALSE)
p_arq_clean <- make_marker_figure("ARQ", show_stats = FALSE)
p_its_clean <- make_marker_figure("ITS", show_stats = FALSE)

p_bac_sig <- make_marker_figure("BAC", show_stats = TRUE)
p_arq_sig <- make_marker_figure("ARQ", show_stats = TRUE)
p_its_sig <- make_marker_figure("ITS", show_stats = TRUE)


write.csv(descriptive_full, file.path(tab_dir, "01_descriptive_statistics_full_GENUS_Hill.csv"), row.names = FALSE)
write.csv(descriptive_paper, file.path(tab_dir, "02_descriptive_statistics_paper_GENUS_Hill.csv"), row.names = FALSE)
write.csv(assumptions_decision_paper, file.path(tab_dir, "03_assumptions_and_selected_tests_GENUS.csv"), row.names = FALSE)
write.csv(global_tests_full, file.path(tab_dir, "04_global_tests_full_GENUS_ANOVA_or_ART.csv"), row.names = FALSE)
write.csv(global_tests_paper, file.path(tab_dir, "05_global_tests_paper_GENUS.csv"), row.names = FALSE)
write.csv(pairwise_full, file.path(tab_dir, "06_pairwise_Control_vs_Wildfire_by_year_full_GENUS.csv"), row.names = FALSE)
write.csv(pairwise_paper, file.path(tab_dir, "07_pairwise_Control_vs_Wildfire_by_year_paper_GENUS.csv"), row.names = FALSE)
write.csv(pairwise_significant_paper, file.path(tab_dir, "08_pairwise_significant_only_GENUS.csv"), row.names = FALSE)
write.csv(alpha_long, file.path(tab_dir, "09_alpha_values_long_GENUS_Hill.csv"), row.names = FALSE)

if (file.exists(genus_qc_marker_file)) {
  genus_qc_marker <- read.csv(genus_qc_marker_file, stringsAsFactors = FALSE)
  write.csv(genus_qc_marker, file.path(supp_dir, "S1_genus_level_QC_by_marker.csv"), row.names = FALSE)
} else {
  genus_qc_marker <- data.frame(Message = "QC marker file not found")
}

if (file.exists(genus_qc_sample_file)) {
  genus_qc_sample <- read.csv(genus_qc_sample_file, stringsAsFactors = FALSE)
  write.csv(genus_qc_sample, file.path(supp_dir, "S2_genus_level_QC_by_sample.csv"), row.names = FALSE)
} else {
  genus_qc_sample <- data.frame(Message = "QC sample file not found")
}


xlsx_main <- file.path(
  excel_dir,
  "MAIN_TABLES_ALPHA_DIVERSITY_GENUS_FINAL_ANALYSIS.xlsx"
)

wb <- createWorkbook()

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom",
  halign = "center",
  valign = "center",
  wrapText = TRUE
)

body_style <- createStyle(
  valign = "center",
  wrapText = TRUE
)

add_sheet <- function(wb, sheet_name, data) {
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data, headerStyle = header_style)
  
  addStyle(
    wb,
    sheet = sheet_name,
    style = body_style,
    rows = 1:(nrow(data) + 1),
    cols = 1:ncol(data),
    gridExpand = TRUE
  )
  
  freezePane(wb, sheet_name, firstRow = TRUE)
  addFilter(wb, sheet_name, rows = 1, cols = 1:ncol(data))
  setColWidths(wb, sheet_name, cols = 1:ncol(data), widths = "auto")
}

readme_main <- data.frame(
  Section = c(
    "Level",
    "Alpha diversity metrics",
    "Main figure",
    "Parametric analyses",
    "Non-parametric analyses",
    "Pairwise comparisons",
    "Multiple testing correction",
    "Recommended main material",
    "Recommended supplementary material"
  ),
  Description = c(
    "All alpha diversity analyses were performed at the genus level using only identified genera.",
    "D0 q=0 = genus richness; D1 q=1 = exp(Shannon) at genus level; D2 q=2 = inverse Simpson at genus level.",
    "Use CLEAN figures as the main paper figures.",
    "When residual normality and variance homogeneity were met, Type III factorial ANOVA was applied.",
    "When residual normality was not met, aligned rank transform ANOVA was applied.",
    "Control vs Wildfire comparisons were performed within each post-fire year using emmeans for ANOVA cases and Wilcoxon tests for non-parametric cases.",
    "Pairwise p-values were adjusted using the Benjamini-Hochberg method within each marker, region and metric.",
    "Use global tests table and significant pairwise comparisons if needed.",
    "Use assumptions, QC, descriptive statistics and complete pairwise results as supplementary tables."
  )
)

add_sheet(wb, "00_Read_me", readme_main)
add_sheet(wb, "01_Global_tests_paper", global_tests_paper)
add_sheet(wb, "02_Pairwise_sig_only", pairwise_significant_paper)
add_sheet(wb, "03_Pairwise_all_paper", pairwise_paper)
add_sheet(wb, "04_Descriptive_paper", descriptive_paper)
add_sheet(wb, "05_Assumptions_tests", assumptions_decision_paper)
add_sheet(wb, "06_Global_tests_full", global_tests_full)
add_sheet(wb, "07_Pairwise_full", pairwise_full)

saveWorkbook(wb, xlsx_main, overwrite = TRUE)


xlsx_supp <- file.path(
  supp_dir,
  "SUPPLEMENTARY_ALPHA_DIVERSITY_GENUS_TABLES.xlsx"
)

wb2 <- createWorkbook()

readme_supp <- data.frame(
  Supplementary_table = c(
    "Table S1",
    "Table S2",
    "Table S3",
    "Table S4",
    "Table S5",
    "Table S6",
    "Table S7"
  ),
  Content = c(
    "Genus-level QC by marker.",
    "Genus-level QC by sample.",
    "Statistical assumptions and selected tests.",
    "Descriptive statistics of genus-level Hill numbers.",
    "Global tests for Treatment, Year and Treatment × Year.",
    "Complete pairwise Control vs Wildfire comparisons by year.",
    "Raw alpha diversity values at genus level."
  ),
  Suggested_use = c(
    "Important for archaeal interpretation.",
    "Important to document samples without identified genus.",
    "Justifies use of mixed statistical tests.",
    "Supplementary descriptive support.",
    "Can be supplementary or main depending on journal space.",
    "Best placed in supplementary because it is long.",
    "For transparency and reproducibility."
  )
)

add_sheet(wb2, "00_Read_me", readme_supp)
add_sheet(wb2, "S1_Genus_QC_marker", genus_qc_marker)
add_sheet(wb2, "S2_Genus_QC_sample", genus_qc_sample)
add_sheet(wb2, "S3_Assumptions_tests", assumptions_decision_paper)
add_sheet(wb2, "S4_Descriptive", descriptive_paper)
add_sheet(wb2, "S5_Global_tests", global_tests_paper)
add_sheet(wb2, "S6_Pairwise_all", pairwise_paper)
add_sheet(wb2, "S7_Raw_alpha_values", alpha_long)

if (file.exists(alpha_wide_file)) {
  alpha_wide <- read.csv(alpha_wide_file, stringsAsFactors = FALSE)
  add_sheet(wb2, "S8_Raw_alpha_wide", alpha_wide)
}

saveWorkbook(wb2, xlsx_supp, overwrite = TRUE)


methods_text <- c(
  "Alpha diversity was assessed at the genus level using only identified genera. ASV counts were first aggregated by genus, excluding unclassified, uncultured, unidentified, unknown, ambiguous, chloroplast and mitochondrial assignments. Diversity was expressed as Hill numbers, including genus richness D0 (q = 0), Shannon diversity D1 (q = 1; exp[Shannon entropy]) and Simpson diversity D2 (q = 2; inverse Simpson index).",
  "",
  "For each molecular marker, region and diversity metric, statistical tests were selected according to model assumptions. Normality was evaluated using Shapiro-Wilk tests on model residuals, and variance homogeneity was evaluated using Levene's test. When residual normality and variance homogeneity were satisfied, Type III factorial ANOVA was fitted with treatment, post-fire year and their interaction as fixed factors. When residual normality was not satisfied, aligned rank transform ANOVA was used to retain the factorial structure. Pairwise comparisons between control and wildfire-affected soils within each post-fire year were performed using estimated marginal means for parametric models and Wilcoxon rank-sum tests for non-parametric cases. Pairwise p-values were adjusted using the Benjamini-Hochberg method."
)

caption_text <- c(
  "Figure X. Genus-level alpha diversity of soil microbial communities expressed as Hill numbers. D0 (q = 0) represents genus richness, D1 (q = 1) represents Shannon effective diversity and D2 (q = 2) represents Simpson effective diversity. Boxes represent the interquartile range, horizontal lines indicate medians, black diamonds indicate means and points represent individual samples."
)

writeLines(methods_text, file.path(text_dir, "METHODS_alpha_diversity_GENUS.txt"))
writeLines(caption_text, file.path(text_dir, "FIGURE_CAPTION_alpha_diversity_GENUS.txt"))


stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

zip_file <- file.path(
  project_root,
  paste0("ALPHA_DIVERSITY_GENUS_FINAL_FOR_PAPER_", stamp, ".zip")
)

old_wd <- getwd()
setwd(project_root)

zip::zipr(
  zipfile = zip_file,
  files = basename(outdir)
)

setwd(old_wd)


cat("\n=====================================================\n")
cat("Genus-level alpha diversity analyses generated\n")
cat("=====================================================\n")

cat("\nCarpeta principal:\n")
cat(outdir, "\n")

cat("\nFiguras:\n")
cat(fig_dir, "\n")

cat("\nExcel principal:\n")
cat(xlsx_main, "\n")

cat("\nExcel suplementario:\n")
cat(xlsx_supp, "\n")

cat("\nTextos para manuscrito:\n")
cat(text_dir, "\n")

cat("\nZIP final:\n")
cat(zip_file, "\n")

cat("\nDescarga en RStudio:\n")
cat(basename(zip_file), "\n")

cat("\nResumen de pruebas globales aplicadas:\n")
print(
  global_tests_full %>%
    count(Marker, Prefix, Region, Metric_label, Global_test)
)

cat("\nPairwise comparison summary by years since fire:\n")
print(
  pairwise_full %>%
    count(Marker, Prefix, Region, Metric_label, Pairwise_test)
)

cat("\nFiguras principales limpias:\n")
cat("BAC_alpha_diversity_GENUS_Hill_D0_D1_D2_CLEAN.png/tiff/pdf/svg\n")
cat("ARQ_alpha_diversity_GENUS_Hill_D0_D1_D2_CLEAN.png/tiff/pdf/svg\n")
cat("ITS_alpha_diversity_GENUS_Hill_D0_D1_D2_CLEAN.png/tiff/pdf/svg\n")

cat("\nFiguras con significancia:\n")
cat("BAC_alpha_diversity_GENUS_Hill_D0_D1_D2_WITH_SIGNIFICANCE.png/tiff/pdf/svg\n")
cat("ARQ_alpha_diversity_GENUS_Hill_D0_D1_D2_WITH_SIGNIFICANCE.png/tiff/pdf/svg\n")
cat("ITS_alpha_diversity_GENUS_Hill_D0_D1_D2_WITH_SIGNIFICANCE.png/tiff/pdf/svg\n")

cat("\nStatistical interpretation:\n")
cat("Usa las figuras CLEAN como principales.\n")
cat("Usa el Excel principal para resultados.\n")
cat("Usa el Excel suplementario para supuestos, QC y tablas completas.\n")

cat("=====================================================\n")

