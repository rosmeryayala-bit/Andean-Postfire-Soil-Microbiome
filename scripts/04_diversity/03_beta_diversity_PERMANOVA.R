# Beta diversity
# Bray-Curtis dissimilarity, PCoA, PERMANOVA and dispersion tests.
# Input: project-specific phyloseq objects, site metadata, or soil data tables.
# Output: analysis tables and figures as defined by the project directory settings.


rm(list = ls())
options(stringsAsFactors = FALSE)


suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(vegan)
  library(ape)
  library(zip)
  library(svglite)
  library(ggrepel)
})

cat("\n=====================================================\n")
cat("BETA DIVERSITY AT GENUS LEVEL\n")
cat("Bray-Curtis | PCoA | PERMANOVA | betadisper\n")
cat("=====================================================\n")


project_root <- path.expand("~/Rosmery/Tesis Rosmery-Metabarcoding")

PERMUTATIONS <- 999
MIN_GENUS_TOTAL_READS <- 1
MIN_GENUS_PREVALENCE <- 1

outdir <- file.path(project_root, "RESULTADOS_BETA_DIVERSITY_GENUS_FINAL_PAPER")

fig_dir   <- file.path(outdir, "01_FIGURES")
tab_dir   <- file.path(outdir, "02_TABLES_CSV")
excel_dir <- file.path(outdir, "03_EXCEL_TABLES")
supp_dir  <- file.path(outdir, "04_SUPPLEMENTARY")
text_dir  <- file.path(outdir, "05_TEXT_FOR_MANUSCRIPT")
data_dir  <- file.path(outdir, "06_GENUS_MATRICES")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(excel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)


find_existing_path <- function(paths) {
  paths <- path.expand(paths)
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

fix_duplicate_names <- function(df) {
  rn <- rownames(df)
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- make.unique(colnames(df))
  rownames(df) <- rn
  df
}

find_col <- function(names_vec, candidates) {
  hit <- candidates[candidates %in% names_vec][1]
  if (length(hit) == 0 || is.na(hit)) {
    hit2 <- names_vec[tolower(names_vec) %in% tolower(candidates)][1]
    if (length(hit2) == 0 || is.na(hit2)) return(NA_character_)
    return(hit2)
  }
  hit
}

normalize_region <- function(x) {
  x <- trimws(as.character(x))
  ifelse(
    grepl("amazon", x, ignore.case = TRUE),
    "Amazonas",
    ifelse(grepl("cusco|cuzco", x, ignore.case = TRUE), "Cusco", x)
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
  x <- sub("_(BAC|ARQ|FUN|ITS)$", "", x)
  x <- sub("_(BAC|ARQ|FUN|ITS)\\.raw$", "", x)
  x <- gsub("_", "", x)
  x
}

is_bad_genus <- function(x) {
  x0 <- trimws(as.character(x))
  xl <- tolower(x0)
  
  is.na(x0) |
    x0 == "" |
    xl %in% c(
      "na", "n/a", "none", "null",
      "unclassified", "unassigned", "uncultured",
      "unknown", "unidentified", "metagenome",
      "norank", "no rank", "ambiguous"
    ) |
    grepl(
      "uncultured|unclassified|unassigned|unknown|unidentified|metagenome|norank|no rank|ambiguous|incertae|chloroplast|mitochondria",
      xl,
      ignore.case = TRUE
    )
}

remove_contaminant_taxa <- function(ps) {
  
  tax_df <- tryCatch(
    as.data.frame(tax_table(ps), stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(tax_df)) return(ps)
  
  tax_df <- fix_duplicate_names(tax_df)
  
  bad_taxa <- apply(
    tax_df,
    1,
    function(z) {
      any(grepl("chloroplast|mitochondria", paste(z, collapse = " "), ignore.case = TRUE))
    }
  )
  
  bad_taxa[is.na(bad_taxa)] <- FALSE
  
  if (any(bad_taxa)) {
    cat("\nTaxa removed as chloroplast/mitochondria:", sum(bad_taxa), "\n")
    ps <- prune_taxa(!bad_taxa, ps)
  }
  
  ps
}

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

get_first_col <- function(df, patterns) {
  for (pat in patterns) {
    hit <- grep(pat, colnames(df), value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

safe_num <- function(df, row_i, col_i) {
  if (length(row_i) == 0 || is.na(col_i) || !(col_i %in% colnames(df))) return(NA_real_)
  suppressWarnings(as.numeric(df[[col_i]][row_i[1]]))
}

safe_term <- function(df, row_i) {
  if (length(row_i) == 0) return(NA_character_)
  as.character(df$Term[row_i[1]])
}

extract_adonis_row <- function(adonis_obj, term_pattern = "Treatment") {
  
  ad_tab <- as.data.frame(adonis_obj)
  ad_tab$Term <- rownames(ad_tab)
  rownames(ad_tab) <- NULL
  
  df_col <- get_first_col(ad_tab, c("^Df$"))
  ss_col <- get_first_col(ad_tab, c("SumOfSqs", "Sum Sq", "SumsOfSqs"))
  r2_col <- get_first_col(ad_tab, c("^R2$"))
  f_col  <- get_first_col(ad_tab, c("^F$", "F.Model", "F value", "F.value"))
  p_col  <- get_first_col(ad_tab, c("Pr\\(>F\\)", "^Pr", "Pr\\."))
  
  row_i <- grep(term_pattern, ad_tab$Term, ignore.case = TRUE)
  
  if (length(row_i) == 0) {
    row_i <- which(!grepl("Residual|Residuals|Total", ad_tab$Term, ignore.case = TRUE))
  }
  
  if (length(row_i) == 0) {
    return(
      data.frame(
        Term_found = NA_character_,
        Df = NA_real_,
        SumOfSqs = NA_real_,
        R2 = NA_real_,
        F_value = NA_real_,
        p_value = NA_real_,
        Status_extract = "No model term found",
        stringsAsFactors = FALSE
      )
    )
  }
  
  data.frame(
    Term_found = safe_term(ad_tab, row_i),
    Df = safe_num(ad_tab, row_i, df_col),
    SumOfSqs = safe_num(ad_tab, row_i, ss_col),
    R2 = safe_num(ad_tab, row_i, r2_col),
    F_value = safe_num(ad_tab, row_i, f_col),
    p_value = safe_num(ad_tab, row_i, p_col),
    Status_extract = "OK",
    stringsAsFactors = FALSE
  )
}

format_adonis_result <- function(adonis_obj,
                                 marker,
                                 prefix,
                                 region,
                                 model_scope,
                                 model_formula,
                                 distance_method = "Bray-Curtis") {
  
  tab <- as.data.frame(adonis_obj)
  tab$Term <- rownames(tab)
  rownames(tab) <- NULL
  
  df_col <- get_first_col(tab, c("^Df$"))
  ss_col <- get_first_col(tab, c("SumOfSqs", "Sum Sq", "SumsOfSqs"))
  r2_col <- get_first_col(tab, c("^R2$"))
  f_col  <- get_first_col(tab, c("^F$", "F.Model", "F value", "F.value"))
  p_col  <- get_first_col(tab, c("Pr\\(>F\\)", "^Pr", "Pr\\."))
  
  keep_rows <- which(!grepl("Residual|Residuals|Total", tab$Term, ignore.case = TRUE))
  
  if (length(keep_rows) == 0) {
    return(
      data.frame(
        Marker = marker,
        Prefix = prefix,
        Region = region,
        Model_scope = model_scope,
        Model_formula = model_formula,
        Distance = distance_method,
        Term = NA_character_,
        Df = NA_real_,
        SumOfSqs = NA_real_,
        R2 = NA_real_,
        F_value = NA_real_,
        p_value = NA_real_,
        p_value_formatted = NA_character_,
        Significance = "",
        Permutations = PERMUTATIONS,
        Status = "No model terms found",
        stringsAsFactors = FALSE
      )
    )
  }
  
  out <- lapply(
    keep_rows,
    function(i) {
      pval <- safe_num(tab, i, p_col)
      data.frame(
        Marker = marker,
        Prefix = prefix,
        Region = region,
        Model_scope = model_scope,
        Model_formula = model_formula,
        Distance = distance_method,
        Term = as.character(tab$Term[i]),
        Df = safe_num(tab, i, df_col),
        SumOfSqs = safe_num(tab, i, ss_col),
        R2 = safe_num(tab, i, r2_col),
        F_value = safe_num(tab, i, f_col),
        p_value = pval,
        p_value_formatted = p_format(pval),
        Significance = p_stars(pval),
        Permutations = PERMUTATIONS,
        Status = "OK",
        stringsAsFactors = FALSE
      )
    }
  )
  
  bind_rows(out)
}


prepare_genus_beta_marker <- function(marker_label,
                                      base_candidates,
                                      ps_candidates,
                                      file_prefix) {
  
  cat("\n=====================================================\n")
  cat("PROCESSING:", marker_label, "\n")
  cat("=====================================================\n")
  
  base <- find_existing_path(base_candidates)
  if (is.na(base)) stop("Folder not found for: ", marker_label)
  
  ps_file <- find_existing_path(file.path(base, ps_candidates))
  if (is.na(ps_file)) stop("Phyloseq object not found for: ", marker_label)
  
  cat("\nFolder used:\n")
  cat(base, "\n")
  cat("\nPhyloseq object used:\n")
  cat(ps_file, "\n")
  
  ps <- readRDS(ps_file)
  
  cat("\nLoaded object:\n")
  cat("Samples:", nsamples(ps), "\n")
  cat("Initial ASVs:", ntaxa(ps), "\n")
  cat("Initial reads:", sum(sample_sums(ps)), "\n")
  
  meta_df <- data.frame(sample_data(ps), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(meta_df) <- sample_names(ps)
  meta_df <- fix_duplicate_names(meta_df)
  
  region_col <- find_col(colnames(meta_df), c("Region", "region", "Departamento", "departamento"))
  treat_col  <- find_col(colnames(meta_df), c("Treatment", "treatment", "Tratamiento", "tratamiento", "Tratament"))
  year_col   <- find_col(colnames(meta_df), c("Years_since_fire", "Years_since_fire.1", "YearsF", "Year", "year", "Año", "Anio"))
  
  if (is.na(region_col)) stop("Region column not found for: ", marker_label)
  if (is.na(treat_col))  stop("Treatment column not found for: ", marker_label)
  if (is.na(year_col))   stop("Post-fire year column not found for: ", marker_label)
  
  meta_df$SampleID <- rownames(meta_df)
  meta_df$SampleLabel <- clean_sample_label(meta_df$SampleID)
  meta_df$Region <- normalize_region(meta_df[[region_col]])
  meta_df$Treatment <- normalize_treatment(meta_df[[treat_col]])
  meta_df$Year <- normalize_year(meta_df[[year_col]])
  meta_df$Year_num <- suppressWarnings(as.numeric(meta_df$Year))
  
  keep_samples <- meta_df$SampleID[
    meta_df$Region %in% c("Amazonas", "Cusco") &
      meta_df$Treatment %in% c("Control", "Wildfire")
  ]
  
  ps <- prune_samples(keep_samples, ps)
  ps <- remove_contaminant_taxa(ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  
  meta_df <- meta_df[sample_names(ps), , drop = FALSE]
  meta_df <- fix_duplicate_names(meta_df)
  
  cat("\nSamples by region, treatment and post-fire year:\n")
  print(table(meta_df$Region, meta_df$Treatment, meta_df$Year))
  
  tax_df <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE, check.names = FALSE)
  tax_df <- fix_duplicate_names(tax_df)
  tax_df$ASV_ID <- taxa_names(ps)
  
  genus_col <- find_col(
    colnames(tax_df),
    c("Genus", "genus", "GENUS", "Género", "Genero", "Rank6", "rank6")
  )
  
  if (is.na(genus_col)) {
    stop("Genus column not found in taxonomy for: ", marker_label)
  }
  
  tax_df$Genus_clean <- trimws(as.character(tax_df[[genus_col]]))
  tax_df$Genus_clean[is_bad_genus(tax_df$Genus_clean)] <- NA
  
  identified_taxa <- tax_df$ASV_ID[!is.na(tax_df$Genus_clean)]
  
  cat("\nASVs with identified genus:", length(identified_taxa), "\n")
  cat("ASVs without identified genus:", ntaxa(ps) - length(identified_taxa), "\n")
  
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  otu_mat <- otu_mat[, rownames(meta_df), drop = FALSE]
  
  total_reads_all <- colSums(otu_mat, na.rm = TRUE)
  
  otu_identified <- otu_mat[rownames(otu_mat) %in% identified_taxa, , drop = FALSE]
  
  tax_identified <- tax_df %>%
    filter(ASV_ID %in% rownames(otu_identified)) %>%
    select(ASV_ID, Genus_clean)
  
  tax_identified <- tax_identified[match(rownames(otu_identified), tax_identified$ASV_ID), ]
  
  if (nrow(otu_identified) > 0) {
    genus_mat <- rowsum(
      otu_identified,
      group = tax_identified$Genus_clean,
      reorder = TRUE
    )
  } else {
    genus_mat <- matrix(
      0,
      nrow = 0,
      ncol = ncol(otu_mat),
      dimnames = list(character(0), colnames(otu_mat))
    )
  }
  
  genus_mat <- as.matrix(genus_mat)
  
  genus_total_reads <- rowSums(genus_mat, na.rm = TRUE)
  genus_prevalence <- rowSums(genus_mat > 0, na.rm = TRUE)
  
  keep_genera <- genus_total_reads >= MIN_GENUS_TOTAL_READS &
    genus_prevalence >= MIN_GENUS_PREVALENCE
  
  genus_mat <- genus_mat[keep_genera, , drop = FALSE]
  
  total_reads_identified_genus <- colSums(genus_mat, na.rm = TRUE)
  
  percent_reads_identified_genus <- ifelse(
    total_reads_all > 0,
    total_reads_identified_genus / total_reads_all * 100,
    NA_real_
  )
  
  samples_without_identified_genus <- names(total_reads_identified_genus)[total_reads_identified_genus == 0]
  
  if (length(samples_without_identified_genus) > 0) {
    cat("\nWARNING: samples without reads assigned to identified genera:\n")
    print(samples_without_identified_genus)
  }
  
  genus_counts <- t(genus_mat)
  genus_counts <- as.matrix(genus_counts)
  
  evaluable_samples <- rownames(genus_counts)[rowSums(genus_counts, na.rm = TRUE) > 0]
  
  genus_counts_eval <- genus_counts[evaluable_samples, , drop = FALSE]
  meta_eval <- meta_df[evaluable_samples, , drop = FALSE]
  
  meta_eval$SampleID <- rownames(meta_eval)
  
  genus_rel <- sweep(genus_counts_eval, 1, rowSums(genus_counts_eval, na.rm = TRUE), "/")
  genus_rel[!is.finite(genus_rel)] <- 0
  
  genus_counts_df <- as.data.frame(genus_counts_eval, check.names = FALSE)
  genus_counts_df$SampleID <- rownames(genus_counts_df)
  genus_counts_df <- genus_counts_df %>% select(SampleID, everything())
  
  genus_rel_df <- as.data.frame(genus_rel, check.names = FALSE)
  genus_rel_df$SampleID <- rownames(genus_rel_df)
  genus_rel_df <- genus_rel_df %>% select(SampleID, everything())
  
  genus_qc_marker <- data.frame(
    Marker = marker_label,
    Prefix = file_prefix,
    Total_samples_initial = nsamples(ps),
    Samples_used_beta_diversity = nrow(meta_eval),
    Samples_excluded_no_identified_genus = length(samples_without_identified_genus),
    Total_ASVs_after_filtering = ntaxa(ps),
    ASVs_with_identified_genus = length(identified_taxa),
    Identified_genera_used = ncol(genus_rel),
    Mean_percent_reads_with_identified_genus = mean(percent_reads_identified_genus, na.rm = TRUE),
    Min_percent_reads_with_identified_genus = min(percent_reads_identified_genus, na.rm = TRUE),
    Max_percent_reads_with_identified_genus = max(percent_reads_identified_genus, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  sample_qc <- data.frame(
    Marker = marker_label,
    Prefix = file_prefix,
    SampleID = colnames(otu_mat),
    Reads_total_all_ASV = as.numeric(total_reads_all),
    Reads_with_identified_genus = as.numeric(total_reads_identified_genus),
    Percent_reads_with_identified_genus = as.numeric(percent_reads_identified_genus),
    Beta_diversity_status = ifelse(
      total_reads_identified_genus > 0,
      "Included",
      "Excluded: no reads assigned to identified genera"
    ),
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      meta_df %>%
        select(SampleID, SampleLabel, Region, Treatment, Year, Year_num),
      by = "SampleID"
    )
  
  list(
    marker = marker_label,
    prefix = file_prefix,
    meta = meta_eval,
    genus_counts = genus_counts_eval,
    genus_rel = genus_rel,
    genus_counts_df = genus_counts_df,
    genus_rel_df = genus_rel_df,
    genus_qc_marker = genus_qc_marker,
    sample_qc = sample_qc
  )
}


bac <- prepare_genus_beta_marker(
  marker_label = "Bacteria 16S",
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

arq <- prepare_genus_beta_marker(
  marker_label = "Archaea 16S",
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

its <- prepare_genus_beta_marker(
  marker_label = "Fungi ITS",
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

write.csv(bac$genus_counts_df, file.path(data_dir, "BAC_genus_counts_sample_by_genus.csv"), row.names = FALSE)
write.csv(arq$genus_counts_df, file.path(data_dir, "ARQ_genus_counts_sample_by_genus.csv"), row.names = FALSE)
write.csv(its$genus_counts_df, file.path(data_dir, "ITS_genus_counts_sample_by_genus.csv"), row.names = FALSE)

write.csv(bac$genus_rel_df, file.path(data_dir, "BAC_genus_relative_abundance_sample_by_genus.csv"), row.names = FALSE)
write.csv(arq$genus_rel_df, file.path(data_dir, "ARQ_genus_relative_abundance_sample_by_genus.csv"), row.names = FALSE)
write.csv(its$genus_rel_df, file.path(data_dir, "ITS_genus_relative_abundance_sample_by_genus.csv"), row.names = FALSE)

genus_qc_marker_all <- bind_rows(bac$genus_qc_marker, arq$genus_qc_marker, its$genus_qc_marker)
sample_qc_all <- bind_rows(bac$sample_qc, arq$sample_qc, its$sample_qc)


run_dispersion_analysis <- function(dist_mat,
                                    group_factor,
                                    meta,
                                    marker,
                                    prefix,
                                    region,
                                    model_scope,
                                    group_definition) {
  
  group_factor <- factor(group_factor)
  
  if (length(unique(group_factor)) < 2 || length(group_factor) != nrow(meta)) {
    return(
      list(
        test = data.frame(
          Marker = marker,
          Prefix = prefix,
          Region = region,
          Model_scope = model_scope,
          Group_definition = group_definition,
          Dispersion_test = "betadisper + permutest",
          F_value = NA_real_,
          p_value = NA_real_,
          p_value_formatted = NA_character_,
          Dispersion_result = "Not evaluable",
          Permutations = PERMUTATIONS,
          stringsAsFactors = FALSE
        ),
        distances = data.frame()
      )
    )
  }
  
  bd <- tryCatch(
    vegan::betadisper(dist_mat, group_factor),
    error = function(e) NULL
  )
  
  if (is.null(bd)) {
    return(
      list(
        test = data.frame(
          Marker = marker,
          Prefix = prefix,
          Region = region,
          Model_scope = model_scope,
          Group_definition = group_definition,
          Dispersion_test = "betadisper + permutest",
          F_value = NA_real_,
          p_value = NA_real_,
          p_value_formatted = NA_character_,
          Dispersion_result = "Error in betadisper",
          Permutations = PERMUTATIONS,
          stringsAsFactors = FALSE
        ),
        distances = data.frame()
      )
    )
  }
  
  pt <- tryCatch(
    vegan::permutest(bd, permutations = PERMUTATIONS),
    error = function(e) NULL
  )
  
  if (is.null(pt)) {
    f_value <- NA_real_
    p_value <- NA_real_
    dispersion_result <- "Error in permutest"
  } else {
    pt_tab <- as.data.frame(pt$tab)
    f_col <- get_first_col(pt_tab, c("^F$", "F value", "F.value"))
    p_col <- get_first_col(pt_tab, c("Pr\\(>F\\)", "^Pr", "Pr\\."))
    f_value <- safe_num(pt_tab, 1, f_col)
    p_value <- safe_num(pt_tab, 1, p_col)
    dispersion_result <- ifelse(
      is.na(p_value),
      "Not evaluable",
      ifelse(p_value >= 0.05, "Homogeneous dispersion", "Heterogeneous dispersion")
    )
  }
  
  test_out <- data.frame(
    Marker = marker,
    Prefix = prefix,
    Region = region,
    Model_scope = model_scope,
    Group_definition = group_definition,
    Dispersion_test = "betadisper + permutest",
    F_value = f_value,
    p_value = p_value,
    p_value_formatted = p_format(p_value),
    Dispersion_result = dispersion_result,
    Permutations = PERMUTATIONS,
    stringsAsFactors = FALSE
  )
  
  dist_out <- data.frame(
    Marker = marker,
    Prefix = prefix,
    Region = region,
    Model_scope = model_scope,
    SampleID = names(bd$distances),
    Distance_to_centroid = as.numeric(bd$distances),
    Group = as.character(group_factor),
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      meta %>% select(SampleID, SampleLabel, Region, Treatment, Year, Year_num),
      by = "SampleID"
    )
  
  list(test = test_out, distances = dist_out)
}

run_pcoa <- function(dist_mat, meta, marker, prefix) {
  
  pco <- tryCatch(
    ape::pcoa(dist_mat, correction = "cailliez"),
    error = function(e) ape::pcoa(dist_mat)
  )
  
  rel_eig <- pco$values$Relative_eig
  
  if (length(rel_eig) < 2 || all(is.na(rel_eig))) {
    eig <- pco$values$Eigenvalues
    eig_pos <- eig[eig > 0]
    rel_eig <- eig / sum(eig_pos)
  }
  
  coords <- data.frame(
    SampleID = rownames(pco$vectors),
    PCoA1 = pco$vectors[, 1],
    PCoA2 = pco$vectors[, 2],
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      meta %>% select(SampleID, SampleLabel, Region, Treatment, Year, Year_num),
      by = "SampleID"
    ) %>%
    mutate(
      Marker = marker,
      Prefix = prefix,
      PCoA1_percent = round(rel_eig[1] * 100, 2),
      PCoA2_percent = round(rel_eig[2] * 100, 2)
    )
  
  variance <- data.frame(
    Marker = marker,
    Prefix = prefix,
    PCoA1_percent = round(rel_eig[1] * 100, 2),
    PCoA2_percent = round(rel_eig[2] * 100, 2),
    stringsAsFactors = FALSE
  )
  
  list(coords = coords, variance = variance)
}


run_beta_marker <- function(marker_obj) {
  
  marker <- marker_obj$marker
  prefix <- marker_obj$prefix
  meta <- marker_obj$meta
  genus_rel <- marker_obj$genus_rel
  
  meta <- meta[rownames(genus_rel), , drop = FALSE]
  meta$SampleID <- rownames(meta)
  
  dist_all <- vegan::vegdist(genus_rel, method = "bray")
  
  pcoa_all <- run_pcoa(dist_all, meta, marker, prefix)
  
  permanova_list <- list()
  dispersion_test_list <- list()
  dispersion_distance_list <- list()
  pairwise_list <- list()
  

  if (length(unique(meta$Region)) >= 2 && length(unique(meta$Treatment)) >= 2) {
    
    meta_overall <- meta %>%
      mutate(
        Region = factor(Region, levels = c("Amazonas", "Cusco")),
        Treatment = factor(Treatment, levels = c("Control", "Wildfire"))
      )
    rownames(meta_overall) <- meta_overall$SampleID
    
    ad_overall <- tryCatch(
      vegan::adonis2(
        dist_all ~ Region * Treatment,
        data = meta_overall,
        permutations = PERMUTATIONS,
        by = "terms"
      ),
      error = function(e) NULL
    )
    
    if (!is.null(ad_overall)) {
      permanova_list[[length(permanova_list) + 1]] <- format_adonis_result(
        ad_overall,
        marker = marker,
        prefix = prefix,
        region = "All regions",
        model_scope = "Overall marker-level model",
        model_formula = "Bray-Curtis ~ Region * Treatment"
      )
    }
    
    disp_overall <- run_dispersion_analysis(
      dist_mat = dist_all,
      group_factor = interaction(meta_overall$Region, meta_overall$Treatment, drop = TRUE),
      meta = meta_overall,
      marker = marker,
      prefix = prefix,
      region = "All regions",
      model_scope = "Overall marker-level model",
      group_definition = "Region × Treatment"
    )
    
    dispersion_test_list[[length(dispersion_test_list) + 1]] <- disp_overall$test
    dispersion_distance_list[[length(dispersion_distance_list) + 1]] <- disp_overall$distances
  }
  

  for (reg_i in c("Amazonas", "Cusco")) {
    
    idx <- meta$Region == reg_i
    
    if (sum(idx) < 6) next
    
    meta_reg <- meta[idx, , drop = FALSE] %>%
      mutate(
        Treatment = factor(Treatment, levels = c("Control", "Wildfire")),
        Year = factor(Year, levels = as.character(sort(unique(as.numeric(as.character(Year))))))
      )
    
    rownames(meta_reg) <- meta_reg$SampleID
    
    genus_rel_reg <- genus_rel[rownames(meta_reg), , drop = FALSE]
    
    if (length(unique(meta_reg$Treatment)) < 2 || length(unique(meta_reg$Year)) < 2) next
    
    dist_reg <- vegan::vegdist(genus_rel_reg, method = "bray")
    
    ad_reg <- tryCatch(
      vegan::adonis2(
        dist_reg ~ Treatment * Year,
        data = meta_reg,
        permutations = PERMUTATIONS,
        by = "terms"
      ),
      error = function(e) NULL
    )
    
    if (!is.null(ad_reg)) {
      permanova_list[[length(permanova_list) + 1]] <- format_adonis_result(
        ad_reg,
        marker = marker,
        prefix = prefix,
        region = reg_i,
        model_scope = "Region-specific model",
        model_formula = "Bray-Curtis ~ Treatment * Year"
      )
    }
    
    disp_reg <- run_dispersion_analysis(
      dist_mat = dist_reg,
      group_factor = interaction(meta_reg$Treatment, meta_reg$Year, drop = TRUE),
      meta = meta_reg,
      marker = marker,
      prefix = prefix,
      region = reg_i,
      model_scope = "Region-specific model",
      group_definition = "Treatment × Year"
    )
    
    dispersion_test_list[[length(dispersion_test_list) + 1]] <- disp_reg$test
    dispersion_distance_list[[length(dispersion_distance_list) + 1]] <- disp_reg$distances
    

    for (year_i in sort(unique(meta_reg$Year_num))) {
      
      idx_y <- meta_reg$Year_num == year_i
      
      meta_y <- meta_reg[idx_y, , drop = FALSE]
      rownames(meta_y) <- meta_y$SampleID
      
      genus_rel_y <- genus_rel_reg[rownames(meta_y), , drop = FALSE]
      
      n_control <- sum(meta_y$Treatment == "Control")
      n_wildfire <- sum(meta_y$Treatment == "Wildfire")
      
      if (nrow(meta_y) < 4 || length(unique(meta_y$Treatment)) < 2 || min(table(meta_y$Treatment)) < 2) {
        pairwise_list[[length(pairwise_list) + 1]] <- data.frame(
          Marker = marker,
          Prefix = prefix,
          Region = reg_i,
          Year = as.character(year_i),
          Comparison = "Control vs Wildfire",
          Test = "Pairwise PERMANOVA",
          Distance = "Bray-Curtis",
          Term_found = NA_character_,
          Df = NA_real_,
          SumOfSqs = NA_real_,
          R2 = NA_real_,
          F_value = NA_real_,
          p_value = NA_real_,
          p_value_formatted = NA_character_,
          Significance = "",
          Permutations = PERMUTATIONS,
          n_Control = n_control,
          n_Wildfire = n_wildfire,
          Status = "Not evaluable: insufficient samples",
          stringsAsFactors = FALSE
        )
        next
      }
      
      dist_y <- vegan::vegdist(genus_rel_y, method = "bray")
      
      ad_y <- tryCatch(
        vegan::adonis2(
          dist_y ~ Treatment,
          data = meta_y,
          permutations = PERMUTATIONS
        ),
        error = function(e) NULL
      )
      
      if (is.null(ad_y)) {
        pairwise_list[[length(pairwise_list) + 1]] <- data.frame(
          Marker = marker,
          Prefix = prefix,
          Region = reg_i,
          Year = as.character(year_i),
          Comparison = "Control vs Wildfire",
          Test = "Pairwise PERMANOVA",
          Distance = "Bray-Curtis",
          Term_found = NA_character_,
          Df = NA_real_,
          SumOfSqs = NA_real_,
          R2 = NA_real_,
          F_value = NA_real_,
          p_value = NA_real_,
          p_value_formatted = NA_character_,
          Significance = "",
          Permutations = PERMUTATIONS,
          n_Control = n_control,
          n_Wildfire = n_wildfire,
          Status = "Error in pairwise PERMANOVA",
          stringsAsFactors = FALSE
        )
        next
      }
      
      extracted <- extract_adonis_row(ad_y, term_pattern = "Treatment")
      
      pairwise_list[[length(pairwise_list) + 1]] <- data.frame(
        Marker = marker,
        Prefix = prefix,
        Region = reg_i,
        Year = as.character(year_i),
        Comparison = "Control vs Wildfire",
        Test = "Pairwise PERMANOVA",
        Distance = "Bray-Curtis",
        Term_found = extracted$Term_found,
        Df = extracted$Df,
        SumOfSqs = extracted$SumOfSqs,
        R2 = extracted$R2,
        F_value = extracted$F_value,
        p_value = extracted$p_value,
        p_value_formatted = p_format(extracted$p_value),
        Significance = p_stars(extracted$p_value),
        Permutations = PERMUTATIONS,
        n_Control = n_control,
        n_Wildfire = n_wildfire,
        Status = ifelse(extracted$Status_extract == "OK", "OK", extracted$Status_extract),
        stringsAsFactors = FALSE
      )
    }
  }
  
  list(
    pcoa_coords = pcoa_all$coords,
    pcoa_variance = pcoa_all$variance,
    permanova = bind_rows(permanova_list),
    dispersion_tests = bind_rows(dispersion_test_list),
    dispersion_distances = bind_rows(dispersion_distance_list),
    pairwise = bind_rows(pairwise_list)
  )
}

beta_bac <- run_beta_marker(bac)
beta_arq <- run_beta_marker(arq)
beta_its <- run_beta_marker(its)

pcoa_coords_all <- bind_rows(beta_bac$pcoa_coords, beta_arq$pcoa_coords, beta_its$pcoa_coords)
pcoa_variance_all <- bind_rows(beta_bac$pcoa_variance, beta_arq$pcoa_variance, beta_its$pcoa_variance)

permanova_full <- bind_rows(beta_bac$permanova, beta_arq$permanova, beta_its$permanova) %>%
  mutate(
    p_value_formatted = p_format(p_value),
    Significance = p_stars(p_value)
  ) %>%
  arrange(Prefix, Model_scope, Region, Term)

dispersion_tests <- bind_rows(beta_bac$dispersion_tests, beta_arq$dispersion_tests, beta_its$dispersion_tests) %>%
  arrange(Prefix, Model_scope, Region)

dispersion_distances <- bind_rows(beta_bac$dispersion_distances, beta_arq$dispersion_distances, beta_its$dispersion_distances)

pairwise_permanova <- bind_rows(beta_bac$pairwise, beta_arq$pairwise, beta_its$pairwise) %>%
  group_by(Marker, Prefix, Region) %>%
  mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH"),
    p_adjusted_BH_formatted = p_format(p_adjusted_BH),
    Significance_BH = p_stars(p_adjusted_BH)
  ) %>%
  ungroup() %>%
  arrange(Prefix, Region, as.numeric(Year))


table_s1_qc_dispersion <- dispersion_tests %>%
  left_join(
    genus_qc_marker_all %>%
      select(
        Marker,
        Prefix,
        Total_samples_initial,
        Samples_used_beta_diversity,
        Samples_excluded_no_identified_genus,
        ASVs_with_identified_genus,
        Identified_genera_used,
        Mean_percent_reads_with_identified_genus
      ),
    by = c("Marker", "Prefix")
  ) %>%
  transmute(
    Marker,
    Region,
    Model_scope,
    `Group used for dispersion test` = Group_definition,
    `Initial samples` = Total_samples_initial,
    `Samples used in beta diversity` = Samples_used_beta_diversity,
    `Samples excluded due to no identified genus` = Samples_excluded_no_identified_genus,
    `ASVs with identified genus` = ASVs_with_identified_genus,
    `Identified genera used` = Identified_genera_used,
    `Mean reads assigned to identified genus (%)` = round(Mean_percent_reads_with_identified_genus, 2),
    `Dispersion test` = Dispersion_test,
    `Dispersion F` = round(F_value, 4),
    `Dispersion p value` = p_value_formatted,
    `Dispersion result` = Dispersion_result,
    Permutations
  )

table_s2_permanova <- permanova_full %>%
  transmute(
    Marker,
    Region,
    `Model scope` = Model_scope,
    `Model formula` = Model_formula,
    Distance,
    Term,
    Df,
    `Sum of squares` = round(SumOfSqs, 4),
    R2 = round(R2, 4),
    `F value` = round(F_value, 4),
    `p value` = p_value_formatted,
    Significance,
    Permutations,
    Status
  )

table_s3_pairwise <- pairwise_permanova %>%
  transmute(
    Marker,
    Region,
    Year,
    Comparison,
    Test,
    Distance,
    `n Control` = n_Control,
    `n Wildfire` = n_Wildfire,
    Df,
    `Sum of squares` = round(SumOfSqs, 4),
    R2 = round(R2, 4),
    `F value` = round(F_value, 4),
    `p value` = p_value_formatted,
    `p adjusted BH` = p_adjusted_BH_formatted,
    Significance = Significance_BH,
    Status
  )


marker_titles <- c(
  "BAC" = "Bacterial communities",
  "ARQ" = "Archaeal communities",
  "ITS" = "Fungal communities"
)

treatment_colors <- c(
  "Control" = "#9BD7E8",
  "Wildfire" = "#F4A261"
)

year_shape_values <- c(
  "1" = 16,
  "2" = 17,
  "4" = 15,
  "5" = 3,
  "6" = 7,
  "10" = 8,
  "15" = 18,
  "20" = 4
)

make_pcoa_plot <- function(prefix_i, labelled = FALSE) {
  
  df <- pcoa_coords_all %>%
    filter(Prefix == prefix_i) %>%
    mutate(
      Treatment = factor(Treatment, levels = c("Control", "Wildfire")),
      Region = factor(Region, levels = c("Amazonas", "Cusco")),
      Year = as.character(Year)
    )
  
  var_i <- pcoa_variance_all %>%
    filter(Prefix == prefix_i) %>%
    slice(1)
  
  x_lab <- paste0("PCoA1 (", var_i$PCoA1_percent, "%)")
  y_lab <- paste0("PCoA2 (", var_i$PCoA2_percent, "%)")
  
  shapes_i <- year_shape_values[names(year_shape_values) %in% unique(df$Year)]
  
  p <- ggplot(
    df,
    aes(
      x = PCoA1,
      y = PCoA2,
      color = Treatment,
      shape = Year
    )
  ) +
    stat_ellipse(
      aes(group = Treatment, color = Treatment),
      type = "t",
      level = 0.68,
      linewidth = 0.55,
      linetype = "dashed",
      show.legend = FALSE
    ) +
    geom_point(size = 3.2, alpha = 0.95) +
    facet_wrap(~ Region, nrow = 1) +
    scale_color_manual(values = treatment_colors, drop = FALSE) +
    scale_shape_manual(values = shapes_i, drop = FALSE) +
    labs(
      title = marker_titles[[prefix_i]],
      x = x_lab,
      y = y_lab,
      color = NULL,
      shape = "Post-fire year"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10, color = "black"),
      legend.position = "top",
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.55)
    )
  
  if (labelled) {
    p <- p +
      ggrepel::geom_text_repel(
        aes(label = SampleLabel),
        size = 2.4,
        max.overlaps = 40,
        show.legend = FALSE
      )
  }
  
  p
}

save_plot_all_formats <- function(plot_obj, base_file, width, height) {
  ggsave(paste0(base_file, ".png"),  plot = plot_obj, width = width, height = height, dpi = 400)
  ggsave(paste0(base_file, ".tiff"), plot = plot_obj, width = width, height = height, dpi = 400, compression = "lzw")
  ggsave(paste0(base_file, ".pdf"),  plot = plot_obj, width = width, height = height)
  ggsave(paste0(base_file, ".svg"),  plot = plot_obj, width = width, height = height)
}

p_bac_clean <- make_pcoa_plot("BAC", labelled = FALSE)
p_arq_clean <- make_pcoa_plot("ARQ", labelled = FALSE)
p_its_clean <- make_pcoa_plot("ITS", labelled = FALSE)

p_bac_label <- make_pcoa_plot("BAC", labelled = TRUE)
p_arq_label <- make_pcoa_plot("ARQ", labelled = TRUE)
p_its_label <- make_pcoa_plot("ITS", labelled = TRUE)

save_plot_all_formats(
  p_bac_clean,
  file.path(fig_dir, "BAC_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN"),
  width = 9,
  height = 5.2
)

save_plot_all_formats(
  p_arq_clean,
  file.path(fig_dir, "ARQ_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN"),
  width = 9,
  height = 5.2
)

save_plot_all_formats(
  p_its_clean,
  file.path(fig_dir, "ITS_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN"),
  width = 9,
  height = 5.2
)

save_plot_all_formats(
  p_bac_label,
  file.path(supp_dir, "BAC_beta_diversity_GENUS_BrayCurtis_PCoA_LABELLED"),
  width = 10,
  height = 5.8
)

save_plot_all_formats(
  p_arq_label,
  file.path(supp_dir, "ARQ_beta_diversity_GENUS_BrayCurtis_PCoA_LABELLED"),
  width = 10,
  height = 5.8
)

save_plot_all_formats(
  p_its_label,
  file.path(supp_dir, "ITS_beta_diversity_GENUS_BrayCurtis_PCoA_LABELLED"),
  width = 10,
  height = 5.8
)

p_all_clean <- (p_bac_clean / p_arq_clean / p_its_clean) +
  plot_annotation(
    title = "Genus-level beta diversity based on Bray-Curtis dissimilarity",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  ) &
  theme(legend.position = "top")

save_plot_all_formats(
  p_all_clean,
  file.path(fig_dir, "ALL_MARKERS_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN"),
  width = 10,
  height = 14.5
)

dispersion_plot <- dispersion_distances %>%
  mutate(
    Marker = factor(Marker, levels = c("Bacteria 16S", "Archaea 16S", "Fungi ITS")),
    Treatment = factor(Treatment, levels = c("Control", "Wildfire"))
  ) %>%
  ggplot(
    aes(
      x = Group,
      y = Distance_to_centroid,
      fill = Treatment
    )
  ) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, color = "gray25") +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.85) +
  facet_grid(Marker ~ Model_scope, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = treatment_colors, drop = FALSE) +
  labs(
    title = "Homogeneity of multivariate dispersion",
    x = "Group",
    y = "Distance to centroid",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

save_plot_all_formats(
  dispersion_plot,
  file.path(supp_dir, "S_beta_diversity_GENUS_betadisper_distances_to_centroid"),
  width = 14,
  height = 9
)


write.csv(genus_qc_marker_all, file.path(tab_dir, "01_genus_level_QC_by_marker.csv"), row.names = FALSE)
write.csv(sample_qc_all, file.path(tab_dir, "02_genus_level_QC_by_sample.csv"), row.names = FALSE)
write.csv(pcoa_coords_all, file.path(tab_dir, "03_PCoA_coordinates_GENUS_BrayCurtis.csv"), row.names = FALSE)
write.csv(pcoa_variance_all, file.path(tab_dir, "04_PCoA_variance_explained_GENUS_BrayCurtis.csv"), row.names = FALSE)
write.csv(permanova_full, file.path(tab_dir, "05_PERMANOVA_full_GENUS_BrayCurtis.csv"), row.names = FALSE)
write.csv(dispersion_tests, file.path(tab_dir, "06_betadisper_tests_GENUS_BrayCurtis.csv"), row.names = FALSE)
write.csv(dispersion_distances, file.path(tab_dir, "07_betadisper_distances_to_centroid.csv"), row.names = FALSE)
write.csv(pairwise_permanova, file.path(tab_dir, "08_pairwise_PERMANOVA_Control_vs_Wildfire_by_year.csv"), row.names = FALSE)

write.csv(table_s1_qc_dispersion, file.path(supp_dir, "Table_S1_QC_and_dispersion.csv"), row.names = FALSE)
write.csv(table_s2_permanova, file.path(supp_dir, "Table_S2_PERMANOVA_global_tests.csv"), row.names = FALSE)
write.csv(table_s3_pairwise, file.path(supp_dir, "Table_S3_pairwise_PERMANOVA_by_year.csv"), row.names = FALSE)


xlsx_file <- file.path(
  excel_dir,
  "BETA_DIVERSITY_GENUS_TABLES_FOR_PAPER_ENGLISH.xlsx"
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
  
  if (nrow(data) > 0 && ncol(data) > 0) {
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
}

readme <- data.frame(
  Section = c(
    "Analysis level",
    "Distance method",
    "Ordination",
    "Overall model",
    "Region-specific model",
    "Pairwise comparisons",
    "Dispersion assumption",
    "Recommended main figure",
    "Recommended supplementary tables"
  ),
  Description = c(
    "All beta diversity analyses were performed at the genus level using only identified genera.",
    "Bray-Curtis dissimilarity was calculated from genus-level relative abundance matrices.",
    "PCoA was used for visualization; Cailliez correction was applied during PCoA computation.",
    "The overall marker-level model tested Region, Treatment, and Region × Treatment.",
    "Region-specific models tested Treatment, post-fire Year, and Treatment × Year.",
    "Pairwise PERMANOVA compared Control vs Wildfire soils within each post-fire year.",
    "Homogeneity of multivariate dispersion was evaluated using betadisper and permutest.",
    "Use ALL_MARKERS_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN.tiff or individual marker CLEAN figures.",
    "Use Table S1, Table S2, and Table S3 as compact supplementary statistical support."
  )
)

add_sheet(wb, "00_Read_me", readme)
add_sheet(wb, "Table_S1_QC_dispersion", table_s1_qc_dispersion)
add_sheet(wb, "Table_S2_PERMANOVA", table_s2_permanova)
add_sheet(wb, "Table_S3_Pairwise", table_s3_pairwise)
add_sheet(wb, "PCoA_coordinates", pcoa_coords_all)
add_sheet(wb, "PCoA_variance", pcoa_variance_all)
add_sheet(wb, "Genus_QC_marker", genus_qc_marker_all)
add_sheet(wb, "Genus_QC_sample", sample_qc_all)
add_sheet(wb, "Dispersion_distances", dispersion_distances)

saveWorkbook(wb, xlsx_file, overwrite = TRUE)


methods_text <- c(
  "Beta diversity was assessed at the genus level using Bray-Curtis dissimilarity calculated from relative abundance matrices of identified genera. ASV counts were first aggregated by genus, excluding unclassified, uncultured, unidentified, unknown, ambiguous, chloroplast and mitochondrial assignments. Samples without reads assigned to identified genera were excluded from genus-level beta diversity analyses.",
  "",
  "Principal coordinate analysis (PCoA) was used to visualize community dissimilarities. Differences in genus-level community composition were tested using PERMANOVA with the adonis2 function in vegan, using 999 permutations. Because post-fire ages differed between regions, region-specific models were fitted to test the effects of treatment, post-fire year and their interaction. An overall marker-level model was also used to evaluate broad effects of region, treatment and their interaction. Homogeneity of multivariate dispersion was assessed using betadisper followed by permutest. Pairwise PERMANOVA comparisons between control and wildfire-affected soils were performed within each post-fire year, and p-values were adjusted using the Benjamini-Hochberg method."
)

caption_text <- c(
  "Figure X. Genus-level beta diversity of soil microbial communities based on Bray-Curtis dissimilarity and visualized using principal coordinate analysis (PCoA). Points represent individual soil samples, colors indicate treatment, and shapes indicate post-fire year. Dashed ellipses represent 68% confidence regions for treatment groups within each region."
)

writeLines(methods_text, file.path(text_dir, "METHODS_beta_diversity_GENUS.txt"))
writeLines(caption_text, file.path(text_dir, "FIGURE_CAPTION_beta_diversity_GENUS.txt"))


stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

zip_file <- file.path(
  project_root,
  paste0("BETA_DIVERSITY_GENUS_FINAL_FOR_PAPER_", stamp, ".zip")
)

old_wd <- getwd()
setwd(project_root)

zip::zipr(
  zipfile = zip_file,
  files = basename(outdir)
)

setwd(old_wd)


cat("\n=====================================================\n")
cat("DONE: GENUS-LEVEL BETA DIVERSITY ANALYSIS GENERATED\n")
cat("=====================================================\n")

cat("\nMain output folder:\n")
cat(outdir, "\n")

cat("\nFigures folder:\n")
cat(fig_dir, "\n")

cat("\nExcel file:\n")
cat(xlsx_file, "\n")

cat("\nText for manuscript:\n")
cat(text_dir, "\n")

cat("\nZIP final package:\n")
cat(zip_file, "\n")

cat("\nDownload this ZIP from RStudio:\n")
cat(basename(zip_file), "\n")

cat("\nMain figures generated:\n")
cat("ALL_MARKERS_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN.tiff\n")
cat("BAC_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN.tiff\n")
cat("ARQ_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN.tiff\n")
cat("ITS_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN.tiff\n")

cat("\nCompact supplementary tables:\n")
cat("Table S1: QC and dispersion tests\n")
cat("Table S2: PERMANOVA global tests\n")
cat("Table S3: Pairwise PERMANOVA by post-fire year\n")

cat("\nPERMANOVA summary:\n")
print(
  permanova_full %>%
    count(Marker, Prefix, Region, Model_scope, Term)
)

cat("\nDispersion summary:\n")
print(
  dispersion_tests %>%
    select(Marker, Region, Model_scope, Group_definition, p_value_formatted, Dispersion_result)
)

cat("\nPairwise PERMANOVA summary:\n")
print(
  pairwise_permanova %>%
    count(Marker, Prefix, Region, Status)
)

cat("\nIMPORTANT:\n")
cat("The warning about 'complete enumeration' is not an error.\n")
cat("It appears because some pairwise comparisons have few samples and vegan uses all possible permutations.\n")

cat("=====================================================\n")


rm(list = ls())
options(stringsAsFactors = FALSE)


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(vegan)
  library(ape)
  library(openxlsx)
  library(svglite)
  library(ggrepel)
  library(zip)
})

cat("\n=====================================================\n")
cat("REPLACING BETA DIVERSITY OUTPUTS WITH REGION-SPECIFIC PCoA\n")
cat("GENUS LEVEL | BRAY-CURTIS | PCoA | TABLE S4 UPDATED\n")
cat("=====================================================\n")


project_root <- path.expand("~/Rosmery/Tesis Rosmery-Metabarcoding")

beta_dir <- file.path(
  project_root,
  "RESULTADOS_BETA_DIVERSITY_GENUS_FINAL_PAPER"
)

matrix_dir <- file.path(beta_dir, "06_GENUS_MATRICES")
tab_dir    <- file.path(beta_dir, "02_TABLES_CSV")
fig_dir    <- file.path(beta_dir, "01_FIGURES")
supp_dir   <- file.path(beta_dir, "04_SUPPLEMENTARY")
excel_dir  <- file.path(beta_dir, "03_EXCEL_TABLES")
text_dir   <- file.path(beta_dir, "05_TEXT_FOR_MANUSCRIPT")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(excel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)

sample_qc_file <- file.path(tab_dir, "02_genus_level_QC_by_sample.csv")
genus_qc_file  <- file.path(tab_dir, "01_genus_level_QC_by_marker.csv")
disp_file      <- file.path(tab_dir, "06_betadisper_tests_GENUS_BrayCurtis.csv")
permanova_file <- file.path(tab_dir, "05_PERMANOVA_full_GENUS_BrayCurtis.csv")
pairwise_file  <- file.path(tab_dir, "08_pairwise_PERMANOVA_Control_vs_Wildfire_by_year.csv")

required_files <- c(
  sample_qc_file,
  genus_qc_file,
  disp_file,
  permanova_file,
  pairwise_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "These required files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

sample_qc  <- read.csv(sample_qc_file, check.names = FALSE)
genus_qc   <- read.csv(genus_qc_file, check.names = FALSE)
disp_tests <- read.csv(disp_file, check.names = FALSE)
permanova_full <- read.csv(permanova_file, check.names = FALSE)
pairwise_permanova <- read.csv(pairwise_file, check.names = FALSE)


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

clean_sample_label <- function(x) {
  x <- as.character(x)
  x <- sub("\\.raw$", "", x)
  x <- sub("_(BAC|ARQ|FUN|ITS)$", "", x)
  x <- sub("_(BAC|ARQ|FUN|ITS)\\.raw$", "", x)
  x <- gsub("_", "", x)
  x
}

save_plot_all_formats <- function(plot_obj, base_file, width, height) {
  
  ggsave(
    paste0(base_file, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 500
  )
  
  ggsave(
    paste0(base_file, ".tiff"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 500,
    compression = "lzw"
  )
  
  ggsave(
    paste0(base_file, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height
  )
  
  ggsave(
    paste0(base_file, ".svg"),
    plot = plot_obj,
    width = width,
    height = height
  )
}


marker_info <- data.frame(
  Prefix = c("BAC", "ARQ", "ITS"),
  Marker = c("Bacteria 16S", "Archaea 16S", "Fungi ITS"),
  Matrix_file = c(
    "BAC_genus_relative_abundance_sample_by_genus.csv",
    "ARQ_genus_relative_abundance_sample_by_genus.csv",
    "ITS_genus_relative_abundance_sample_by_genus.csv"
  ),
  Plot_title = c(
    "Bacterial communities",
    "Archaeal communities",
    "Fungal communities"
  ),
  stringsAsFactors = FALSE
)

treatment_colors <- c(
  "Control" = "#9BD7E8",
  "Wildfire" = "#F4A261"
)

year_shape_values <- c(
  "1" = 16,
  "2" = 17,
  "4" = 15,
  "5" = 3,
  "6" = 7,
  "10" = 8,
  "15" = 18,
  "20" = 4
)


if (!"p_value_formatted" %in% colnames(disp_tests)) {
  disp_tests$p_value_formatted <- p_format(disp_tests$p_value)
}

if (!"Significance" %in% colnames(permanova_full)) {
  permanova_full$Significance <- p_stars(permanova_full$p_value)
}

if (!"p_value_formatted" %in% colnames(permanova_full)) {
  permanova_full$p_value_formatted <- p_format(permanova_full$p_value)
}

if (!"p_adjusted_BH" %in% colnames(pairwise_permanova)) {
  pairwise_permanova <- pairwise_permanova %>%
    group_by(Marker, Prefix, Region) %>%
    mutate(
      p_adjusted_BH = p.adjust(p_value, method = "BH")
    ) %>%
    ungroup()
}

pairwise_permanova <- pairwise_permanova %>%
  mutate(
    p_value_formatted = p_format(p_value),
    p_adjusted_BH_formatted = p_format(p_adjusted_BH),
    Significance_BH = p_stars(p_adjusted_BH)
  )


run_region_specific_pcoa <- function(prefix_i, region_i) {
  
  info_i <- marker_info %>% filter(Prefix == prefix_i)
  
  mat_file <- file.path(matrix_dir, info_i$Matrix_file)
  
  if (!file.exists(mat_file)) {
    stop("Matrix file not found: ", mat_file)
  }
  
  rel_df <- read.csv(mat_file, check.names = FALSE)
  
  if (!"SampleID" %in% colnames(rel_df)) {
    stop("SampleID column not found in: ", mat_file)
  }
  
  meta_i <- sample_qc %>%
    filter(Prefix == prefix_i) %>%
    mutate(
      SampleID = as.character(SampleID),
      SampleLabel = clean_sample_label(SampleID),
      Region = as.character(Region),
      Treatment = as.character(Treatment),
      Year = as.character(Year),
      Year_num = suppressWarnings(as.numeric(Year))
    ) %>%
    filter(
      Region == region_i,
      Beta_diversity_status == "Included"
    ) %>%
    select(SampleID, SampleLabel, Region, Treatment, Year, Year_num)
  
  rel_region <- rel_df %>%
    mutate(SampleID = as.character(SampleID)) %>%
    inner_join(meta_i, by = "SampleID")
  
  if (nrow(rel_region) < 4) {
    stop("Too few samples for region-specific PCoA: ", prefix_i, " - ", region_i)
  }
  
  genus_cols <- setdiff(
    colnames(rel_region),
    c("SampleID", "SampleLabel", "Region", "Treatment", "Year", "Year_num")
  )
  
  genus_mat <- rel_region[, genus_cols, drop = FALSE] %>%
    as.matrix()
  
  rownames(genus_mat) <- rel_region$SampleID
  
  genus_mat[is.na(genus_mat)] <- 0
  
  genus_mat <- genus_mat[
    rowSums(genus_mat, na.rm = TRUE) > 0,
    colSums(genus_mat, na.rm = TRUE) > 0,
    drop = FALSE
  ]
  
  meta_used <- rel_region %>%
    filter(SampleID %in% rownames(genus_mat)) %>%
    select(SampleID, SampleLabel, Region, Treatment, Year, Year_num)
  
  rownames(meta_used) <- meta_used$SampleID
  
  genus_mat <- genus_mat[meta_used$SampleID, , drop = FALSE]
  
  dist_i <- vegan::vegdist(genus_mat, method = "bray")
  
  pco_i <- tryCatch(
    ape::pcoa(dist_i, correction = "cailliez"),
    error = function(e) ape::pcoa(dist_i)
  )
  
  rel_eig <- pco_i$values$Relative_eig
  
  if (length(rel_eig) < 2 || all(is.na(rel_eig))) {
    eig <- pco_i$values$Eigenvalues
    eig_pos <- eig[eig > 0]
    rel_eig <- eig / sum(eig_pos)
  }
  
  pcoa1_percent <- round(rel_eig[1] * 100, 2)
  pcoa2_percent <- round(rel_eig[2] * 100, 2)
  
  coords <- data.frame(
    SampleID = rownames(pco_i$vectors),
    PCoA1 = pco_i$vectors[, 1],
    PCoA2 = pco_i$vectors[, 2],
    stringsAsFactors = FALSE
  ) %>%
    left_join(meta_used, by = "SampleID") %>%
    mutate(
      Marker = info_i$Marker,
      Prefix = prefix_i,
      PCoA_scope = "Region-specific PCoA",
      PCoA1_percent = pcoa1_percent,
      PCoA2_percent = pcoa2_percent
    )
  
  variance <- data.frame(
    Marker = info_i$Marker,
    Prefix = prefix_i,
    Region = region_i,
    PCoA_scope = "Region-specific PCoA",
    n_samples = nrow(coords),
    PCoA1_percent = pcoa1_percent,
    PCoA2_percent = pcoa2_percent,
    Bray_Curtis_matrix = "Calculated independently within each marker and region",
    stringsAsFactors = FALSE
  )
  
  list(
    coords = coords,
    variance = variance
  )
}


make_pcoa_plot <- function(coords_i, labelled = FALSE) {
  
  prefix_i <- unique(coords_i$Prefix)
  region_i <- unique(coords_i$Region)
  
  info_i <- marker_info %>% filter(Prefix == prefix_i)
  
  coords_i <- coords_i %>%
    mutate(
      Treatment = factor(Treatment, levels = c("Control", "Wildfire")),
      Year = as.character(Year)
    )
  
  shapes_i <- year_shape_values[names(year_shape_values) %in% unique(coords_i$Year)]
  
  x_lab <- paste0("PCoA1 (", unique(coords_i$PCoA1_percent), "%)")
  y_lab <- paste0("PCoA2 (", unique(coords_i$PCoA2_percent), "%)")
  
  p <- ggplot(
    coords_i,
    aes(
      x = PCoA1,
      y = PCoA2,
      color = Treatment,
      shape = Year
    )
  ) +
    stat_ellipse(
      aes(group = Treatment, color = Treatment),
      type = "t",
      level = 0.68,
      linewidth = 0.55,
      linetype = "dashed",
      show.legend = FALSE
    ) +
    geom_point(
      size = 3.3,
      alpha = 0.95
    ) +
    scale_color_manual(values = treatment_colors, drop = FALSE) +
    scale_shape_manual(values = shapes_i, drop = FALSE) +
    labs(
      title = paste0(info_i$Plot_title, " - ", region_i),
      x = x_lab,
      y = y_lab,
      color = NULL,
      shape = "Post-fire year"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10.5, face = "bold"),
      axis.text = element_text(size = 9, color = "black"),
      legend.position = "bottom",
      legend.text = element_text(size = 8.5),
      legend.title = element_text(size = 8.5, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.55)
    )
  
  if (labelled) {
    p <- p +
      ggrepel::geom_text_repel(
        aes(label = SampleLabel),
        size = 2.5,
        max.overlaps = 60,
        show.legend = FALSE
      )
  }
  
  return(p)
}


results_list <- list()
plot_clean_list <- list()
plot_label_list <- list()

for (prefix_i in c("BAC", "ARQ", "ITS")) {
  
  for (region_i in c("Amazonas", "Cusco")) {
    
    cat("\nRunning region-specific PCoA:", prefix_i, "-", region_i, "\n")
    
    key_i <- paste(prefix_i, region_i, sep = "_")
    
    res_i <- run_region_specific_pcoa(prefix_i, region_i)
    
    results_list[[key_i]] <- res_i
    
    plot_clean_list[[key_i]] <- make_pcoa_plot(res_i$coords, labelled = FALSE)
    plot_label_list[[key_i]] <- make_pcoa_plot(res_i$coords, labelled = TRUE)
  }
}

pcoa_coords_region_specific <- bind_rows(lapply(results_list, function(x) x$coords))
pcoa_variance_region_specific <- bind_rows(lapply(results_list, function(x) x$variance))


p_bac_clean <- plot_clean_list$BAC_Amazonas + plot_clean_list$BAC_Cusco +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Bacterial communities",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  ) &
  theme(legend.position = "bottom")

p_arq_clean <- plot_clean_list$ARQ_Amazonas + plot_clean_list$ARQ_Cusco +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Archaeal communities",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  ) &
  theme(legend.position = "bottom")

p_its_clean <- plot_clean_list$ITS_Amazonas + plot_clean_list$ITS_Cusco +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Fungal communities",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  ) &
  theme(legend.position = "bottom")

p_all_clean <-
  (plot_clean_list$BAC_Amazonas + plot_clean_list$BAC_Cusco) /
  (plot_clean_list$ARQ_Amazonas + plot_clean_list$ARQ_Cusco) /
  (plot_clean_list$ITS_Amazonas + plot_clean_list$ITS_Cusco) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Region-specific genus-level beta diversity based on Bray-Curtis dissimilarity",
    theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5))
  ) &
  theme(legend.position = "bottom")


save_plot_all_formats(
  p_bac_clean,
  file.path(fig_dir, "BAC_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN"),
  width = 11.5,
  height = 5.3
)

save_plot_all_formats(
  p_arq_clean,
  file.path(fig_dir, "ARQ_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN"),
  width = 11.5,
  height = 5.3
)

save_plot_all_formats(
  p_its_clean,
  file.path(fig_dir, "ITS_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN"),
  width = 11.5,
  height = 5.3
)

save_plot_all_formats(
  p_all_clean,
  file.path(fig_dir, "ALL_MARKERS_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN"),
  width = 11.8,
  height = 13.8
)


p_bac_label <- plot_label_list$BAC_Amazonas + plot_label_list$BAC_Cusco +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Bacterial communities",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  ) &
  theme(legend.position = "bottom")

p_arq_label <- plot_label_list$ARQ_Amazonas + plot_label_list$ARQ_Cusco +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Archaeal communities",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  ) &
  theme(legend.position = "bottom")

p_its_label <- plot_label_list$ITS_Amazonas + plot_label_list$ITS_Cusco +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Fungal communities",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  ) &
  theme(legend.position = "bottom")

p_all_label <-
  (plot_label_list$BAC_Amazonas + plot_label_list$BAC_Cusco) /
  (plot_label_list$ARQ_Amazonas + plot_label_list$ARQ_Cusco) /
  (plot_label_list$ITS_Amazonas + plot_label_list$ITS_Cusco) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Region-specific genus-level beta diversity based on Bray-Curtis dissimilarity",
    theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5))
  ) &
  theme(legend.position = "bottom")

save_plot_all_formats(
  p_bac_label,
  file.path(supp_dir, "BAC_beta_diversity_GENUS_BrayCurtis_PCoA_LABELLED"),
  width = 12.2,
  height = 5.8
)

save_plot_all_formats(
  p_arq_label,
  file.path(supp_dir, "ARQ_beta_diversity_GENUS_BrayCurtis_PCoA_LABELLED"),
  width = 12.2,
  height = 5.8
)

save_plot_all_formats(
  p_its_label,
  file.path(supp_dir, "ITS_beta_diversity_GENUS_BrayCurtis_PCoA_LABELLED"),
  width = 12.2,
  height = 5.8
)

save_plot_all_formats(
  p_all_label,
  file.path(supp_dir, "ALL_MARKERS_beta_diversity_GENUS_BrayCurtis_PCoA_LABELLED"),
  width = 12.6,
  height = 14.8
)


sample_counts_region <- sample_qc %>%
  filter(Beta_diversity_status == "Included") %>%
  group_by(Marker, Prefix, Region) %>%
  summarise(
    `Samples used in beta diversity` = n(),
    .groups = "drop"
  )

disp_region_specific <- disp_tests %>%
  filter(Model_scope == "Region-specific model") %>%
  select(
    Marker,
    Prefix,
    Region,
    Model_scope,
    Group_definition,
    Dispersion_test,
    F_value,
    p_value_formatted,
    Dispersion_result,
    Permutations
  )

table_s4_updated <- pcoa_variance_region_specific %>%
  left_join(
    genus_qc %>%
      select(
        Marker,
        Prefix,
        Total_samples_initial,
        Samples_excluded_no_identified_genus,
        ASVs_with_identified_genus,
        Identified_genera_used,
        Mean_percent_reads_with_identified_genus
      ),
    by = c("Marker", "Prefix")
  ) %>%
  left_join(
    sample_counts_region,
    by = c("Marker", "Prefix", "Region")
  ) %>%
  left_join(
    disp_region_specific,
    by = c("Marker", "Prefix", "Region")
  ) %>%
  transmute(
    Marker,
    Region,
    `PCoA scope` = PCoA_scope,
    `Samples used in beta diversity` = `Samples used in beta diversity`,
    `Samples excluded due to no identified genus` = Samples_excluded_no_identified_genus,
    `ASVs with identified genus` = ASVs_with_identified_genus,
    `Identified genera used` = Identified_genera_used,
    `Mean reads assigned to identified genus (%)` = round(Mean_percent_reads_with_identified_genus, 2),
    `PCoA1 explained variance (%)` = PCoA1_percent,
    `PCoA2 explained variance (%)` = PCoA2_percent,
    `Dispersion group` = Group_definition,
    `Dispersion test` = Dispersion_test,
    `Dispersion F` = round(F_value, 4),
    `Dispersion p value` = p_value_formatted,
    `Dispersion result` = Dispersion_result,
    Permutations
  ) %>%
  arrange(Marker, Region)


table_s5_permanova <- permanova_full %>%
  transmute(
    Marker,
    Region,
    `Model scope` = Model_scope,
    `Model formula` = Model_formula,
    Distance,
    Term,
    Df,
    `Sum of squares` = round(SumOfSqs, 4),
    R2 = round(R2, 4),
    `F value` = round(F_value, 4),
    `p value` = p_value_formatted,
    Significance,
    Permutations,
    Status
  ) %>%
  arrange(Marker, Region, `Model scope`, Term)

table_s6_pairwise <- pairwise_permanova %>%
  transmute(
    Marker,
    Region,
    Year,
    Comparison,
    Test,
    Distance,
    `n Control` = n_Control,
    `n Wildfire` = n_Wildfire,
    Df,
    `Sum of squares` = round(SumOfSqs, 4),
    R2 = round(R2, 4),
    `F value` = round(F_value, 4),
    `p value` = p_value_formatted,
    `p adjusted BH` = p_adjusted_BH_formatted,
    Significance = Significance_BH,
    Status
  ) %>%
  arrange(Marker, Region, as.numeric(Year))


old_supp_table_files <- c(
  file.path(supp_dir, "Table_S1_QC_and_dispersion.csv"),
  file.path(supp_dir, "Table_S2_PERMANOVA_global_tests.csv"),
  file.path(supp_dir, "Table_S3_pairwise_PERMANOVA_by_year.csv")
)

unlink(old_supp_table_files)

write.csv(
  table_s4_updated,
  file.path(supp_dir, "Table_S4_QC_PCoA_region_specific_and_dispersion.csv"),
  row.names = FALSE
)

write.csv(
  table_s5_permanova,
  file.path(supp_dir, "Table_S5_PERMANOVA_global_tests.csv"),
  row.names = FALSE
)

write.csv(
  table_s6_pairwise,
  file.path(supp_dir, "Table_S6_pairwise_PERMANOVA_by_post_fire_year.csv"),
  row.names = FALSE
)

write.csv(
  pcoa_coords_region_specific,
  file.path(tab_dir, "09_REGION_SPECIFIC_PCoA_coordinates_GENUS_BrayCurtis.csv"),
  row.names = FALSE
)

write.csv(
  pcoa_variance_region_specific,
  file.path(tab_dir, "10_REGION_SPECIFIC_PCoA_variance_GENUS_BrayCurtis.csv"),
  row.names = FALSE
)


xlsx_file <- file.path(
  excel_dir,
  "BETA_DIVERSITY_GENUS_TABLES_FOR_PAPER_ENGLISH.xlsx"
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
  
  if (nrow(data) > 0 && ncol(data) > 0) {
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
}

readme <- data.frame(
  Section = c(
    "Analysis level",
    "Distance method",
    "Updated ordination",
    "What changed",
    "What did not change",
    "Main figure to use",
    "Supplementary tables",
    "Important interpretation note"
  ),
  Description = c(
    "All beta diversity analyses were performed at the genus level using only identified genera.",
    "Bray-Curtis dissimilarity was calculated from genus-level relative abundance matrices.",
    "PCoA ordinations were recalculated independently for each marker and region.",
    "The PCoA figure and the PCoA1/PCoA2 variance values in Table S4 were updated.",
    "PERMANOVA and pairwise PERMANOVA results were not recalculated or changed statistically.",
    "Use ALL_MARKERS_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN.tiff as the main beta diversity figure.",
    "Use Table S4, Table S5, and Table S6 as compact supplementary support.",
    "Because PCoA was computed independently by region, PCoA axes are interpreted within each region and should not be compared directly between Amazonas and Cusco."
  ),
  stringsAsFactors = FALSE
)

add_sheet(wb, "00_Read_me", readme)
add_sheet(wb, "Table_S4_QC_PCoA_disp", table_s4_updated)
add_sheet(wb, "Table_S5_PERMANOVA", table_s5_permanova)
add_sheet(wb, "Table_S6_Pairwise", table_s6_pairwise)
add_sheet(wb, "PCoA_variance_region", pcoa_variance_region_specific)
add_sheet(wb, "PCoA_coordinates_region", pcoa_coords_region_specific)

saveWorkbook(wb, xlsx_file, overwrite = TRUE)


methods_text <- c(
  "Beta diversity was assessed at the genus level using Bray-Curtis dissimilarity calculated from relative abundance matrices of identified genera. ASV counts were first aggregated by genus, excluding unclassified, uncultured, unidentified, unknown, ambiguous, chloroplast and mitochondrial assignments. Samples without reads assigned to identified genera were excluded from genus-level beta diversity analyses.",
  "",
  "Principal coordinate analysis (PCoA) was used to visualize community dissimilarities. For the main visualization, PCoA ordinations were computed independently for each molecular marker and region to emphasize within-region patterns. Differences in genus-level community composition were tested using PERMANOVA with the adonis2 function in vegan, using 999 permutations. Because post-fire ages differed between regions, region-specific models were fitted to test the effects of treatment, post-fire year and their interaction. An overall marker-level model was also used to evaluate broad effects of region, treatment and their interaction. Homogeneity of multivariate dispersion was assessed using betadisper followed by permutest. Pairwise PERMANOVA comparisons between control and wildfire-affected soils were performed within each post-fire year, and p-values were adjusted using the Benjamini-Hochberg method."
)

caption_main <- c(
  "Figure 4. Region-specific genus-level beta diversity of soil microbial communities based on Bray-Curtis dissimilarity and visualized using principal coordinate analysis (PCoA). PCoA ordinations were computed independently for each region and molecular marker to visualize within-region community patterns. Points represent individual soil samples, colors indicate treatment, and shapes indicate post-fire year. Dashed ellipses represent 68% confidence regions for treatment groups."
)

caption_table_s4 <- c(
  "Table S4. Genus-level quality control, region-specific PCoA variance, and multivariate dispersion tests for beta diversity analyses. PCoA1 and PCoA2 percentages correspond to ordinations computed independently for each marker and region. Homogeneity of multivariate dispersion was evaluated using betadisper followed by permutest."
)

caption_table_s5 <- c(
  "Table S5. Global PERMANOVA tests for genus-level beta diversity based on Bray-Curtis dissimilarity. The table reports the effects of region, treatment, post-fire year, and their interactions according to the corresponding model."
)

caption_table_s6 <- c(
  "Table S6. Pairwise PERMANOVA comparisons between control and wildfire-affected soils by post-fire year. Pairwise p-values were adjusted using the Benjamini-Hochberg method."
)

writeLines(methods_text, file.path(text_dir, "METHODS_beta_diversity_GENUS.txt"))
writeLines(caption_main, file.path(text_dir, "FIGURE_CAPTION_beta_diversity_GENUS.txt"))
writeLines(caption_table_s4, file.path(text_dir, "CAPTION_Table_S4_beta_diversity.txt"))
writeLines(caption_table_s5, file.path(text_dir, "CAPTION_Table_S5_beta_diversity.txt"))
writeLines(caption_table_s6, file.path(text_dir, "CAPTION_Table_S6_beta_diversity.txt"))


final_dir <- file.path(
  beta_dir,
  "FINAL_MANUSCRIPT_READY_REPLACED"
)

if (dir.exists(final_dir)) unlink(final_dir, recursive = TRUE, force = TRUE)

final_main_fig_dir <- file.path(final_dir, "MAIN_FIGURE")
final_supp_tab_dir <- file.path(final_dir, "SUPPLEMENTARY_TABLES")
final_supp_fig_dir <- file.path(final_dir, "SUPPLEMENTARY_FIGURES")
final_text_dir     <- file.path(final_dir, "METHODS_AND_CAPTIONS")

dir.create(final_main_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_supp_tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_text_dir, recursive = TRUE, showWarnings = FALSE)


file.copy(
  from = list.files(
    fig_dir,
    pattern = "^ALL_MARKERS_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN\\.(tiff|png|pdf|svg)$",
    full.names = TRUE
  ),
  to = final_main_fig_dir,
  overwrite = TRUE
)


file.copy(
  from = c(
    file.path(supp_dir, "Table_S4_QC_PCoA_region_specific_and_dispersion.csv"),
    file.path(supp_dir, "Table_S5_PERMANOVA_global_tests.csv"),
    file.path(supp_dir, "Table_S6_pairwise_PERMANOVA_by_post_fire_year.csv"),
    xlsx_file
  ),
  to = final_supp_tab_dir,
  overwrite = TRUE
)


heatmap_dir <- file.path(
  supp_dir,
  "BRAY_CURTIS_HEATMAPS_FIXED_SAMPLE_ORDER"
)

if (dir.exists(heatmap_dir)) {
  heatmap_files <- list.files(
    heatmap_dir,
    pattern = "ALL_MARKERS.*HEATMAP.*\\.(tiff|png|pdf)$",
    full.names = TRUE
  )
  file.copy(heatmap_files, final_supp_fig_dir, overwrite = TRUE)
}

betadisper_files <- list.files(
  supp_dir,
  pattern = "^S_beta_diversity_GENUS_betadisper_distances_to_centroid\\.(tiff|png|pdf|svg)$",
  full.names = TRUE
)

if (length(betadisper_files) > 0) {
  file.copy(betadisper_files, final_supp_fig_dir, overwrite = TRUE)
}


labelled_pcoa_files <- list.files(
  supp_dir,
  pattern = "^ALL_MARKERS_beta_diversity_GENUS_BrayCurtis_PCoA_LABELLED\\.(tiff|png|pdf|svg)$",
  full.names = TRUE
)

if (length(labelled_pcoa_files) > 0) {
  file.copy(labelled_pcoa_files, final_supp_fig_dir, overwrite = TRUE)
}


file.copy(
  from = list.files(
    text_dir,
    pattern = "METHODS_beta_diversity_GENUS|FIGURE_CAPTION_beta_diversity_GENUS|CAPTION_Table_S4|CAPTION_Table_S5|CAPTION_Table_S6",
    full.names = TRUE
  ),
  to = final_text_dir,
  overwrite = TRUE
)


stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

zip_file <- file.path(
  beta_dir,
  paste0("BETA_DIVERSITY_GENUS_REGION_SPECIFIC_FINAL_REPLACED_", stamp, ".zip")
)

old_wd <- getwd()
setwd(beta_dir)

zip::zipr(
  zipfile = zip_file,
  files = "FINAL_MANUSCRIPT_READY_REPLACED"
)

setwd(old_wd)


cat("\n=====================================================\n")
cat("DONE: BETA DIVERSITY OUTPUTS REPLACED AND UPDATED\n")
cat("=====================================================\n")

cat("\nMain PCoA figures were overwritten here:\n")
cat(fig_dir, "\n")

cat("\nUse this main figure:\n")
cat("ALL_MARKERS_beta_diversity_GENUS_BrayCurtis_PCoA_CLEAN.tiff\n")

cat("\nUpdated Excel file:\n")
cat(xlsx_file, "\n")

cat("\nUpdated supplementary tables:\n")
cat(file.path(supp_dir, "Table_S4_QC_PCoA_region_specific_and_dispersion.csv"), "\n")
cat(file.path(supp_dir, "Table_S5_PERMANOVA_global_tests.csv"), "\n")
cat(file.path(supp_dir, "Table_S6_pairwise_PERMANOVA_by_post_fire_year.csv"), "\n")

cat("\nClean manuscript-ready folder:\n")
cat(final_dir, "\n")

cat("\nClean ZIP package:\n")
cat(zip_file, "\n")

cat("\nIMPORTANT:\n")
cat("Table S4 was updated because PCoA1 and PCoA2 are now region-specific.\n")
cat("Table S5 and Table S6 were re-exported but their PERMANOVA statistics did not change.\n")
cat("Heatmaps and betadisper figures were not recalculated here; they were copied as supplementary support if available.\n")

cat("\nPCoA variance by Marker × Region:\n")
print(pcoa_variance_region_specific)

cat("=====================================================\n")


rm(list = ls())
options(stringsAsFactors = FALSE)


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(vegan)
  library(zip)
})

cat("\n=====================================================\n")
cat("GENUS-LEVEL BRAY-CURTIS HEATMAPS - FIXED SAMPLE ORDER\n")
cat("=====================================================\n")


project_root <- path.expand("~/Rosmery/Tesis Rosmery-Metabarcoding")

beta_dir <- file.path(
  project_root,
  "RESULTADOS_BETA_DIVERSITY_GENUS_FINAL_PAPER"
)

matrix_dir <- file.path(beta_dir, "06_GENUS_MATRICES")
qc_file <- file.path(beta_dir, "02_TABLES_CSV", "02_genus_level_QC_by_sample.csv")

outdir <- file.path(
  beta_dir,
  "04_SUPPLEMENTARY",
  "BRAY_CURTIS_HEATMAPS_FIXED_SAMPLE_ORDER"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


amazonas_order <- c(
  "ACYU1", "ACYU2", "ACYU3",
  "AYU1",  "AYU2",  "AYU3",
  "ACKU1", "ACKU2", "ACKU3",
  "AKU1",  "AKU2",  "AKU3",
  "ACLE1", "ACLE2", "ACLE3",
  "ALE1",  "ALE2",  "ALE3",
  "ACTI1", "ACTI2", "ACTI3",
  "ATI1",  "ATI2",  "ATI3"
)

cusco_order <- c(
  "CCPAT1", "CCPAT2", "CCPAT3",
  "CPAT1",  "CPAT2",  "CPAT3",
  "CCPUC1", "CCPUC2", "CCPUC3",
  "CPUC1",  "CPUC2",  "CPUC3",
  "CCMAT1", "CCMAT2", "CCMAT3",
  "CMAT1",  "CMAT2",  "CMAT3",
  "CCHUC1", "CCHUC2", "CCHUC3",
  "CHUC1",  "CHUC2",  "CHUC3",
  "CCHAC1", "CCHAC2", "CCHAC3",
  "CHAC1",  "CHAC2",  "CHAC3"
)


clean_sample_label <- function(x) {
  x <- as.character(x)
  x <- sub("\\.raw$", "", x)
  x <- sub("_(BAC|ARQ|FUN|ITS)$", "", x)
  x <- sub("_(BAC|ARQ|FUN|ITS)\\.raw$", "", x)
  x <- gsub("_", "", x)
  x
}

marker_info <- data.frame(
  Prefix = c("BAC", "ARQ", "ITS"),
  Marker = c("Bacteria 16S", "Archaea 16S", "Fungi ITS"),
  Matrix_file = c(
    "BAC_genus_relative_abundance_sample_by_genus.csv",
    "ARQ_genus_relative_abundance_sample_by_genus.csv",
    "ITS_genus_relative_abundance_sample_by_genus.csv"
  ),
  Plot_title = c(
    "Bacterial communities",
    "Archaeal communities",
    "Fungal communities"
  ),
  stringsAsFactors = FALSE
)

save_plot_all_formats <- function(plot_obj, base_file, width, height) {
  
  ggsave(
    paste0(base_file, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 500
  )
  
  ggsave(
    paste0(base_file, ".tiff"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 500,
    compression = "lzw"
  )
  
  ggsave(
    paste0(base_file, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height
  )
}


make_bray_matrix_fixed_order <- function(prefix_i, region_i) {
  
  info_i <- marker_info %>% filter(Prefix == prefix_i)
  
  mat_file <- file.path(matrix_dir, info_i$Matrix_file)
  
  if (!file.exists(mat_file)) {
    stop("Matrix file not found: ", mat_file)
  }
  
  rel_df <- read.csv(mat_file, check.names = FALSE)
  
  if (!"SampleID" %in% colnames(rel_df)) {
    stop("SampleID column not found in: ", mat_file)
  }
  
  rel_df$SampleLabel <- clean_sample_label(rel_df$SampleID)
  

  desired_order <- if (region_i == "Amazonas") amazonas_order else cusco_order
  

  present_labels <- intersect(desired_order, rel_df$SampleLabel)
  
  if (length(present_labels) < 2) {
    warning("Less than two samples found for ", prefix_i, " - ", region_i)
    
    empty_mat <- matrix(
      NA_real_,
      nrow = length(desired_order),
      ncol = length(desired_order),
      dimnames = list(desired_order, desired_order)
    )
    
    return(empty_mat)
  }
  
  rel_region <- rel_df %>%
    filter(SampleLabel %in% present_labels) %>%
    mutate(SampleLabel = factor(SampleLabel, levels = desired_order)) %>%
    arrange(SampleLabel)
  
  row_labels <- as.character(rel_region$SampleLabel)
  
  genus_mat <- rel_region %>%
    select(-SampleID, -SampleLabel) %>%
    as.matrix()
  
  rownames(genus_mat) <- row_labels
  

  genus_mat <- genus_mat[, colSums(genus_mat, na.rm = TRUE) > 0, drop = FALSE]
  

  bray <- as.matrix(vegan::vegdist(genus_mat, method = "bray"))
  


  full_bray <- matrix(
    NA_real_,
    nrow = length(desired_order),
    ncol = length(desired_order),
    dimnames = list(desired_order, desired_order)
  )
  
  full_bray[rownames(bray), colnames(bray)] <- bray
  
  return(full_bray)
}


plot_bray_heatmap <- function(bray_mat,
                              prefix_i,
                              region_i,
                              show_sample_names = TRUE) {
  
  info_i <- marker_info %>% filter(Prefix == prefix_i)
  
  sample_order <- rownames(bray_mat)
  
  df <- as.data.frame(as.table(bray_mat))
  colnames(df) <- c("Sample_row", "Sample_col", "Bray_Curtis")
  
  df$Sample_row <- factor(df$Sample_row, levels = rev(sample_order))
  df$Sample_col <- factor(df$Sample_col, levels = sample_order)
  


  block_lines <- if (region_i == "Amazonas") {
    c(6.5, 12.5, 18.5)
  } else {
    c(6.5, 12.5, 18.5, 24.5)
  }
  
  p <- ggplot(df, aes(x = Sample_col, y = Sample_row, fill = Bray_Curtis)) +
    geom_tile(color = NA) +
    geom_vline(xintercept = block_lines, color = "black", linewidth = 0.55) +
    geom_hline(
      yintercept = length(sample_order) - block_lines + 1,
      color = "black",
      linewidth = 0.55
    ) +
    scale_fill_gradient2(
      low = "#08306B",
      mid = "white",
      high = "#B2182B",
      midpoint = 0.5,
      limits = c(0, 1),
      na.value = "grey80",
      name = "Bray-Curtis\ndissimilarity"
    ) +
    coord_fixed() +
    labs(
      title = paste0(info_i$Plot_title, " - ", region_i),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = ifelse(show_sample_names, 6, 0),
        color = "black"
      ),
      axis.text.y = element_text(
        size = ifelse(show_sample_names, 6, 0),
        color = "black"
      ),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  if (!show_sample_names) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
      )
  }
  
  return(p)
}


all_combined_plots <- list()

for (prefix_i in c("BAC", "ARQ", "ITS")) {
  
  cat("\nGenerating Bray-Curtis heatmap for:", prefix_i, "\n")
  
  bray_amz <- make_bray_matrix_fixed_order(prefix_i, "Amazonas")
  bray_cus <- make_bray_matrix_fixed_order(prefix_i, "Cusco")
  

  write.csv(
    bray_amz,
    file.path(outdir, paste0(prefix_i, "_Amazonas_BrayCurtis_distance_matrix_fixed_order.csv"))
  )
  
  write.csv(
    bray_cus,
    file.path(outdir, paste0(prefix_i, "_Cusco_BrayCurtis_distance_matrix_fixed_order.csv"))
  )
  

  p_amz_label <- plot_bray_heatmap(
    bray_amz,
    prefix_i = prefix_i,
    region_i = "Amazonas",
    show_sample_names = TRUE
  )
  
  p_cus_label <- plot_bray_heatmap(
    bray_cus,
    prefix_i = prefix_i,
    region_i = "Cusco",
    show_sample_names = TRUE
  )
  
  p_combined_label <- p_amz_label + p_cus_label +
    plot_layout(widths = c(1, 1.2), guides = "collect") &
    theme(legend.position = "right")
  
  save_plot_all_formats(
    p_combined_label,
    file.path(outdir, paste0(prefix_i, "_beta_diversity_GENUS_BrayCurtis_HEATMAP_FIXED_ORDER_LABELLED")),
    width = 13,
    height = 6.2
  )
  

  p_amz_clean <- plot_bray_heatmap(
    bray_amz,
    prefix_i = prefix_i,
    region_i = "Amazonas",
    show_sample_names = FALSE
  )
  
  p_cus_clean <- plot_bray_heatmap(
    bray_cus,
    prefix_i = prefix_i,
    region_i = "Cusco",
    show_sample_names = FALSE
  )
  
  p_combined_clean <- p_amz_clean + p_cus_clean +
    plot_layout(widths = c(1, 1.2), guides = "collect") &
    theme(legend.position = "right")
  
  save_plot_all_formats(
    p_combined_clean,
    file.path(outdir, paste0(prefix_i, "_beta_diversity_GENUS_BrayCurtis_HEATMAP_FIXED_ORDER_CLEAN")),
    width = 12,
    height = 5.8
  )
  
  all_combined_plots[[prefix_i]] <- p_combined_clean
}


p_all_heatmaps <- all_combined_plots$BAC /
  all_combined_plots$ARQ /
  all_combined_plots$ITS +
  plot_annotation(
    title = "Genus-level Bray-Curtis dissimilarity heatmaps",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

save_plot_all_formats(
  p_all_heatmaps,
  file.path(outdir, "ALL_MARKERS_beta_diversity_GENUS_BrayCurtis_HEATMAP_FIXED_ORDER_CLEAN"),
  width = 12,
  height = 16
)


caption_text <- c(
  "Figure Sx. Genus-level Bray-Curtis dissimilarity heatmaps of soil microbial communities. Heatmaps show pairwise Bray-Curtis dissimilarities among samples calculated from relative abundance matrices of identified genera. Blue indicates lower dissimilarity, whereas red indicates higher dissimilarity. Samples are shown in the same fixed order used for taxonomic composition figures, with black lines separating post-fire chronosequence groups."
)

writeLines(
  caption_text,
  file.path(outdir, "CAPTION_BrayCurtis_heatmaps_GENUS.txt")
)


stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

zip_file <- file.path(
  beta_dir,
  paste0("BRAY_CURTIS_HEATMAPS_GENUS_FIXED_ORDER_", stamp, ".zip")
)

old_wd <- getwd()
setwd(beta_dir)

zip::zipr(
  zipfile = zip_file,
  files = file.path(
    "04_SUPPLEMENTARY",
    "BRAY_CURTIS_HEATMAPS_FIXED_SAMPLE_ORDER"
  )
)

setwd(old_wd)


cat("\n=====================================================\n")
cat("DONE: BRAY-CURTIS HEATMAPS GENERATED\n")
cat("=====================================================\n")

cat("\nOutput folder:\n")
cat(outdir, "\n")

cat("\nZIP file:\n")
cat(zip_file, "\n")

cat("\nMain heatmap files generated:\n")
cat("BAC_beta_diversity_GENUS_BrayCurtis_HEATMAP_FIXED_ORDER_CLEAN.tiff\n")
cat("ARQ_beta_diversity_GENUS_BrayCurtis_HEATMAP_FIXED_ORDER_CLEAN.tiff\n")
cat("ITS_beta_diversity_GENUS_BrayCurtis_HEATMAP_FIXED_ORDER_CLEAN.tiff\n")
cat("ALL_MARKERS_beta_diversity_GENUS_BrayCurtis_HEATMAP_FIXED_ORDER_CLEAN.tiff\n")

cat("\nLabelled versions were also generated.\n")
cat("Use labelled versions for checking sample order.\n")
cat("Use CLEAN versions for the article or supplementary figures.\n")

cat("\nCaption saved as:\n")
cat(file.path(outdir, "CAPTION_BrayCurtis_heatmaps_GENUS.txt"), "\n")

cat("=====================================================\n")
