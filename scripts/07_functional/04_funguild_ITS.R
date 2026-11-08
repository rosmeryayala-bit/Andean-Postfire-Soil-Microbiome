#!/usr/bin/env Rscript
# Fungal ITS2 functional guild annotation (FUNGuild through microeco).
# Usage: Rscript 04_funguild_ITS.R <project_root>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 04_funguild_ITS.R <project_root>")
root <- normalizePath(path.expand(args[[1L]]), mustWork = TRUE)
for (pkg in c("phyloseq", "microeco")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install the required R package: ", pkg)
}
library(phyloseq)
library(microeco)
base <- file.path(root, "RESULT FUNGAL")
ps_file <- file.path(base, "data", "processed_data", "phyloseq", "ps_raw_ITS.rds")
if (!file.exists(ps_file)) stop("Missing ITS2 phyloseq object: ", ps_file)
ps <- readRDS(ps_file)
# Preserve observed ASV abundances; no additional rarefaction in functional analysis.
otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)
tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
meta <- as.data.frame(sample_data(ps), stringsAsFactors = FALSE)
if (!setequal(rownames(otu), rownames(tax))) stop("ASV identifiers do not match taxonomy")
if (!setequal(colnames(otu), rownames(meta))) stop("Sample identifiers do not match metadata")
tax <- tax[rownames(otu), , drop = FALSE]
meta <- meta[colnames(otu), , drop = FALSE]
mt <- microtable$new(otu_table = as.data.frame(otu, check.names = FALSE),
                     sample_table = meta, tax_table = tax)
mt$tidy_dataset()
tf <- trans_func$new(mt)
tf$for_what <- "fungi"
tf$cal_func(fungi_database = "FUNGuild",
            FUNGuild_confidence = c("Highly Probable", "Probable"))
tf$cal_func_FR(abundance_weighted = TRUE)
if (is.null(tf$res_func_FR)) stop("No FUNGuild abundance-weighted profiles returned")
outdir <- file.path(base, "results", "functional_funguild")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(tf$res_func_FR, file.path(outdir, "ITS_FUNGuild_abundance_weighted_percent.csv"))
if (!is.null(tf$res_func)) {
  write.csv(tf$res_func, file.path(outdir, "ITS_FUNGuild_ASV_annotations.csv"))
}
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_FUNGuild.txt"))
cat("FUNGuild annotation completed for ITS2\n")
