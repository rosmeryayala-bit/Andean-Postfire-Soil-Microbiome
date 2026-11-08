#!/usr/bin/env Rscript
# Summarize PICRUSt2 functional predictions as counts per million (CPM).
# Usage: Rscript 03_picrust2_KO_pathways.R bacteria|archaea PROJECT_ROOT
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L || !args[1] %in% c('bacteria','archaea')) {
  stop('Usage: Rscript 03_picrust2_KO_pathways.R bacteria|archaea <project_root>')
}
marker <- args[1]; root <- normalizePath(path.expand(args[2]), mustWork=TRUE)
folder <- if (marker=='bacteria') 'RESULT BACTERIA' else 'RESULT ARQUEA'
base <- file.path(root, folder)
indir <- file.path(base,'results/picrust2_out')
outdir <- file.path(base,'results/functional_picrust2')
dir.create(outdir, recursive=TRUE,showWarnings=FALSE)
read_picrust <- function(path) {
  if (!file.exists(path)) stop('Missing PICRUSt2 output: ',path)
  tab <- read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,quote='',comment.char='')
  if (ncol(tab)<2L) stop('Invalid functional table: ',path)
  names(tab)[1] <- 'Feature'
  mat <- as.matrix(tab[,-1,drop=FALSE]); storage.mode(mat) <- 'numeric'
  rownames(mat) <- tab$Feature
  if (anyNA(mat)) stop('Non-numeric functional abundance values')
  mat
}
normalize_cpm <- function(mat) {
  depth <- colSums(mat)
  if (any(depth<=0)) stop('Samples with zero functional abundance: ',paste(colnames(mat)[depth<=0],collapse=', '))
  sweep(mat,2,depth,'/')*1e6
}
ko_file <- file.path(indir,'KO_metagenome_out/pred_metagenome_unstrat.tsv.gz')
path_file <- file.path(indir,'pathways_out/path_abun_unstrat.tsv.gz')
ko <- normalize_cpm(read_picrust(ko_file))
# The KO abundance cutoff is applied to the mean CPM across samples.
ko <- ko[rowMeans(ko)>=1000,,drop=FALSE]
write.csv(ko,file.path(outdir,'KO_CPM_mean_ge_1000.csv'))
pathways <- normalize_cpm(read_picrust(path_file))
write.csv(pathways,file.path(outdir,'Pathways_CPM.csv'))
write.csv(data.frame(Feature=rownames(pathways), MeanCPM=rowMeans(pathways)),
          file.path(outdir,'Pathways_mean_CPM.csv'),row.names=FALSE)
cat(marker,':',nrow(ko),'KO terms after filtering;',nrow(pathways),'pathways\n')
