#!/usr/bin/env Rscript
# Export a non-rarefied ASV table and its representative 16S sequences for PICRUSt2.
# Usage: Rscript 01_prepare_16s_picrust2.R bacteria|archaea PROJECT_ROOT
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[1] %in% c('bacteria','archaea')) {
  stop('Usage: Rscript 01_prepare_16s_picrust2.R bacteria|archaea <project_root>')
}
marker <- args[1]
root <- normalizePath(path.expand(args[2]), mustWork = TRUE)
if (!requireNamespace('phyloseq', quietly=TRUE) || !requireNamespace('Biostrings', quietly=TRUE)) {
  stop('Packages phyloseq and Biostrings are required')
}
library(phyloseq)
library(Biostrings)
folder <- if (marker == 'bacteria') 'RESULT BACTERIA' else 'RESULT ARQUEA'
ps_filename <- if (marker == 'bacteria') 'ps_sin_prevalencia.rds' else 'ps_sin_prevalencia_ARQ.rds'
base <- file.path(root, folder)
ps_path <- file.path(base, 'data/processed_data/phyloseq', ps_filename)
seqtab_path <- file.path(base,'data/processed_data/nochim/seqtab_nochim.rds')
if (!all(file.exists(c(ps_path,seqtab_path)))) stop('Missing phyloseq or DADA2 ASV table')
ps <- readRDS(ps_path)
seqtab <- readRDS(seqtab_path)
if (!all(taxa_names(ps) %in% colnames(seqtab))) {
  stop('The phyloseq taxon identifiers differ from the DADA2 ASV sequence table')
}
otu <- as(otu_table(ps), 'matrix')
if (!taxa_are_rows(ps)) otu <- t(otu)
stopifnot(identical(rownames(otu), taxa_names(ps)))
keep <- rowSums(otu) > 0
otu <- otu[keep,,drop=FALSE]
ids <- paste0('ASV',seq_len(nrow(otu)))
sequences <- DNAStringSet(rownames(otu))
names(sequences) <- ids
rownames(otu) <- ids
outdir <- file.path(base, 'results/picrust2_ready')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
writeXStringSet(sequences, file.path(outdir,'asv_sequences_picrust2.fna'))
write.table(data.frame(ASV=rownames(otu), otu, check.names=FALSE),
            file.path(outdir,'asv_abundance_picrust2.tsv'),
            sep='\t', row.names=FALSE, quote=FALSE)
meta <- as.data.frame(sample_data(ps), stringsAsFactors=FALSE)
meta$SampleID <- rownames(meta)
if (!setequal(colnames(otu), meta$SampleID)) stop('Sample mismatch between metadata and ASV counts')
meta <- meta[match(colnames(otu), meta$SampleID),,drop=FALSE]
write.csv(meta,file.path(outdir,'metadata_picrust2.csv'),row.names=FALSE)
write.csv(data.frame(ASV=ids,Sequence=as.character(sequences)),
          file.path(outdir,'asv_sequence_map.csv'),row.names=FALSE)
cat('Prepared',marker,':',nrow(otu),'ASVs and',ncol(otu),'samples\n')
