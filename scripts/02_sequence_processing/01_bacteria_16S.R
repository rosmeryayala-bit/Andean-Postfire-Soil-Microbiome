#!/usr/bin/env Rscript
# DADA2: BACTERIA
# Usage: Rscript 01_bacteria_16S.R /data/project /references/training_set.fa.gz /data/Metadata.csv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop('Usage: Rscript 01_bacteria_16S.R <project_root> <taxonomy_training_fasta> <metadata.csv>')
self <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(self) != 1L) stop('Run using Rscript')
engine <- file.path(dirname(normalizePath(sub('^--file=', '', self))), '00_dada2_all_markers.R')
status <- system2(file.path(R.home('bin'), 'Rscript'), shQuote(c(engine,'bacteria',args)), stdout='', stderr='')
if (!identical(as.integer(status), 0L)) stop('DADA2 process failed with exit status ', status)
