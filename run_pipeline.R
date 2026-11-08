#!/usr/bin/env Rscript
# Main entry point for sequence-processing workflows.
# Downstream statistical analyses are organized as separate R scripts.
# Example:
# Rscript run_pipeline.R --list
# Rscript run_pipeline.R --preflight /path/to/project
# Rscript run_pipeline.R --preprocess bacteria /path/to/project /path/to/train.fa.gz /path/to/Metadata.csv
args <- commandArgs(trailingOnly=TRUE)
scriptarg <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
if (length(scriptarg) != 1L) stop("Use Rscript run_pipeline.R ...")
repo <- dirname(normalizePath(sub("^--file=","",scriptarg)))
if (!length(args) || args[[1L]] == "--list") {
  tab <- read.csv(file.path(repo,"PIPELINE_ORDER.csv"),stringsAsFactors=FALSE)
  print(tab[,c("step","domain","analysis","script")],row.names=FALSE,right=FALSE)

  quit(status=0L)
}
if (args[[1L]] == "--preflight") {
  if (length(args)!=2L) stop("Provide the original project root")
  root <- normalizePath(args[[2L]],mustWork=TRUE)
  for (folder in c("RESULT BACTERIA","RESULT ARQUEA","RESULT FUNGAL")) {
    ff <- file.path(root,folder,"data/rawdata")
    cat(if (dir.exists(ff)) "[OK] " else "[MISSING] ",ff,"\n",sep="")
  }
  for (f in c("config/marker_parameters.R","scripts/01_quality_control/01_fastqc_multiqc.sh","scripts/02_sequence_processing/00_dada2_all_markers.R")) {
    cat(if (file.exists(file.path(repo,f))) "[OK] " else "[MISSING] ",f,"\n",sep="")
  }
  quit(status=0L)
}
if (args[[1L]] == "--preprocess") {
  if (length(args) != 5L) stop("Usage: --preprocess bacteria|archaea|fungi <project_root> <training_fasta> <metadata.csv>")
  marker <- args[[2L]]
  if (!marker %in% c("bacteria","archaea","fungi")) stop("Invalid marker")
  paths <- normalizePath(args[3:5],mustWork=TRUE)
  script <- file.path(repo,"scripts/02_sequence_processing/00_dada2_all_markers.R")
  cat("DADA2 preprocessing only for marker:",marker,"\n")
  cat("FastQC/MultiQC is a separate upstream step (use bash script).\n")
  rscript <- file.path(R.home("bin"),"Rscript")
  status <- system2(rscript, args=shQuote(c(script,marker,paths)), stdout="", stderr="")
  if (status != 0L) stop("Preprocessing returned non-zero exit status: ",status)
  quit(status=0L)
}
stop("Unknown option. Use --list, --preflight, or --preprocess")
