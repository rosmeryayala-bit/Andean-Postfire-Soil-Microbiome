#!/usr/bin/env Rscript
# Usage: Rscript scripts/02_sequence_processing/00_dada2_all_markers.R \
#          bacteria /PATH/TO/ORIGINAL_PROJECT /PATH/TO/SILVA_TRAINSET /PATH/TO/Metadata.csv
# Repeat for archaea and fungi with the appropriate taxonomy training file.
# Taxonomy is assigned using the supplied marker-specific training FASTA.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: 00_dada2_all_markers.R <bacteria|archaea|fungi> <project_root> <taxonomy_train_fasta> <metadata.csv>")
}
marker <- args[[1L]]
repo_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(repo_file) != 1L) stop("Invoke using Rscript <script_path> ...")
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", repo_file)), "..", ".."), mustWork=TRUE)
source(file.path(repo_root, "config", "marker_parameters.R"))
if (!marker %in% names(marker_configs)) stop("Unknown marker: ",marker)
cfg <- marker_configs[[marker]]
project_root <- normalizePath(path.expand(args[[2L]]), mustWork=TRUE)
training_file <- normalizePath(path.expand(args[[3L]]), mustWork=TRUE)
metadata_file <- normalizePath(path.expand(args[[4L]]), mustWork=TRUE)
base <- file.path(project_root,cfg$folder)
if (!dir.exists(base)) stop("Missing marker directory: ",base)
if (!requireNamespace("dada2", quietly=TRUE)) stop("Install Bioconductor dada2 first")
if (!requireNamespace("phyloseq", quietly=TRUE)) stop("Install Bioconductor phyloseq first")
if (!requireNamespace("Biostrings", quietly=TRUE)) stop("Install Bioconductor Biostrings first")
library(dada2); library(phyloseq); library(Biostrings)
processed <- file.path(base,"data/processed_data")
outputs <- file.path(base,"results")
for (d in c("filtered_F","filtered_R","error_models","asv_objects","merge","seqtab","nochim","taxonomy","phyloseq")) {
  dir.create(file.path(processed,d),recursive=TRUE,showWarnings=FALSE)
}
for (d in c("tracking","tables","figures","asv_sequences")) {
  dir.create(file.path(outputs,d),recursive=TRUE,showWarnings=FALSE)
}
# Verify both filenames and sample pairs: sorting separately without checking IDs can mispair reads.
raw <- file.path(base,"data/rawdata")
forward <- sort(list.files(raw, pattern=cfg$read_pattern, full.names=TRUE))
reverse <- sort(list.files(raw, pattern=cfg$reverse_pattern, full.names=TRUE))
if (!length(forward) || length(forward) != length(reverse)) stop("Missing/unequal FASTQ pairs for ",marker)
remove_suffix <- function(paths,suffix) sub(suffix,"",basename(paths))
fids <- remove_suffix(forward,cfg$read_pattern)
rids <- remove_suffix(reverse,cfg$reverse_pattern)
if (anyDuplicated(fids) || anyDuplicated(rids) || !setequal(fids,rids)) stop("FASTQ sample identifiers do not pair")
reverse <- reverse[match(fids,rids)]
stopifnot(identical(fids,remove_suffix(reverse,cfg$reverse_pattern)))
cat("Marker:",marker,"; matched FASTQ pairs:",length(forward),"\n")
filtFs <- file.path(processed,"filtered_F",paste0(fids,"_filt_1.fastq.gz"))
filtRs <- file.path(processed,"filtered_R",paste0(fids,"_filt_2.fastq.gz"))
params <- list(fwd=forward, filt=filtFs, rev=reverse, filt.rev=filtRs,
  maxN=0, maxEE=cfg$maxEE, truncQ=2, rm.phix=TRUE, compress=TRUE,
  multithread=FALSE, matchIDs=TRUE, verbose=TRUE)
if (!is.null(cfg$truncLen)) params$truncLen <- cfg$truncLen
if (!is.null(cfg$minLen)) params$minLen <- cfg$minLen
filtered <- do.call(dada2::filterAndTrim,params)
write.csv(filtered,file.path(outputs,"tracking",paste0("filtering_reads_",marker,".csv")))
saveRDS(filtered,file.path(outputs,"tracking",paste0("filtering_reads_",marker,".rds")))
# Read quality-dependent error models, no pooling in ASV inference.
errF <- learnErrors(filtFs,nbases=1e8,multithread=FALSE,randomize=TRUE,verbose=TRUE)
errR <- learnErrors(filtRs,nbases=1e8,multithread=FALSE,randomize=TRUE,verbose=TRUE)
saveRDS(errF,file.path(processed,"error_models","errF.rds"))
saveRDS(errR,file.path(processed,"error_models","errR.rds"))
pdf(file.path(outputs,"figures",paste0("error_model_",marker,".pdf")),width=12,height=6)
par(mfrow=c(1,2)); plotErrors(errF,nominalQ=TRUE); plotErrors(errR,nominalQ=TRUE); dev.off()
derepFs <- derepFastq(filtFs,verbose=TRUE)
derepRs <- derepFastq(filtRs,verbose=TRUE)
names(derepFs) <- fids; names(derepRs) <- fids
dadaFs <- dada(derepFs,err=errF,multithread=FALSE,pool=FALSE)
dadaRs <- dada(derepRs,err=errR,multithread=FALSE,pool=FALSE)
for (nm in c("derepFs","derepRs","dadaFs","dadaRs")) saveRDS(get(nm),file.path(processed,"asv_objects",paste0(nm,".rds")))
mergers <- mergePairs(dadaFs,derepFs,dadaRs,derepRs,verbose=TRUE)
saveRDS(mergers,file.path(processed,"merge","mergers.rds"))
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab,file.path(processed,"seqtab","seqtab.rds"))
seqtab.nochim <- removeBimeraDenovo(seqtab,method="consensus",multithread=FALSE,verbose=TRUE)
saveRDS(seqtab.nochim,file.path(processed,"nochim","seqtab_nochim.rds"))
write.csv(data.frame(Sample=rownames(seqtab.nochim),Reads=rowSums(seqtab.nochim)),file.path(outputs,"tables",paste0("read_totals_",marker,".csv")),row.names=FALSE)
# Assign taxonomy; the exact reference dataset must be supplied for replication.
taxa <- assignTaxonomy(seqtab.nochim,training_file,multithread=FALSE,tryRC=TRUE,verbose=TRUE)
saveRDS(taxa,file.path(processed,"taxonomy",paste0("taxa_",marker,".rds")))
write.csv(taxa,file.path(outputs,"tables",paste0("taxonomy_",marker,".csv")))
# Filename compatibility with original downstream commands.
legacy_taxa <- switch(marker,bacteria="taxa_silva.rds",archaea="taxa_silva_ARQ.rds",fungi="taxa_unite.rds")
saveRDS(taxa,file.path(processed,"taxonomy",legacy_taxa))
# Build ASV FASTA with matching IDs; save map to the original sequence-header column.
asv_seqs <- colnames(seqtab.nochim)
asv_ids <- paste0("ASV",seq_along(asv_seqs))
seqs <- DNAStringSet(asv_seqs); names(seqs) <- asv_ids
fasta_name <- switch(marker,bacteria="ASVs.fa",archaea="ASVs_ARQ.fa",fungi="ASVs_ITS.fa")
writeXStringSet(seqs,file.path(outputs,"asv_sequences",fasta_name))
write.csv(data.frame(ASV_ID=asv_ids,Sequence=asv_seqs),file.path(outputs,"tables",paste0("ASV_ID_mapping_",marker,".csv")),row.names=FALSE)
# Reconcile metadata. Raw read labels may contain _BAC.raw, _ARQ, _FUN and A_/C_ prefix.
md <- read.csv(metadata_file,sep=";",header=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
if (!ncol(md) || !nrow(md)) stop("Metadata file empty/invalid")
clean_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.raw$","",x)
  x <- sub("_(BAC|ARQ|FUN)$","",x)
  x <- sub("^A_","A",x)
  x <- sub("^C_","C",x)
  x
}
md_id <- clean_id(md[[1]])
seq_id <- clean_id(rownames(seqtab.nochim))
if (anyDuplicated(md_id) || anyDuplicated(seq_id)) stop("Duplicate sample ID after normalization")
if (!all(seq_id %in% md_id)) stop("Metadata missing these samples: ",paste(setdiff(seq_id,md_id),collapse=", "))
md <- md[match(seq_id,md_id),,drop=FALSE]
rownames(md) <- rownames(seqtab.nochim)
stopifnot(nrow(md) == nrow(seqtab.nochim))
if (!setequal(colnames(seqtab.nochim),rownames(taxa))) stop("Taxonomy and ASV table disagree")
taxa <- taxa[colnames(seqtab.nochim),,drop=FALSE]
ps <- phyloseq(otu_table(seqtab.nochim,taxa_are_rows=FALSE),tax_table(taxa),sample_data(md))
phylo_out <- file.path(processed,"phyloseq")
saveRDS(ps,file.path(phylo_out,cfg$raw_phyloseq))
if (marker == "bacteria") {
  ranks <- rank_names(ps); keep <- rep(TRUE,ntaxa(ps))
  if ("Order" %in% ranks) keep <- keep & (is.na(tax_table(ps)[,"Order"]) | tax_table(ps)[,"Order"] != "Chloroplast")
  if ("Family" %in% ranks) keep <- keep & (is.na(tax_table(ps)[,"Family"]) | tax_table(ps)[,"Family"] != "Mitochondria")
  curated <- prune_taxa(keep,ps)
} else if (marker == "archaea") {
  if (!"Kingdom" %in% rank_names(ps)) stop("No Kingdom column for Archaea-domain filtering")
  curated <- subset_taxa(ps,Kingdom == "Archaea")
} else curated <- ps
curated <- prune_taxa(taxa_sums(curated) > 0,curated)
saveRDS(curated,file.path(phylo_out,cfg$curated_phyloseq))
if (marker == "bacteria") saveRDS(curated,file.path(phylo_out,"ps_sin_prevalencia.rds"))
if (marker == "archaea") saveRDS(curated,file.path(phylo_out,"ps_only_archaea_ARQ.rds"))
getN <- function(x) sum(getUniques(x))
counts <- data.frame(Sample=fids, Filtered=filtered[,2],
  DenoisedF=sapply(dadaFs,getN),DenoisedR=sapply(dadaRs,getN),
  Merged=sapply(mergers,function(x) if (nrow(x)) sum(x$abundance) else 0),
  Nonchim=rowSums(seqtab.nochim)[fids],row.names=NULL)
write.csv(counts,file.path(outputs,"tracking",paste0("read_tracking_",marker,".csv")),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(outputs,"tracking",paste0("sessionInfo_",marker,".txt")))
cat("DADA2 completed for",marker,"; samples",nsamples(curated),"; ASVs",ntaxa(curated),"\n")
