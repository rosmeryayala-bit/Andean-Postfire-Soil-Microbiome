# Marker-specific sequence filtering and taxonomy settings
marker_configs <- list(
  bacteria = list(
    folder = 'RESULT BACTERIA',
    read_pattern = '_1.fastq.gz$', reverse_pattern = '_2.fastq.gz$',
    truncLen = c(250L,250L), maxEE = c(3,4), minLen = NULL,
    database = 'SILVA NR99 v138.2',
    raw_phyloseq = 'ps_raw.rds', curated_phyloseq = 'ps_noplast.rds'
  ),
  archaea = list(
    folder = 'RESULT ARQUEA',
    read_pattern = '_ARQ_1.fastq.gz$', reverse_pattern = '_ARQ_2.fastq.gz$',
    truncLen = c(250L,250L), maxEE = c(3,4), minLen = NULL,
    database = 'SILVA NR99 v138.2',
    raw_phyloseq = 'ps_raw_ARQ.rds', curated_phyloseq = 'ps_sin_prevalencia_ARQ.rds'
  ),
  fungi = list(
    folder = 'RESULT FUNGAL',
    read_pattern = '_1.fastq.gz$', reverse_pattern = '_2.fastq.gz$',
    truncLen = NULL, maxEE = c(2,2), minLen = 50L,
    database = 'UNITE Fungi v10.0',
    raw_phyloseq = 'ps_raw_ITS.rds', curated_phyloseq = 'ps_raw_ITS.rds'
  )
)

marker_parameters <- marker_configs
