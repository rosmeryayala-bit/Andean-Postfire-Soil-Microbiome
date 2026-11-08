# Data inputs

This repository contains code for FASTQ processing, microbial community analyses, soil chemical analyses and functional inference. Input datasets, taxonomic reference training FASTA files and intermediate R objects are stored separately.

## Raw sequencing inputs

| Marker | Paired-end FASTQ directory | Expected suffixes |
|---|---|---|
| Bacteria | `RESULT BACTERIA/data/rawdata/` | `_1.fastq.gz`, `_2.fastq.gz` |
| Archaea | `RESULT ARQUEA/data/rawdata/` | `_ARQ_1.fastq.gz`, `_ARQ_2.fastq.gz` |
| Fungi | `RESULT FUNGAL/data/rawdata/` | `_1.fastq.gz`, `_2.fastq.gz` |

The metadata file `Metadata.csv` is semicolon-delimited and contains sample identifiers, region, treatment and years since fire. Each marker has an associated training reference sequence file: SILVA NR99 v138.2 for bacterial and archaeal 16S and UNITE Fungi v10.0 for fungal ITS2.

## Downstream files

The biological community analyses use three processed phyloseq objects:

- `RESULT BACTERIA/data/processed_data/phyloseq/ps_sin_prevalencia.rds`
- `RESULT ARQUEA/data/processed_data/phyloseq/ps_sin_prevalencia_ARQ.rds`
- `RESULT FUNGAL/data/processed_data/phyloseq/ps_raw_ITS.rds`

Soil chemical analyses and genus–soil associations use an 18-unit site–treatment table, typically stored at `RESULTADOS_SOIL_UNIQUE_UNITS_FINAL/TABLES/soil_metadata_unique_units.csv`. The alpha diversity module reads values and model assumption summaries from the `RESULTADOS_ALPHA_DIVERSITY_GENUS_SUPUESTOS` directory. Relevant output filenames and output directories are defined in each analysis script.

Sequence reads can be deposited with an appropriate sequence repository, with corresponding accession numbers cited separately. The repository does not distribute raw FASTQ or large intermediate files.
