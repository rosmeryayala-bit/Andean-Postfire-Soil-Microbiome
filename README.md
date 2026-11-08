# Soil microbiome responses to wildfires in the Peruvian Andes

Bioinformatics and statistical analysis scripts for high-Andean grassland soil microbiomes from Amazonas and Cusco, Peru. The study includes bacterial 16S rRNA V3–V4, archaeal 16S rRNA V4–V5, and fungal ITS2 amplicon datasets.

The sampling design comprises **54 soil samples for molecular analyses** (27 wildfire-affected and 27 unburned reference samples, with three biological replicates per site–treatment combination) and **18 site–treatment soil units** for chemical characterisation. Years since fire are 2, 4, 6 and 10 in Amazonas and 1, 5, 10, 15 and 20 in Cusco.

## Workflow

| Stage | Scripts | Description |
|---|---|---|
| 1 | `scripts/01_quality_control/` | Forward and reverse FASTQ quality checks with FastQC and MultiQC |
| 2 | `scripts/02_sequence_processing/` | Marker-specific DADA2 filtering, error learning, dereplication, ASV inference, paired-read merging and chimera removal |
| 3 | `scripts/02_sequence_processing/` | Taxonomic classification with SILVA or UNITE and phyloseq object generation |
| 4 | `scripts/03_taxonomy/` | Taxonomic profiles and shared/exclusive genera |
| 5 | `scripts/04_diversity/` | Rarefaction/extrapolation, genus-level Hill diversity, Bray–Curtis, PCoA, PERMANOVA and dispersion tests |
| 6 | `scripts/05_soil/` | Soil chemical properties and edaphic statistical analyses |
| 7 | `scripts/06_associations/` | Region-specific genus–soil associations using MaAsLin2 |
| 8 | `scripts/07_functional/` | 16S PICRUSt2 gene/pathway prediction, FUNGuild and FungalTraits analyses of ITS2 |

The ordered script inventory is available in [`PIPELINE_ORDER.csv`](PIPELINE_ORDER.csv).

## Sequence processing parameters

| Marker | Truncation length | Maximum expected errors | Taxonomic reference |
|---|---|---|---|
| Bacterial 16S | 250 bp forward; 250 bp reverse | `c(3, 4)` | SILVA NR99 v138.2 |
| Archaeal 16S | 250 bp forward; 250 bp reverse | `c(3, 4)` | SILVA NR99 v138.2 |
| Fungal ITS2 | No fixed-length truncation; minimum read length 50 bp | `c(2, 2)` | UNITE Fungi v10.0 |

The DADA2 parameters additionally include `maxN = 0`, `truncQ = 2`, removal of PhiX sequences, separate forward/reverse error models with up to 100 million training bases, ASV inference without pooling (`pool = FALSE`), and consensus chimera removal. Full marker settings are defined in [`config/marker_parameters.R`](config/marker_parameters.R).

## Quick start

The code uses R and Bash. From the root of this repository, with the sample data and required reference training FASTA files available:

```bash
# Show the processing stages
Rscript run_pipeline.R --list

# Examine input folder availability
Rscript run_pipeline.R --preflight /path/to/project

# Paired-end read quality control
bash scripts/01_quality_control/01_fastqc_multiqc.sh /path/to/project bacteria
bash scripts/01_quality_control/01_fastqc_multiqc.sh /path/to/project archaea
bash scripts/01_quality_control/01_fastqc_multiqc.sh /path/to/project fungi

# Sequence processing and taxonomy
Rscript run_pipeline.R --preprocess bacteria /path/to/project /path/to/SILVA_training.fa.gz /path/to/Metadata.csv
Rscript run_pipeline.R --preprocess archaea /path/to/project /path/to/SILVA_training.fa.gz /path/to/Metadata.csv
Rscript run_pipeline.R --preprocess fungi /path/to/project /path/to/UNITE_training.fa.gz /path/to/Metadata.csv
```

The three marker-specific wrappers can also be run from `scripts/02_sequence_processing/`. Scripts for diversity, edaphic data and genus–soil associations follow the sequence in `PIPELINE_ORDER.csv`. The original analysis modules use the project input and output directory conventions described in [`DATA_INPUTS.md`](DATA_INPUTS.md).

## Functional inference

### Bacterial and archaeal 16S

PICRUSt2 is used with non-rarefied ASV abundance tables and corresponding representative sequences. Predicted KO profiles are normalised to counts per million (CPM), and low-abundance KO terms with mean CPM below 1,000 are excluded. Pathway summaries are exported separately.

```bash
Rscript scripts/07_functional/01_prepare_16s_picrust2.R bacteria /path/to/project
bash scripts/07_functional/02_run_16s_picrust2.sh bacteria /path/to/project
Rscript scripts/07_functional/03_picrust2_KO_pathways.R bacteria /path/to/project

Rscript scripts/07_functional/01_prepare_16s_picrust2.R archaea /path/to/project
bash scripts/07_functional/02_run_16s_picrust2.sh archaea /path/to/project
Rscript scripts/07_functional/03_picrust2_KO_pathways.R archaea /path/to/project
```

### Fungal ITS2

Fungal ecological guilds are annotated with FUNGuild using the `microeco` workflow, retaining the `Probable` and `Highly Probable` confidence categories. FungalTraits annotations complement functional interpretation. Abundance-weighted profiles are exported for downstream comparisons.

```bash
Rscript scripts/07_functional/04_funguild_ITS.R /path/to/project
Rscript scripts/07_functional/05_fungaltraits_ITS.R /path/to/project
```

## Software and data

Methods describe R v4.4.2, DADA2 v1.34.0, phyloseq v1.50.0, vegan v2.7-3, MaAsLin2 v1.20.0, FastQC v0.12.1, MultiQC v1.32, PICRUSt2 v2.0.3b0, and the indicated reference database releases. Additional R packages are loaded by the corresponding analysis scripts. Their exact environment can be recorded with `sessionInfo()`.

FASTQ files, database FASTA files, sample metadata, processed phyloseq objects and soil chemical measurements are external inputs and are not included in this code repository. See [`DATA_INPUTS.md`](DATA_INPUTS.md) for expected paths and file formats. To reproduce the numerical results, the scripts must be run with the corresponding input datasets, training reference files and software environment.
