#!/usr/bin/env bash
# FastQC 0.12.1 + MultiQC 1.32; one run per marker.
# Usage: bash scripts/preprocessing/00_fastqc_multiqc.sh /path/to/project_root bacteria [threads]
set -euo pipefail
ROOT=${1:?Provide original project root}; MARKER=${2:?bacteria|archaea|fungi}; THREADS=${3:-8}
case "$MARKER" in
  bacteria) FOLDER='RESULT BACTERIA'; PATTERN='*_1.fastq.gz'; REV='*_2.fastq.gz' ;;
  archaea) FOLDER='RESULT ARQUEA'; PATTERN='*_ARQ_1.fastq.gz'; REV='*_ARQ_2.fastq.gz' ;;
  fungi) FOLDER='RESULT FUNGAL'; PATTERN='*_1.fastq.gz'; REV='*_2.fastq.gz' ;;
  *) echo 'Marker must be bacteria|archaea|fungi' >&2; exit 2 ;;
esac
command -v fastqc >/dev/null || { echo 'FastQC not installed' >&2; exit 1; }
command -v multiqc >/dev/null || { echo 'MultiQC not installed' >&2; exit 1; }
BASE="$ROOT/$FOLDER"; RAW="$BASE/data/rawdata"; mkdir -p "$BASE/results/qc_forward" "$BASE/results/qc_reverse" "$BASE/results/multiqc_forward" "$BASE/results/multiqc_reverse" "$BASE/results/multiqc_total"
shopt -s nullglob
FW=("$RAW"/$PATTERN); RV=("$RAW"/$REV)
if (( ${#FW[@]}==0 || ${#FW[@]}!=${#RV[@]} )); then echo 'Missing or unequal FASTQ pairs' >&2; exit 1; fi
printf 'Marker %s: %s forward / %s reverse\n' "$MARKER" "${#FW[@]}" "${#RV[@]}"
fastqc "${FW[@]}" -o "$BASE/results/qc_forward" -t "$THREADS"
fastqc "${RV[@]}" -o "$BASE/results/qc_reverse" -t "$THREADS"
multiqc "$BASE/results/qc_forward" -o "$BASE/results/multiqc_forward" --force
multiqc "$BASE/results/qc_reverse" -o "$BASE/results/multiqc_reverse" --force
multiqc "$BASE/results/qc_forward" "$BASE/results/qc_reverse" -o "$BASE/results/multiqc_total" --force
