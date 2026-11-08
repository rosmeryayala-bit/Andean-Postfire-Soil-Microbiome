#!/usr/bin/env bash
# PICRUSt2 v2.0.3b0 workflow for bacterial or archaeal 16S data.
# Usage: bash 02_run_16s_picrust2.sh bacteria|archaea /path/to/project_root [threads]
set -euo pipefail
MARKER=${1:?bacteria or archaea}; ROOT=${2:?project root}; THREADS=${3:-12}
case "$MARKER" in
  bacteria) FOLDER='RESULT BACTERIA' ;;
  archaea) FOLDER='RESULT ARQUEA' ;;
  *) echo 'Invalid marker' >&2; exit 2 ;;
esac
command -v picrust2_pipeline.py >/dev/null || { echo 'PICRUSt2 executable not found' >&2; exit 1; }
BASE="$ROOT/$FOLDER"
IN="$BASE/results/picrust2_ready"
OUT="$BASE/results/picrust2_out"
mkdir -p "$BASE/results/logs"
test -f "$IN/asv_sequences_picrust2.fna"
test -f "$IN/asv_abundance_picrust2.tsv"
picrust2_pipeline.py \
  -s "$IN/asv_sequences_picrust2.fna" \
  -i "$IN/asv_abundance_picrust2.tsv" \
  -o "$OUT" -p "$THREADS" --stratified --per_sequence_contrib \
  2>&1 | tee "$BASE/results/logs/picrust2_${MARKER}.log"
