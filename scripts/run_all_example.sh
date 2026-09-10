#!/usr/bin/env bash
set -euo pipefail

# Do not commit real data paths or API keys.
# Set this in your shell before running:
#   export ALPHA_GENOME_API_KEY="..."

INPUT_TSV="/absolute/path/to/all_482_for_avi.tsv"
OUTDIR="/absolute/path/to/alphagenome_avi_outputs"

python scripts/score_avi_all482.py \
  --input-tsv "$INPUT_TSV" \
  --outdir "$OUTDIR" \
  --api-key-env ALPHA_GENOME_API_KEY \
  --scorer AVI_SCORE \
  --resume

python scripts/analyze_avi_results.py \
  --scored-tsv "$OUTDIR/01_all482_with_avi.tsv" \
  --outdir "$OUTDIR"
