#!/bin/bash
# Pangolin splice-site scores (dilated residual CNN, tissue-specific P(site used)).
# Per transcript_filter (via TRANSCRIPT_FILTER env, PLAN §4.12e):
# Inputs:  processed_data/pangolin_events_${TRANSCRIPT_FILTER}.csv  (written by 03)
#          data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# Output:  processed_data/pangolin_scores_${TRANSCRIPT_FILTER}.csv
# GPU optional; CPU ~40 min–2 hr at default batch_size=256.
set -euo pipefail

: "${TRANSCRIPT_FILTER:?TRANSCRIPT_FILTER env var must be set}"
export TRANSCRIPT_FILTER
module load miniforge3/24.7.1
export PYTHONUNBUFFERED=1
# Env selectable so the rule can run the titan-compatible torch (cu118) env on
# the idle gpu01 titans instead of the a40-only cu128 `ihec-as`. Defaults to
# ihec-as for standalone/interactive use.
PANGOLIN_ENV="${PANGOLIN_ENV:-ihec-as}"
mamba run --no-capture-output -n "${PANGOLIN_ENV}" python scripts/compute_pangolin_scores.py "$@"
