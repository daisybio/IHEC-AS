#!/bin/bash
# Pangolin splice-site scores (dilated residual CNN, tissue-specific P(site used)).
# Inputs:  processed_data/pangolin_events.csv  (written by 03)
#          data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# Output:  processed_data/pangolin_scores.csv
# GPU optional; CPU ~40 min–2 hr at default batch_size=256.
set -euo pipefail

module load miniforge3/24.7.1
export PYTHONUNBUFFERED=1
mamba run --no-capture-output -n ihec-as python scripts/compute_pangolin_scores.py "$@"
