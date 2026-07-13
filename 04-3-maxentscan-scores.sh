#!/bin/bash
# MaxEntScan 5' and 3' splice-site scores (per transcript_filter).
# Inputs:  processed_data/{5ss,5ss_up,3ss,3ss_down}_${TRANSCRIPT_FILTER}.fasta  (03)
# Outputs: processed_data/{5scores,5up_scores,3scores,3down_scores}_${TRANSCRIPT_FILTER}.txt
set -euo pipefail

tf="${TRANSCRIPT_FILTER:?TRANSCRIPT_FILTER env var must be set}"

maxentscan_score5.pl processed_data/5ss_${tf}.fasta      > processed_data/5scores_${tf}.txt
maxentscan_score5.pl processed_data/5ss_up_${tf}.fasta   > processed_data/5up_scores_${tf}.txt
maxentscan_score3.pl processed_data/3ss_${tf}.fasta      > processed_data/3scores_${tf}.txt
maxentscan_score3.pl processed_data/3ss_down_${tf}.fasta > processed_data/3down_scores_${tf}.txt
