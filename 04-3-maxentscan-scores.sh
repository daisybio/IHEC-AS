#!/bin/bash
# MaxEntScan 5' and 3' splice-site scores.
# Inputs:  processed_data/{5ss,5ss_up,3ss,3ss_down}.fasta  (written by 03)
# Outputs: processed_data/{5scores,5up_scores,3scores,3down_scores}.txt
set -euo pipefail

maxentscan_score5.pl processed_data/5ss.fasta      > processed_data/5scores.txt
maxentscan_score5.pl processed_data/5ss_up.fasta   > processed_data/5up_scores.txt
maxentscan_score3.pl processed_data/3ss.fasta      > processed_data/3scores.txt
maxentscan_score3.pl processed_data/3ss_down.fasta > processed_data/3down_scores.txt
