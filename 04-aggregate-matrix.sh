#!/bin/bash
#SBATCH -J AGG_WGBS
#SBATCH --output logs/%x.%A_%a.%N.%j.txt
#SBATCH -p exbio-cpu
#SBATCH -c 20
#SBATCH --mem-per-cpu 16G
#SBATCH -t 0-1:0:00 # days-hh:mm:ss

Rscript 04-aggregate_matrix.R