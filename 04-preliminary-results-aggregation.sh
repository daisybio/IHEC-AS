#!/bin/bash
#SBATCH -J AGG
#SBATCH --output %x.%j.%N.txt
#SBATCH -p exbio-cpu
#SBATCH -c 20
#SBATCH --mem 36G
#SBATCH -t 1-00:00:00 # days-hh:mm:ss


srun Rscript 04-preliminary-results-aggregation.R
