#!/bin/bash
#
#SBATCH -J ML-local
#SBATCH --output %x.%j.%N.txt
#SBATCH -p exbio-cpu
#SBATCH -c 40
#SBATCH --mem 50G
#SBATCH -t 2-00:00:00 # days-hh:mm:ss

srun Rscript 06-ml-local.R
