#!/bin/bash
#
#SBATCH -J ML-local
#SBATCH --output %x.%j.%N.txt
#SBATCH -p exbio-cpu
#SBATCH -c 10
#SBATCH --mem-per-cpu 15G
#SBATCH -t 4-00:00:00 # days-hh:mm:ss

Rscript 06-ml-local.R
