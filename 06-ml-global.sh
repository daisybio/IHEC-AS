#!/bin/bash
#
#SBATCH -J ML-global
#SBATCH --output %x.%j.%N.txt
#SBATCH -p exbio-cpu
#SBATCH -c 10
#SBATCH --mem 100G
#SBATCH -t 2-12:00:00 # days-hh:mm:ss

Rscript 06-ml-global.R
