#!/bin/bash
#SBATCH -J TF
#SBATCH --output %x.%j.%N.txt
#SBATCH -p exbio-cpu
#SBATCH -c 1
#SBATCH --mem 2G
#SBATCH -t 4-00:00:00 # days-hh:mm:ss

srun Rscript 07-event-specific-get-tf.R
