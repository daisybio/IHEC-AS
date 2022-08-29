#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=GET-TF
#SBATCH --output=%x.%j.txt
#SBATCH --mem-per-cpu=2G

Rscript 07-event-specific-get-tf.R
