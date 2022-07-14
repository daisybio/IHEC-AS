#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=IHEC-IA-AS
#SBATCH --output=%x.%j.txt
#SBATCH --mem-per-cpu=30G

Rscript 04-preliminary-results-aggregation.R
