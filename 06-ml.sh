#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=ML-build
#SBATCH --output=%x.%j.txt
#SBATCH --mem-per-cpu=20G

Rscript 06-ml.R
