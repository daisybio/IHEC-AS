#!/bin/bash
#
#SBATCH -c 20
#SBATCH -J AGG
#SBATCH --output=%x.%j.%N.txt
#SBATCH --mem=380000
#SBATCH -t 2-00:00:00 # days-hh:mm:ss

Rscript 04-preliminary-results-aggregation.R
