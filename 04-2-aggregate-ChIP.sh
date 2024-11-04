#!/bin/bash
#SBATCH -J AGG_CHIP
#SBATCH --error sample_dts/error_logs/%x.%A_%a.%N.%j.txt
#SBATCH --output sample_dts/logs/%x.%A_%a.%N.%j.txt
#SBATCH -p shared-cpu
#SBATCH -c 1
#SBATCH --mem 10G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quirin.manz@tum.de
#SBATCH -t 0-1:00:00 # days-hh:mm:ss
#SBATCH --array=0-2429

# array for chip:
readarray -t my_array < <(ls /nfs/data3/IHEC/ChIP-Seq/*pval*)
# make sure that you run this script in the conda environment, such that bigWigAverageOverBed is installed
bigWigAverageOverBed ${my_array[${SLURM_ARRAY_TASK_ID}]} aggregateOver.bed sample_dts/$(basename ${my_array[${SLURM_ARRAY_TASK_ID}]}).tab -minMax && gzip -f sample_dts/$(basename ${my_array[${SLURM_ARRAY_TASK_ID}]}).tab
