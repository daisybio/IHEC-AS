#!/bin/bash
#SBATCH -J AGG_CHIP
#SBATCH --output logs/%x.%A_%a.%N.%j.txt
#SBATCH -p shared-cpu
#SBATCH -c 1
#SBATCH --mem 1G
#SBATCH -t 0-0:30:00 # days-hh:mm:ss
#SBATCH --array=0-2429
# --array=0-809
# array for chip:
readarray -t my_array < <(ls /nfs/data3/IHEC/ChIP-Seq/*pval*)
#readarray -t my_array < <(ls /nfs/data/IHEC/RNAseq/WGBS/*.bw)

srun bigWigAverageOverBed ${my_array[${SLURM_ARRAY_TASK_ID}]} aggregateOver.bed sample_dts/$(basename ${my_array[${SLURM_ARRAY_TASK_ID}]}).tab -minMax && gzip -f sample_dts/$(basename ${my_array[${SLURM_ARRAY_TASK_ID}]}).tab
