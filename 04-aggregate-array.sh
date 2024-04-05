#!/bin/bash
#SBATCH -J AGG_CHIP
#SBATCH --error error_logs/%x.%A_%a.%N.%j.txt
#SBATCH --output logs/%x.%A_%a.%N.%j.txt
#SBATCH -p shared-cpu
#SBATCH -c 1
#SBATCH --mem 10G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quirin.manz@tum.de
#SBATCH -t 0-0:30:00 # days-hh:mm:ss
#SBATCH --array=0-2429
# --array=0-809
# array for chip:
readarray -t my_array < <(ls /nfs/data3/IHEC/ChIP-Seq/*pval*)
#readarray -t my_array < <(ls /nfs/data/IHEC/RNAseq/WGBS/*.bw)
# my_array=("/nfs/data3/IHEC/ChIP-Seq/ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000001.4.02b605e5-4373-4009-a925-f6e4537147f5.pval.signal.bigwig")

bigWigAverageOverBed ${my_array[${SLURM_ARRAY_TASK_ID}]} aggregateOver.bed sample_dts/$(basename ${my_array[${SLURM_ARRAY_TASK_ID}]}).tab -minMax && gzip -f sample_dts/$(basename ${my_array[${SLURM_ARRAY_TASK_ID}]}).tab
