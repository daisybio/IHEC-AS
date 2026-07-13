#!/bin/bash
#SBATCH --job-name=event_glmnet
#SBATCH --output=event_glmnet_logs/%x_%a_%A_%N.log
#SBATCH --error=event_glmnet_logs/%x_%a_%A_%N.err
#SBATCH -p shared-cpu
#SBATCH --mem=4G
#SBATCH -c 1
#SBATCH --time=12:00:00
#SBATCH --exclude=compms-gpu-1.exbio.wzw.tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quirin.manz@tum.de

IDS_FILE="$1"
CFG_FILE="$2"
PROJECT_DIR="$3"

source /etc/profile.d/modules.sh
module load r/4.2.1

# SLURM_ARRAY_TASK_ID is 0-indexed; sed line numbers are 1-indexed
id=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$IDS_FILE")

cd "$PROJECT_DIR"
Rscript "${PROJECT_DIR}/09zz-ml-event-glmnet-tidymodels.R" "$CFG_FILE" "$id"
