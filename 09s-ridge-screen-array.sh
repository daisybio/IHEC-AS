#!/bin/bash
#SBATCH --job-name=ridge_screen
#SBATCH --output=event_glmnet_logs/%x_%a_%A_%N.log
#SBATCH --error=event_glmnet_logs/%x_%a_%A_%N.err
#SBATCH -p shared-cpu
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=04:00:00
#SBATCH --exclude=compms-gpu-1.exbio.wzw.tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quirin.manz@tum.de

# Tier-1 fit-free ridge screen — one event per array task.
# Usage: sbatch --array=0-<N-1>%<throttle> 09s-ridge-screen-array.sh IDS CFG PROJECT [OFFSET]
# OFFSET (default 0) lets 09s-dispatch.R chunk an id list larger than SLURM's
# MaxArraySize into several array submissions: global line = OFFSET + task + 1.
IDS_FILE="$1"
CFG_FILE="$2"
PROJECT_DIR="$3"
OFFSET="${4:-0}"

source /etc/profile.d/modules.sh
module load r/4.2.1

# SLURM_ARRAY_TASK_ID is 0-indexed; sed line numbers are 1-indexed
line=$((OFFSET + SLURM_ARRAY_TASK_ID + 1))
id=$(sed -n "${line}p" "$IDS_FILE")
[ -z "$id" ] && { echo "No id at line $line (past end of $IDS_FILE) — skipping"; exit 0; }

cd "$PROJECT_DIR"
Rscript "${PROJECT_DIR}/09s-ridge-screen.R" "$CFG_FILE" "$id" 4
