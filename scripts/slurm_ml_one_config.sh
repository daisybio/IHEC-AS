#!/bin/bash
#SBATCH -J ML_CFG
#SBATCH --error splicing_ml/output/ml_error_logs/%x.%A_%a.%N.%j.txt
#SBATCH --output splicing_ml/output/ml_logs/%x.%A_%a.%N.%j.txt
#SBATCH -p shared-gpu
#SBATCH --gres=gpu:1
#SBATCH --qos=limitgpus
#SBATCH -c 8
#SBATCH --mem 12G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quirin.manz@tum.de
#SBATCH -t 1-00:00:00
#SBATCH --array=0-0%20

set -euo pipefail

# Usage (must be run inside environment with splicing_ml installed):
# sbatch --array=0-(N-1)%20 scripts/slurm_ml_one_config.sh <config_tsv> <data_path> <output_root> [env_name]
#
# config_tsv format (tab-separated, no header):
#   event_type\ttranscript_filter\tvariability\tgroup_col

CONFIG_TSV="${1:?missing config TSV path}"
DATA_PATH="${2:?missing data path}"
OUTPUT_ROOT="${3:?missing output root}"
ENV_NAME="${4:-}"

if [[ -n "$ENV_NAME" ]]; then
  if ! command -v mamba >/dev/null 2>&1; then
    echo "[ERROR] mamba command not found on compute node while ENV_NAME was provided" >&2
    exit 127
  fi
  eval "$(mamba shell hook --shell bash)"
  mamba activate "$ENV_NAME"
fi

readarray -t CONFIG_LINES < "$CONFIG_TSV"
CFG_LINE="${CONFIG_LINES[${SLURM_ARRAY_TASK_ID}]}"
IFS=$'\t' read -r EVENT_TYPE TRANSCRIPT_FILTER VARIABILITY GROUP_COL <<< "$CFG_LINE"

if [[ -z "${GROUP_COL:-}" ]]; then
  echo "[ERROR] Missing group_col in config TSV line for task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

OUT_DIR="${OUTPUT_ROOT}/${EVENT_TYPE}_${TRANSCRIPT_FILTER}_${VARIABILITY}_${GROUP_COL}"
mkdir -p "$OUT_DIR"

echo "[SLURM] task=${SLURM_ARRAY_TASK_ID} config=${EVENT_TYPE}/${TRANSCRIPT_FILTER}/${VARIABILITY}/${GROUP_COL}"
echo "[SLURM] output=${OUT_DIR}"

# Build command as an array to avoid line-continuation/comment parsing pitfalls.
CMD=(
  python -u run_splicing_ml.py
  --debug
  --data-path "$DATA_PATH"
  --output-dir "$OUT_DIR"
  --only-event-type "$EVENT_TYPE"
  --only-transcript-filter "$TRANSCRIPT_FILTER"
  --only-variability "$VARIABILITY"
  --only-group "$GROUP_COL"
)

# Use one fewer core than allocated to leave a small overhead margin.
if [[ -n "${SLURM_CPUS_PER_TASK:-}" && "${SLURM_CPUS_PER_TASK}" -gt 1 ]]; then
  MAX_CORES=$((SLURM_CPUS_PER_TASK - 1))
  CMD+=(--max-cores "$MAX_CORES")
  echo "[SLURM] max_cores=${MAX_CORES} (from SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK})"
fi

echo "[SLURM] running: ${CMD[*]}"
"${CMD[@]}"
