#!/bin/bash
#SBATCH -J ML_CFG
#SBATCH --error splicing_ml/output/ml_error_logs/%x.%A_%a.%N.%j.txt
#SBATCH --output splicing_ml/output/ml_logs/%x.%A_%a.%N.%j.txt
#SBATCH -p shared-gpu
#SBATCH --gres=gpu:a40:1
#SBATCH --qos=limitgpus
# QOS limitgpus limits: gres/gpu=4, cpu=80 across all running jobs.
# Effective per-tier concurrency (CPU-bound): L=2 (2x32=64), M=3 (3x24=72), S=4 (4x16=64, GPU-bound).
#SBATCH -c 24
#SBATCH --mem 24G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quirin.manz@tum.de
#SBATCH -t 0-08:00:00
#SBATCH --array=0-0%20
#
# Resource tiers (submit_ml_config_array.sh overrides -c/--mem/--gres/--time per tier):
#   Tier L  SE + both variability    : --mem=40G  -c 32  --gres=gpu:a40:1  --time=0-12:00:00
#   Tier M  SE + High/Low | RI + both: --mem=24G  -c 24  --gres=gpu:a40:1  --time=0-08:00:00  (defaults above)
#   Tier S  RI + High/Low            : --mem=16G  -c 16  --gres=gpu:a40:1  --time=0-06:00:00
# Note: only A40 (48 GB VRAM) is used; Titan GPUs on gpu01 have ~12 GB VRAM and 92 GB node RAM,
# which is insufficient for cuML SVM on even the smallest tier's dataset sizes.
# (outer_splits=10 original estimates: L=2-00:00:00, M=1-12:00:00, S=1-00:00:00)
#
# Rationale:
#   Every job loads the full dataset (~6.4 GB pandas RSS).  Larger subsets add
#   up to ~1.4 GB of feature matrix which gets copied ~3-4x during nested CV.
#   SE/both configs reach 3.9 M rows; RI/High configs stay under 130 K rows.
#   Walltime driven by LogisticRegressionCV (elasticnet, CPU-only): with
#   outer_splits=5 (inner_splits=4, 5 l1_ratios -> 20 tasks) ~10 min/fold at 16c.
#   MLP excluded from default models.

set -euo pipefail

# Usage (must be run inside environment with splicing_ml installed):
# sbatch --array=0-(N-1)%20 scripts/slurm_ml_one_config.sh <config_tsv> <data_path> <output_root> [env_name] [wandb_project]
#
# config_tsv format (tab-separated, no header):
#   event_type\ttranscript_filter\tvariability\tgroup_col\t[task]
#
# task: classification | regression | both (default: both)
#
# Pass a wandb_project name to enable W&B tracking. The script will call
# `wandb login` so stored credentials are loaded in the job environment.
# Optional env toggles:
#   WANDB_REQUIRE_AUTH=0           -> add --wandb-no-require-auth
#   WANDB_FOLD_SUBRUNS=1           -> add --wandb-fold-subruns
#   WANDB_NO_FOLD_TABLE=1          -> add --wandb-no-fold-table
#   WANDB_NO_TUNING_DETAILS=1      -> add --wandb-no-tuning-details
#   WANDB_NO_BASELINE_METRICS=1    -> add --wandb-no-baseline-metrics

CONFIG_TSV="${1:?missing config TSV path}"
DATA_PATH="${2:?missing data path}"
OUTPUT_ROOT="${3:?missing output root}"
ENV_NAME="${4:-}"
WANDB_PROJECT="${5:-}"

if [[ -n "$ENV_NAME" ]]; then
  # SLURM does not source .bashrc, so load the environment module that
  # provides mamba before trying to activate the conda environment.
  module load miniforge3/24.7.1
  eval "$(conda shell.bash hook)"
  conda activate "$ENV_NAME"
fi

readarray -t CONFIG_LINES < "$CONFIG_TSV"
CFG_LINE="${CONFIG_LINES[${SLURM_ARRAY_TASK_ID}]}"
IFS=$'\t' read -r EVENT_TYPE TRANSCRIPT_FILTER VARIABILITY GROUP_COL TASK <<< "$CFG_LINE"
TASK="${TASK:-both}"

if [[ -z "${GROUP_COL:-}" ]]; then
  echo "[ERROR] Missing group_col in config TSV line for task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

OUT_DIR="${OUTPUT_ROOT}/${EVENT_TYPE}_${TRANSCRIPT_FILTER}_${VARIABILITY}_${GROUP_COL}"
mkdir -p "$OUT_DIR"

echo "[SLURM] task=${SLURM_ARRAY_TASK_ID} config=${EVENT_TYPE}/${TRANSCRIPT_FILTER}/${VARIABILITY}/${GROUP_COL} ml_task=${TASK}"
echo "[SLURM] output=${OUT_DIR}"

# Build command as an array to avoid line-continuation/comment parsing pitfalls.
CMD=(
  python -u run_splicing_ml.py
  --debug
  --tune-threshold
  --optuna
  --data-path "$DATA_PATH"
  --output-dir "$OUT_DIR"
  --only-event-type "$EVENT_TYPE"
  --only-transcript-filter "$TRANSCRIPT_FILTER"
  --only-variability "$VARIABILITY"
  --only-group "$GROUP_COL"
)

case "$TASK" in
  classification) CMD+=(--skip-regression) ;;
  regression)     CMD+=(--skip-classification) ;;
  both)           ;;  # run both — no skip flag needed
  *) echo "[ERROR] Unknown task '${TASK}'; expected classification|regression|both" >&2; exit 1 ;;
esac
echo "[SLURM] ml_task=${TASK} (skip flags applied if any)"

# Use one fewer core than allocated to leave a small overhead margin.
if [[ -n "${SLURM_CPUS_PER_TASK:-}" && "${SLURM_CPUS_PER_TASK}" -gt 1 ]]; then
  MAX_CORES="${SLURM_CPUS_PER_TASK}" #((SLURM_CPUS_PER_TASK - 1))
  CMD+=(--max-cores "$MAX_CORES")
  echo "[SLURM] max_cores=${MAX_CORES} (from SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK})"
fi

if [[ -n "$WANDB_PROJECT" ]]; then
  echo "[SLURM] Loading W&B stored credentials"
  wandb login >/dev/null
  CMD+=(--wandb --wandb-project "$WANDB_PROJECT")
  if [[ "${WANDB_REQUIRE_AUTH:-1}" == "0" ]]; then
    CMD+=(--wandb-no-require-auth)
  fi
  if [[ "${WANDB_FOLD_SUBRUNS:-0}" == "1" ]]; then
    CMD+=(--wandb-fold-subruns)
  fi
  if [[ "${WANDB_NO_FOLD_TABLE:-0}" == "1" ]]; then
    CMD+=(--wandb-no-fold-table)
  fi
  if [[ "${WANDB_NO_TUNING_DETAILS:-0}" == "1" ]]; then
    CMD+=(--wandb-no-tuning-details)
  fi
  if [[ "${WANDB_NO_BASELINE_METRICS:-0}" == "1" ]]; then
    CMD+=(--wandb-no-baseline-metrics)
  fi
  echo "[SLURM] W&B enabled: project=${WANDB_PROJECT}"
fi

echo "[SLURM] running: ${CMD[*]}"
"${CMD[@]}"
