#!/bin/bash
set -euo pipefail

# Smart submit helper for one-job-per-config execution.
#
# Generates only the subset configs that are actually present in the data, then
# submits SLURM arrays tiered by estimated resource requirements.
#
# Resource tiers (based on measured data characteristics):
#   Tier L  SE + both variability     2.7-3.9 M rows  --mem=40G  -c 16  24 h
#   Tier M  SE + High/Low | RI + both 0.4-2.5 M rows  --mem=24G  -c 12  12 h
#   Tier S  RI + High/Low             55K-840K rows    --mem=16G  -c  8   6 h
#
# Every job loads the full ~6.4 GB dataset (pandas), so 16 GB is the floor.
#
# Usage:
#   bash scripts/submit_ml_config_array.sh [data_path] [output_root] [env_name] [max_parallel] [wandb_project]
#
# Omit wandb_project (or pass "") to disable W&B tracking.
# When provided, this script runs `wandb login` after activating the
# environment so stored credentials are loaded.

DATA_PATH="${1:-processed_data/aggregated_dt_filtered.csv.gz}"
OUTPUT_ROOT="${2:-processed_data/slurm_ml_outputs}"
ENV_NAME="${3:-ihec-as}"
MAX_PARALLEL="${4:-20}"
WANDB_PROJECT="${5:-}"

mkdir -p "$OUTPUT_ROOT"

CONFIG_L="scripts/configs_tier_L.tsv"
CONFIG_M="scripts/configs_tier_M.tsv"
CONFIG_S="scripts/configs_tier_S.tsv"

# Activate mamba environment on submit node before generating config TSV and sbatch.
if ! command -v mamba >/dev/null 2>&1; then
  echo "[ERROR] mamba command not found in current shell." >&2
  echo "[ERROR] Initialize mamba in your shell, then rerun this script." >&2
  exit 127
fi

eval "$(mamba shell hook --shell bash)"
mamba activate "$ENV_NAME"

if [[ -n "$WANDB_PROJECT" ]]; then
  echo "[W&B] Loading stored credentials for project '${WANDB_PROJECT}'"
  wandb login >/dev/null
fi

export DATA_PATH CONFIG_L CONFIG_M CONFIG_S

python - <<'PY'
import os
from splicing_ml.data import generate_subset_configs, load_dataset

DATA_PATH = os.environ["DATA_PATH"]
CONFIG_L   = os.environ["CONFIG_L"]
CONFIG_M   = os.environ["CONFIG_M"]
CONFIG_S   = os.environ["CONFIG_S"]

df = load_dataset(DATA_PATH)
configs = generate_subset_configs(df)

def tier(cfg):
    if cfg.event_type == "SE" and cfg.variability == "both":
        return "L"
    if cfg.event_type == "RI" and cfg.variability not in {"both"}:
        return "S"
    return "M"

buckets = {"L": [], "M": [], "S": []}
for cfg in configs:
    buckets[tier(cfg)].append(cfg)

for key, path in [("L", CONFIG_L), ("M", CONFIG_M), ("S", CONFIG_S)]:
    with open(path, "w", encoding="utf-8") as f:
        for cfg in buckets[key]:
            for task in ("classification", "regression"):
                f.write(f"{cfg.event_type}\t{cfg.transcript_filter}\t{cfg.variability}\t{cfg.group_col}\t{task}\n")
    n_rows = len(buckets[key]) * 2
    print(f"Tier {key}: {n_rows} jobs ({len(buckets[key])} configs x 2 tasks) -> {path}")
PY

submit_tier() {
  local tier_label="$1"
  local config_tsv="$2"
  local mem="$3"
  local cpus="$4"
  local walltime="$5"

  if [[ ! -s "$config_tsv" ]]; then
    echo "[SKIP] Tier ${tier_label}: no configs found"
    return
  fi
  local n_cfg
  n_cfg=$(wc -l < "$config_tsv")

  local array_spec="0-$((n_cfg - 1))%${MAX_PARALLEL}"
  echo "[SUBMIT] Tier ${tier_label}: ${n_cfg} configs, array=${array_spec}, mem=${mem}, cpus=${cpus}, time=${walltime}"

  sbatch \
    --export=ALL \
    --mem="${mem}" \
    -c "${cpus}" \
    --time="${walltime}" \
    --array="${array_spec}" \
    --job-name="ML_${tier_label}" \
    scripts/slurm_ml_one_config.sh \
    "$config_tsv" "$DATA_PATH" "$OUTPUT_ROOT" "$ENV_NAME" "$WANDB_PROJECT"
}

submit_tier "L" "$CONFIG_L" "40G" "16" "0-24:00:00"
submit_tier "M" "$CONFIG_M" "24G" "12" "0-12:00:00"
submit_tier "S" "$CONFIG_S" "16G"  "8" "0-06:00:00"
