#!/bin/bash
set -euo pipefail

# Smart submit helper for one-job-per-config execution.
#
# Generates only the subset configs that are actually present in the data, then
# submits SLURM arrays tiered by estimated resource requirements.
#
# Resource tiers (based on measured data characteristics, MLP excluded from default models):
#   Tier L  SE + both variability     2.7-3.9 M rows  --mem=40G  -c 32  12 h
#   Tier M  SE + High/Low | RI + both 0.4-2.5 M rows  --mem=24G  -c 24   8 h
#   Tier S  RI + High/Low             55K-840K rows    --mem=16G  -c 16   6 h
# (outer_splits=10 original estimates: L=48h, M=36h, S=24h)
#
# Walltime is driven by LogisticRegressionCV (elasticnet, CPU-only) which
# parallelises over inner_folds x n_l1_ratios tasks.  With outer_splits=5
# (default): inner_splits=4, 5 l1_ratios -> 20 tasks.  Measured at 477k train
# samples: ~10 min/fold at 16c (2 batches of 20).  Both tasks x 5 folds gives
# the estimates above.  Adding --models mlp roughly doubles runtime.
#
# Every job loads the full ~6.4 GB dataset (pandas), so 16 GB is the floor.
#
# Usage:
#   bash scripts/submit_ml_config_array.sh [data_path] [output_root] [env_name] [wandb_project]
#
# Omit wandb_project (or pass "") to disable W&B tracking.
# Override per-tier concurrency via env vars: MAX_PARALLEL_L, MAX_PARALLEL_M, MAX_PARALLEL_S.
# When provided, this script runs `wandb login` after activating the
# environment so stored credentials are loaded.

DATA_PATH="${1:-processed_data/aggregated_dt_filtered.csv.gz}"
OUTPUT_ROOT="${2:-splicing_ml/output/slurm_ml_outputs}"
ENV_NAME="${3:-ihec-as}"
WANDB_PROJECT="${4:-}" # splicing-ml
# QOS limitgpus: max 4 GPUs and 80 CPUs concurrently across all running jobs.
# Per-tier limits are CPU-bound: L=2 (2x32c), M=3 (3x24c), S=4 (4x16c, GPU-bound).
# These assume only one tier is active; reduce if submitting multiple tiers simultaneously.
MAX_PARALLEL_L="${MAX_PARALLEL_L:-2}"
MAX_PARALLEL_M="${MAX_PARALLEL_M:-3}"
MAX_PARALLEL_S="${MAX_PARALLEL_S:-4}"

mkdir -p "$OUTPUT_ROOT"

CONFIG_L="scripts/configs_tier_L.tsv"
CONFIG_M="scripts/configs_tier_M.tsv"
CONFIG_S="scripts/configs_tier_S.tsv"

# Activate mamba environment on submit node before generating config TSV and sbatch.
if [[ -n "$ENV_NAME" ]]; then
  # SLURM does not source .bashrc, so load the environment module that
  # provides mamba before trying to activate the conda environment.
  module load miniforge3/24.7.1
  eval "$(conda shell.bash hook)"
  conda activate "$ENV_NAME"
fi

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

df = load_dataset(DATA_PATH, reader_backend = "polars")
configs = generate_subset_configs(df)

def tier(cfg):
    if cfg.event_type == "SE" and cfg.variability == "both":
        return "L"
    if cfg.event_type == "RI" and cfg.variability != "both":
        return "S"
    return "M"

buckets = {"L": [], "M": [], "S": []}
for cfg in configs:
    buckets[tier(cfg)].append(cfg)

for key, path in [("L", CONFIG_L), ("M", CONFIG_M), ("S", CONFIG_S)]:
    with open(path, "w", encoding="utf-8") as f:
        for cfg in buckets[key]:
            if cfg.transcript_filter != "biotype_filtered" or cfg.variability != "both":
                continue
            f.write(f"{cfg.event_type}\t{cfg.transcript_filter}\t{cfg.variability}\t{cfg.group_col}\tclassification\n")
    n_rows = sum(
        1 for cfg in buckets[key]
        if cfg.transcript_filter != "biotype_filtered" and cfg.variability != "both"
    )
    print(f"Tier {key}: {n_rows} jobs -> {path}")
PY

submit_tier() {
  local tier_label="$1"
  local config_tsv="$2"
  local mem="$3"
  local cpus="$4"
  local walltime="$5"
  local max_parallel="$6"

  if [[ ! -s "$config_tsv" ]]; then
    echo "[SKIP] Tier ${tier_label}: no configs found"
    return
  fi
  local n_cfg
  n_cfg=$(wc -l < "$config_tsv")

  local array_spec="0-$((n_cfg - 1))%${max_parallel}"
  echo "[SUBMIT] Tier ${tier_label}: ${n_cfg} configs, array=${array_spec}, mem=${mem}, cpus=${cpus}, time=${walltime}"

  sbatch \
    --export=ALL \
    --mem="${mem}" \
    -c "${cpus}" \
    --time="${walltime}" \
    --array="${array_spec}" \
    --job-name="ML_${tier_label}" \
    --exclude=compms-gpu-1.exbio.wzw.tum.de \
    scripts/slurm_ml_one_config.sh \
    "$config_tsv" "$DATA_PATH" "$OUTPUT_ROOT" "$ENV_NAME" "$WANDB_PROJECT"
}

# Original (outer_splits=10): L=2-00:00:00, M=1-12:00:00, S=1-00:00:00
# GPU: all tiers use A40 (48 GB VRAM). The only other GPU type available is Titan (~12 GB VRAM,
# 92 GB node RAM), which is insufficient for cuML SVM on any tier's dataset sizes.
submit_tier "S" "$CONFIG_S" "32G" "8"  "0-06:00:00" "$MAX_PARALLEL_S" && \
submit_tier "M" "$CONFIG_M" "64G" "8"  "0-12:00:00" "$MAX_PARALLEL_M" && \
submit_tier "L" "$CONFIG_L" "128G" "12" "0-24:00:00" "$MAX_PARALLEL_L"