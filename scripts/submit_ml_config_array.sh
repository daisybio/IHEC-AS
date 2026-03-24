#!/bin/bash
set -euo pipefail

# Smart submit helper for one-job-per-config execution.
#
# Generates only the subset configs that are actually present in the data, then
# submits a SLURM array where each task processes exactly one config.
#
# Usage:
#   bash scripts/submit_ml_config_array.sh [data_path] [output_root] [env_name] [max_parallel]

DATA_PATH="${1:-processed_data/aggregated_dt_filtered.csv.gz}"
OUTPUT_ROOT="${2:-processed_data/slurm_ml_outputs}"
ENV_NAME="${3:-ihec-as}"
MAX_PARALLEL="${4:-20}"

CONFIG_TSV="processed_data/slurm_ml_config_list.tsv"
mkdir -p "$(dirname "$CONFIG_TSV")"
mkdir -p "$OUTPUT_ROOT"

# Activate mamba environment on submit node before generating config TSV and sbatch.
if ! command -v mamba >/dev/null 2>&1; then
  echo "[ERROR] mamba command not found in current shell." >&2
  echo "[ERROR] Initialize mamba in your shell, then rerun this script." >&2
  exit 127
fi

eval "$(mamba shell hook --shell bash)"
mamba activate "$ENV_NAME"

export DATA_PATH
export CONFIG_TSV

python - <<'PY'
import os
from splicing_ml.data import generate_subset_configs, load_dataset

DATA_PATH = os.environ["DATA_PATH"]
CONFIG_TSV = os.environ["CONFIG_TSV"]

df = load_dataset(DATA_PATH)
configs = generate_subset_configs(df)

with open(CONFIG_TSV, "w", encoding="utf-8") as f:
    for cfg in configs:
        f.write(
            f"{cfg.event_type}\t{cfg.transcript_filter}\t{cfg.variability}\t{cfg.group_col}\n"
        )

print(f"Wrote {len(configs)} configs to {CONFIG_TSV}")
PY

N_CFG=$(wc -l < "$CONFIG_TSV")
if [[ "$N_CFG" -lt 1 ]]; then
  echo "No configs found in $DATA_PATH"
  exit 1
fi

ARRAY_SPEC="0-$((N_CFG-1))%${MAX_PARALLEL}"

echo "Submitting ${N_CFG} configs with array=${ARRAY_SPEC}"

sbatch --array="$ARRAY_SPEC" scripts/slurm_ml_one_config.sh "$CONFIG_TSV" "$DATA_PATH" "$OUTPUT_ROOT" "$ENV_NAME"
