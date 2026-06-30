#!/bin/bash
# Snakemake --cluster submission wrapper.
# Called as: slurm_submit.sh MEM_MB THREADS RUNTIME GPU RULE JOBID JOBSCRIPT
# Snakemake appends the jobscript as the last positional argument.

MEM_MB="${1:?mem_mb required}"
THREADS="${2:?threads required}"
RUNTIME="${3:?runtime required}"
GPU="${4:?gpu flag required}"
RULE="${5:?rule required}"
JOBID="${6:?jobid required}"
SCRIPT="${@: -1}"  # jobscript is always the last arg

mkdir -p logs/slurm reports

CMD=(
    sbatch --parsable
    "--mem=${MEM_MB}M"
    "--cpus-per-task=${THREADS}"
    "--time=${RUNTIME}"
    "--job-name=smk-${RULE}"
    "--output=logs/slurm/${RULE}_${JOBID}.out"
    "--error=logs/slurm/${RULE}_${JOBID}.err"
    "--mail-type=FAIL"
    "--mail-user=quirin.manz@tum.de"
)

if [[ "${GPU}" == "1" ]]; then
    CMD+=(
        "--partition=shared-gpu"
        "--gres=gpu:a40:1"
        "--qos=limitgpus"
        "--exclude=compms-gpu-1.exbio.wzw.tum.de"
    )
else
    CMD+=("--partition=exbio-interactive")
fi

CMD+=("$SCRIPT")
exec "${CMD[@]}"
