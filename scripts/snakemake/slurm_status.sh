#!/bin/bash
# Snakemake cluster-status script for SLURM.
# Called with the SLURM job ID; prints "running", "success", or "failed".

JOBID="$1"
STATUS=$(sacct -j "$JOBID" --noheader --format=State --parsable2 2>/dev/null | head -1)

case "$STATUS" in
    COMPLETED)          echo "success" ;;
    FAILED|CANCELLED*|TIMEOUT|NODE_FAIL|OUT_OF_MEMORY) echo "failed" ;;
    PENDING|RUNNING|COMPLETING) echo "running" ;;
    *)                  echo "running" ;;  # unknown → assume still running
esac
