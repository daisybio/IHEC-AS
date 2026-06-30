# IHEC-AS

This is the code repository for [Revisiting Evidence for Epigenetic Control of Alternative Splicing](https://doi.org/10.1101/2024.08.30.610315). Since the code in this repository conducts genome-wide analyses of 405 reference epigenomes, there is no small demo dataset available, but after running the first script, a file with the respective paths on the EpiATLAS FTP server is written, which can be subsequently downloaded.

Most analyses are done in R using the .Rmd files in this repo. To keep track of R package versions, we use `renv`. To restore this project's versions, which are documented in the [renv.lock](renv.lock), use `renv::restore()`.
Some analyses use other languages. For those, we have a mamba/conda environment with the documented versions in [env.yml](env.yml) that you can restore using `mamba env create -f env.yml`.

## Running the pipeline

The full pipeline is managed by [Snakemake](https://snakemake.readthedocs.io/) (installed in the `ihec-as` mamba environment). Steps are numbered 01–10; `08-*` global ML models are replaced by the `splicing_ml` Python package.

### Dry run (check DAG, no execution)

```bash
mamba run -n ihec-as snakemake --profile profiles/slurm -n
```

### Full run on SLURM

```bash
mamba run -n ihec-as snakemake --profile profiles/slurm
```

This submits all jobs to SLURM automatically. CPU and GPU rules are dispatched to the correct partitions via `scripts/snakemake/slurm_submit.sh`.

### Run a single transcript filter

```bash
mamba run -n ihec-as snakemake --profile profiles/slurm \
    --config transcript_filters='["biotype_filtered"]'
```

### Run a specific target

```bash
# Re-run correlation for one filter only
mamba run -n ihec-as snakemake --profile profiles/slurm \
    processed_data/correlation_intrinsic_biotype_filtered.csv.gz

# Re-run all ML configs for one filter
mamba run -n ihec-as snakemake --profile profiles/slurm \
    --config transcript_filters='["biotype_filtered"]' \
    $(mamba run -n ihec-as snakemake --profile profiles/slurm \
        --config transcript_filters='["biotype_filtered"]' -n --quiet 2>/dev/null \
        | grep splicing_ml | awk '{print $NF}')
```

### Force re-run a rule

```bash
mamba run -n ihec-as snakemake --profile profiles/slurm \
    --forcerun create_aggregated_dt
```

### Configuration

Pipeline parameters (paths, resource budgets, filter lists) are in [`config/snakemake_config.yaml`](config/snakemake_config.yaml).  
Local machine parameters go in [`.Rprofile`](.Rprofile).

## Environment Setup

### CPU servers

Use the CPU environment file:

```bash
mamba env create -f env.cpu.yml
mamba run -n ihec-as ./run_splicing_ml.py --check-gpu
```

Expected check output on non-GPU hosts:

```text
XGBoost : cpu (no GPU)
cuML SVM: cpu (no GPU)
PyTorch : cpu
```

### GPU servers (cuML + torch compatibility)

When creating the GPU environment from [env.yml](env.yml), pip installs torch CUDA wheels that can shadow RAPIDS CUDA libraries needed by cuML. The following post-create cleanup is required.

```bash
mamba env create -f env.yml
mamba run -n ihec-as python -m pip uninstall -y nvidia-cublas-cu12 nvidia-nvjitlink-cu12 nvidia-cusparse-cu12
mamba run -n ihec-as python -c "import cuml; print('cuml', cuml.__version__)"
mamba run -n ihec-as python run_splicing_ml.py --check-gpu
```

Verified expected output after cleanup on a GPU host:

```text
cuml 26.04.00a122

XGBoost : cuda
cuML SVM: cuda
PyTorch : cuda
```

## W&B Tracking

The ML pipeline logs to Weights & Biases when `--wandb` is enabled.

Default hierarchy:

- 1 parent run per pipeline launch (`job_type=orchestrator`)
- 1 child run per subset config and task (`job_type=subset-task`)

Default visibility improvements:

- all outer-fold/model metrics are logged
- fold-level comparison table in each child run (`folds/by_model`)
- structured tuning summaries (best score and parameter summary)
- baseline and delta-vs-baseline metrics

Strict authentication is enabled by default for reproducibility:

- if `--wandb` is used and `WANDB_API_KEY` is missing, the run fails fast
- disable strict behavior with `--wandb-no-require-auth` (falls back to NullTracker)

Example:

```bash
mamba run -n ihec-as python run_splicing_ml.py \
	--data-path processed_data/aggregated_dt_filtered.csv.gz \
	--output-dir splicing_ml/output/ml_splicing_outputs \
	--wandb --wandb-project splicing-ml
```

Useful W&B detail controls:

- `--wandb-no-fold-table`: disable fold comparison table logging
- `--wandb-no-tuning-details`: disable structured tuning logs
- `--wandb-no-baseline-metrics`: disable baseline and delta logs
- `--wandb-fold-subruns`: optional per-fold sub-runs (higher run/API volume)
