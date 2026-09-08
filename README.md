# IHEC-AS

This is the code repository for [Revisiting Evidence for Epigenetic Control of Alternative Splicing](https://doi.org/10.1101/2024.08.30.610315). Since the code in this repository conducts genome-wide analyses of 379 reference epigenomes (415 RNA-seq experiments), there is no small demo dataset available, but after running the first script, a file with the respective paths on the EpiATLAS FTP server is written, which can be subsequently downloaded.

Most analyses are done in R using the .Rmd files in this repo. To keep track of R package versions, we use `renv`. To restore this project's versions, which are documented in the [renv.lock](renv.lock), use `renv::restore()`.
Some analyses use other languages. For those, we have a mamba/conda environment with the documented versions in [env.yml](env.yml) that you can restore using `mamba env create -f env.yml`.

## Running the pipeline

The full pipeline is managed by [Snakemake](https://snakemake.readthedocs.io/) (installed in the `ihec-as` mamba environment). Steps are numbered 01–11; `08-*` global ML models are replaced by the `splicing_ml` Python package.
Lettered stages (`02-3b`, `03b`, `05c`, `09p`, …) are not appendices — they sit at their true position
in the DAG, and the letter only records that they were added after the integer numbering was fixed.

**Two cohort sizes appear in this repository and they are not interchangeable.** The **modelled** cohort
is **379 reference epigenomes / 415 RNA-seq experiments** — quote this for anything about what was
analysed. The **QC/metadata** cohort is **405 / 441**, a clean superset: 26 experiments have QC metadata
but were never modelled (absent from the rMATS quantification). Some QC outputs name the 405/441 scope
in their own headers, correctly for those files.

Stage order (restructured 2026-07-06):

```
01-gather-data
  → 02-1-transcript-filters      (transcript filters; its SUPPA2 chunks are REFERENCE-ONLY --
                                  rMATS-turbo is the PSI-quantification tool, at 02-3)
  → 02-2-rnaseq-normalisation    (GeTMM + DESeq2 vst gene expression)
  → 02-3-rmats-event-filtering   (Procedure-2 + VST gene-expression gate)
      ├─ 02-3b-atlas-summary         (resource counts: 1,522 uuids / 745,545 events)
      ├─ 02-3c-atlas-build-full      (the five-class data release)
      │    → 02-3d-atlas-landscape   (small summaries the resource figure reads)
      └─ 02-4-irfinder-concordance   (independent RI caller, supplement)
  → 03-prepare-aggregation
      └─ 03b-psi-variance-decomposition  (fit-free PSI variance budget; reads no epigenetic data)
  → 04-*  (WGBS / ChIP / MaxEntScan / Pangolin)
  → 05-create-aggregated-dt
      ├─ 05b-feature-pca-sanity
      └─ 05c-cohort-counts  → 05d-supplementary-qc-table
  → 06-correlation / 06c-paired-protocol / 07-2-ml-global-comparison  (+ splicing_ml)
  → 09  event models, two tiers:
      09-1 (feature tables) → 09s-ridge-screen (one job per event, ~34k)
        → 09s-aggregate (FDR)
             ├─ 09f-floor-sample → 09f-floor-score  (detection-power sweep)
             ├─ 09p-permute-groups                  (negative control)
             └─ 09zz (elastic net, hits only) → 09-2 (report)
  → 10-experimental-events
  → 11-paper-figures   TERMINAL: the sole source for every figure and number the manuscript quotes
```

### Per-transcript-filter architecture

Every stage from `02-2` onward runs **once per `transcript_filter`**. Snakemake supplies the
active filter as the `{transcript_filter}` wildcard (exported to R/Python scripts as the
`TRANSCRIPT_FILTER` env var; interactive runs fall back to
`getOption("EpiATLAS_AS_PRIMARY_FILTER", "biotype_filtered")`). All per-filter outputs carry a
`_{transcript_filter}` suffix (e.g. `processed_data/aggregated_dt_filtered_biotype_filtered.csv.gz`).

Because each filter is an independent DAG branch, the filters **build in parallel** and
rebuilding one filter leaves the others untouched. Only the cohort-wide roots (`01`, `02-1`) run
once and are shared across filters.

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

Four keys govern how much ML work `rule all` commits to — worth checking before a long run:

| key | default | effect |
|---|---|---|
| `ml_tasks` | `[classification]` | Tasks built by `rule all`. Regression is off by default; nothing downstream consumes it (`07-2` reads only `*_results_classification.json.gz`). The per-task rules still exist, so a regression output remains buildable as an explicit target — add `regression` here to restore it to `rule all`. |
| `variabilities` | `[High, Low, both]` | Strata for the **main** ML battery. |
| `ablation_variabilities` | `[both]` | Strata for the **ablation** battery only, deliberately narrower — the battery asks "sequence-only vs epigenetics-only", which does not need the High/Low split. `group_col` is not an ablation axis at all: those output paths hardcode `seqnames`, because `ontology` grouping leaks event identity. |
| `shap_max_samples` | `100000` | Outer-test rows per fold fed to SHAP. Cost is linear in rows; measured over 110 real TreeExplainer folds, the median fold costs ~13 s but the worst reaches ~30 min at 50k, so the tail governs. `0` means every row (exhaustive) — roughly 20 h/fold on an SE fold, so reserve it for cases where exact per-sample values are required. Linear-model SHAP is unaffected: it is closed-form and costs ~0 regardless. |

`output_level` must stay `diagnostics` (its default) for SHAP, per-fold predictions and transformed-feature
samples to be stored at all; `compact` discards them regardless of the settings above.

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
	--data-path processed_data/aggregated_dt_filtered_biotype_filtered.csv.gz \
	--output-dir splicing_ml/output/ml_splicing_outputs \
	--wandb --wandb-project splicing-ml
```

Useful W&B detail controls:

- `--wandb-no-fold-table`: disable fold comparison table logging
- `--wandb-no-tuning-details`: disable structured tuning logs
- `--wandb-no-baseline-metrics`: disable baseline and delta logs
- `--wandb-fold-subruns`: optional per-fold sub-runs (higher run/API volume)

## Local run artifacts (written with or without W&B)

Two artifacts are written to disk regardless of tracking, so a run remains inspectable with no W&B
account and after `--no-wandb`:

- `logs/hp_search_{config_key}.jsonl` — one line per (fold, model) with scores, train scores and the
  selected hyperparameters. Keyed by config (`{event_type}-{transcript_filter}-{variability}-{group_col}-{task}`,
  plus a `-feat_…` suffix for ablation runs), the same key used for the W&B child-run name, so a local
  line can be matched to its hosted run.
- `qc/model_input_sanity_{config_key}.json` — what the model actually receives *after* preprocessing:
  shape, per-feature min/max/mean/sd, non-finite counts, one-hot column count vs expected, and
  zero-variance columns. NaN or Inf here is treated as a failure rather than a data quirk, because the
  numeric branch median-imputes and the categorical branch one-hot-encodes, so nothing missing should
  survive; zero variance is only a warning.

HTML reports additionally render a SHAP importance panel per model (xgb/lgbm via `TreeExplainer`,
linear via `LinearExplainer`) plus a cross-model **rank**-agreement table. SHAP magnitudes are
comparable only within a model — different tree families assign systematically different absolute
magnitudes to the same ranking — so only rankings are compared across models.
