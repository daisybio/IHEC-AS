"""
IHEC-AS Snakemake workflow.

Run (dry-run):
    mamba run -n ihec-as snakemake --profile profiles/slurm -n

Run (SLURM):
    mamba run -n ihec-as snakemake --profile profiles/slurm

Run a single filter / several in parallel — each transcript_filter is an
independent DAG branch (shared roots 01, 02-1 run once). Default target =
primary_filter; build others by requesting their targets or extending
config["transcript_filters"]:
    mamba run -n ihec-as snakemake --profile profiles/slurm \
        processed_data/aggregated_dt_filtered_transcripts.csv.gz

Steps (everything from 02-2 on is per-transcript_filter, suffix _{tf}):
    01   gather-data              →  processed_data/file_table.csv.gz          (cohort-wide)
    02-1 transcript-filters       →  splicing_analysis/filtered_transcript_ids.rds, gencode.v29.{tf}.gtf  (cohort-wide)
    02-2 rnaseq-normalisation     →  processed_data/gene_expression_normalised_{tf}.csv.gz  (getmm+vst)
    02-3 rmats-event-filtering    →  splicing_analysis/rmats/{tf}/event_{et}.{psi,jc.csv.gz} (Procedure 2 + VST gate), qc/vst_expression_cutoff_{tf}.csv
    03   prepare-aggregation      →  processed_data/{aggregateOver,keep_rows_manual,sample_cols,ijc_sjc_dt,event_annotations_dt,psi_long_dt,pangolin_events,*ss}_{tf}.*
    04-1 aggregate-WGBS           →  sample_dts/WGBS_agg_{tf}.csv.gz
    04-2 aggregate-ChIP           →  sample_dts/chip_agg_{tf}/… , chip_agg_{tf}.done
    04-3 maxentscan               →  processed_data/{3,3down,5,5up}scores_{tf}.txt
    04-4 pangolin                 →  processed_data/pangolin_scores_{tf}.csv
    04-5 rbp-binding-sites        →  processed_data/rbp_per_event.rds  (dependency of create_aggregated_dt)
    05   create-aggregated-dt     →  processed_data/aggregated_dt_filtered_{tf}.csv.gz
    05b  feature-pca-sanity       →  reports/05b-feature-pca-sanity_{primary}.html
    06   correlation              →  processed_data/correlation_intrinsic_{tf}.csv.gz  (primary filter)
    08   splicing_ml_{classification,regression}          →  splicing_ml/output/{et}_{tf}_{var}_{gc}/ (parallel per task)
    08b  splicing_ml_ablation_{classification,regression} →  splicing_ml/output_ablation/{fg}/{et}_{primary}_{var}_seqnames/ (battery, primary filter, parallel per task)
    09-1 event-models             →  processed_data/event_models/{tf}/.done
    09-2 ml-analysis              →  reports/09-2-ml-local-new_{primary}.html
    10   experimental-events      →  reports/10-experimental-events_{tf}.html

Rule → SLURM: the slurm executor plugin submits each job as its own `sbatch`,
translating threads→--cpus-per-task, resources.mem_mb→--mem, runtime→--time,
slurm_partition/qos/gres→flags, slurm_extra→appended verbatim. The rule shell
(wrapped by shell.prefix below) is the job script. aggregate_chip_one is the
one exception: it is collapsed into a single `sbatch --array` via the profile's
`slurm-array-jobs` setting (see that rule).
"""

from itertools import product
from glob import glob

configfile: "config/snakemake_config.yaml"

# splicing_ml source files: entry script + whole package. Declared as `input:`
# of the splicing_ml rule so any code edit enters Snakemake's rerun-triggers
# (the shell string alone is opaque to change detection).
SPLICING_ML_SRC = ["run_splicing_ml.py"] + sorted(glob("splicing_ml/**/*.py", recursive=True))

CONDA_BASE = "/nfs/data/cluster/software/miniforge3/24.7.1/miniforge3"
shell.prefix(
    "source /usr/share/modules/init/bash 2>/dev/null || true && "
    # pandoc is needed by every rmarkdown::render rule and is absent on compute
    # nodes (only login has /usr/bin/pandoc) — load the module so it is on PATH
    # for all SLURM jobs (.Rprofile then finds it). Guarded so a node without
    # the module doesn't abort the shell.
    "module load pandoc/3.6.2 2>/dev/null || true && "
    # Load R 4.2.1 for ALL rules here (bash, where `module` is defined). Nodes
    # differ in default R (compms = 3.6.3, others = 4.2.1); the renv library is
    # built for 4.2.1, so a 3.6.3 node sees "no packages installed". Doing it in
    # the prefix (not `sh -c 'module load …'`, where `module` is undefined and
    # silently no-ops — the real cause of the compms WGBS failures) guarantees
    # 4.2.1 + a complete renv on every node for bare-Rscript rules.
    "module load r/4.2.1 2>/dev/null || true && "
    # Cap multithreaded OpenBLAS/OMP to the SLURM allocation (read at BLAS load,
    # before R). Unset, openblas-pthread uses the whole node's cores, ignoring
    # the cgroup → oversubscription on matrix-heavy stages (vst, rowSds, joins).
    "export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK && "
    f"source {CONDA_BASE}/etc/profile.d/conda.sh && "
    f"source {CONDA_BASE}/etc/profile.d/mamba.sh && "
)

localrules: all, clean

# ── Wildcard values ────────────────────────────────────────────────────────────
# ALL_* = universe of valid filters (wildcard constraints → any is buildable by
# explicit target). TRANSCRIPT_FILTERS = default build set (rule all / ML space).
ALL_TRANSCRIPT_FILTERS = config.get("all_transcript_filters", config["transcript_filters"])
TRANSCRIPT_FILTERS = config["transcript_filters"]
EVENT_TYPES        = config["event_types"]
VARIABILITIES      = config["variabilities"]
GROUP_COLS         = config["group_cols"]
PRIMARY            = config["primary_filter"]

# Pangolin device toggle (no file edit): `--config pangolin_use_gpu=true` for GPU,
# default CPU. Derive the pangolin_scores gpu flag + threads so the existing
# _partition()/_gres()/_qos()/R() helpers route it correctly. CPU -> 16 threads
# (encoding parallelises over cores); GPU -> 4.
_pangolin_gpu = str(config.get("pangolin_use_gpu", True)).lower() in ("true", "1", "yes")
config["resources"]["pangolin_scores"]["gpu"] = 1 if _pangolin_gpu else 0
config["resources"]["pangolin_scores"]["threads"] = 4 if _pangolin_gpu else 16
# Pangolin GPU gres/env. Default = idle gpu01 titans via the cu121 pangolin-titan
# env (Titan V sm_70 verified). a40 path uses ihec-as (cu128). See config.
PANGOLIN_GRES = config.get("pangolin_gpu_gres", "gpu:a40:1") if _pangolin_gpu else None
_pangolin_titan = bool(_pangolin_gpu) and "titan" in (PANGOLIN_GRES or "")
PANGOLIN_ENV = config.get(
    "pangolin_env", "pangolin-titan" if _pangolin_titan else "ihec-as"
)
# batch: 4096 on the 48G a40; 1024 on the 12G titan or on CPU (host-RAM fit).
PANGOLIN_BATCH = 4096 if (_pangolin_gpu and not _pangolin_titan) else 1024

# W&B tracking toggle (no file edit): `--config wandb_enabled=true`. Off by
# default -- matches run_splicing_ml.py's own CLI default (--wandb not passed
# unless asked). Old scripts/slurm_ml_one_config.sh always ran `wandb login`
# and passed --wandb/--wandb-project when a project name was given; mirrored
# here as a single flat toggle rather than scripts/'s 5th positional CLI arg.
WANDB_ENABLED = str(config.get("wandb_enabled", False)).lower() in ("true", "1", "yes")
WANDB_PROJECT = config.get("wandb_project", "splicing-ml")
WANDB_NO_REQUIRE_AUTH = str(config.get("wandb_no_require_auth", False)).lower() in (
    "true", "1", "yes",
)


def wandb_login_cmd():
    # Loads stored credentials (e.g. ~/.netrc from a prior `wandb login`) into
    # the job environment, same as the old script did unconditionally when a
    # project name was passed.
    return "wandb login >/dev/null 2>&1 && " if WANDB_ENABLED else ""


def wandb_cli_args():
    if not WANDB_ENABLED:
        return ""
    args = f"--wandb --wandb-project {WANDB_PROJECT}"
    if WANDB_NO_REQUIRE_AUTH:
        args += " --wandb-no-require-auth"
    return args

# ML configs: all combos minus excluded ones
EXCLUDE_ML = {
    (e["transcript_filter"], e["variability"])
    for e in config.get("exclude_ml_configs", [])
}

def _ml_configs():
    for tf, et, var, gc in product(TRANSCRIPT_FILTERS, EVENT_TYPES, VARIABILITIES, GROUP_COLS):
        if (tf, var) not in EXCLUDE_ML:
            yield et, tf, var, gc

ML_CONFIGS = list(_ml_configs())

# Feature-group ablation battery (§4.15) folded into `rule all` -- PRIMARY
# filter only, seqnames only (ontology excluded: leaks event identity on
# low-variability subsets, see memory project_splicing_ml_ontology_cv_event_leak.md).
# Sequence-only vs epigenetics-only across both event types x all 3 variability
# strata = 12 targets, matching the battery the user requested 2026-07-14.
ABLATION_FEATURE_GROUPS = ["sequence", "histone+dnam"]
ABLATION_CONFIGS = [
    (fg, et, var) for fg in ABLATION_FEATURE_GROUPS
    for et in EVENT_TYPES
    for var in VARIABILITIES
]

wildcard_constraints:
    transcript_filter = "|".join(ALL_TRANSCRIPT_FILTERS),
    event_type        = "|".join(EVENT_TYPES),
    variability       = "|".join(VARIABILITIES),
    group_col         = "|".join(GROUP_COLS),


# ── Helper: resource lookup ────────────────────────────────────────────────────
def R(rule_key, field):
    return config["resources"][rule_key][field]

def ml_resource(et, var, field):
    if et == "SE" and var == "both":
        tier = "splicing_ml_L"
    elif et == "RI" and var in ("High", "Low"):
        tier = "splicing_ml_S"
    else:
        tier = "splicing_ml_M"
    return config["resources"][tier][field]

def _partition(rule_key):
    return (config["slurm_gpu_partition"] if config["resources"][rule_key]["gpu"]
            else config["slurm_cpu_partition"])

def _extra(rule_key):
    base = "--mail-type=FAIL --mail-user=quirin.manz@tum.de"
    if config["resources"][rule_key]["gpu"]:
        return f"--exclude={config['slurm_gpu_exclude']} {base}"
    return base

def _qos(rule_key):
    return config["slurm_gpu_qos"] if config["resources"][rule_key]["gpu"] else None

def _gres(rule_key):
    return config["slurm_gpu_gres"] if config["resources"][rule_key]["gpu"] else None

def _ml_partition(wc):
    gpu = ml_resource(wc.event_type, wc.variability, "gpu")
    return config["slurm_gpu_partition"] if gpu else config["slurm_cpu_partition"]

def _ml_extra(wc):
    base = "--mail-type=FAIL --mail-user=quirin.manz@tum.de"
    if ml_resource(wc.event_type, wc.variability, "gpu"):
        return f"--exclude={config['slurm_gpu_exclude']} {base}"
    return base

def _ml_qos(wc):
    gpu = ml_resource(wc.event_type, wc.variability, "gpu")
    return config["slurm_gpu_qos"] if gpu else None

def _ml_gres(wc):
    gpu = ml_resource(wc.event_type, wc.variability, "gpu")
    return config["slurm_gpu_gres"] if gpu else None


# ── ChIP BigWig file list (for aggregate_chip_one wildcard) ───────────────────
# ChIP tracks come in TWO namings: observed (ihec.chipseq…{uuid}.pval.signal.bigwig,
# ~2295) and ChromImpute imputed (impute_{epirr}_{mark}.pval.bw, ~135). Two globs
# → both aggregated (05 §4.13l flags the observed/imputed source per mark). The
# single observed glob previously missed all 135 imputed tracks.
CHIP_BW,  = glob_wildcards("/nfs/data3/IHEC/ChIP-Seq/{bw}.pval.signal.bigwig")
CHIP_IMP, = glob_wildcards("/nfs/data3/IHEC/ChIP-Seq/impute_{bw}.pval.bw")


# ── Rule all ──────────────────────────────────────────────────────────────────
rule all:
    input:
        # Per-filter aggregated data (default: primary filter only — request
        # other filters' targets or extend config["transcript_filters"] to
        # build them; each is an independent, parallel DAG branch).
        f"processed_data/aggregated_dt_filtered_{PRIMARY}.csv.gz",
        f"processed_data/correlation_intrinsic_{PRIMARY}.csv.gz",
        # ML global models (all configs) -- classification/regression run as
        # separate parallel jobs (see rule splicing_ml_classification/_regression)
        [f"splicing_ml/output/{et}_{tf}_{var}_{gc}/splicing_ml_results_{task}.pkl.gz"
         for et, tf, var, gc in ML_CONFIGS for task in ("classification", "regression")],
        # Feature-group ablation battery (PRIMARY filter, seqnames only)
        [f"splicing_ml/output_ablation/{fg}/{et}_{PRIMARY}_{var}_seqnames/splicing_ml_results_{task}.pkl.gz"
         for fg, et, var in ABLATION_CONFIGS for task in ("classification", "regression")],
        # Fig2B-style comparison plot over the two batteries above (rule ml_global_comparison)
        f"reports/07-2-ml-global-comparison_{PRIMARY}.html",
        # Event-specific models (primary filter only — extend if needed):
        # Tier-1 ridge-screen results + Tier-2 elastic-net on hits.
        f"processed_data/event_models/{PRIMARY}/screen_results.csv.gz",
        f"processed_data/event_models/{PRIMARY}/.done",
        # QC / diagnostic notebooks (previously standalone-only; now part of
        # the default build per user direction 2026-07-14: "all should be all")
        f"reports/04-2-chip-signal-sanity_{PRIMARY}.html",
        f"qc/irfinder_tree_comparison_{PRIMARY}.csv.gz",
        f"reports/05b-feature-pca-sanity_{PRIMARY}.html",
        # Final analysis reports
        f"reports/09-2-ml-local-new_{PRIMARY}.html",
        f"reports/10-experimental-events_{PRIMARY}.html",


# ── Step 01: gather data ──────────────────────────────────────────────────────
rule gather_data:
    input:
        rmd = "01-gather-data.Rmd",
    output:
        file_table      = "processed_data/file_table.csv.gz",
        qc_summary      = "processed_data/qc_summary.csv",
        qc_covariates   = "processed_data/qc_flag_covariates.csv",
    log: "logs/01_gather_data.log"
    threads: R("gather_data", "threads")
    resources:
        mem_mb          = R("gather_data", "mem_mb"),
        runtime         = R("gather_data", "runtime"),
        slurm_partition = _partition("gather_data"),
        slurm_extra     = _extra("gather_data"),
        qos             = _qos("gather_data"),
        gres            = _gres("gather_data"),
    shell:
        """
        Rscript -e "rmarkdown::render('01-gather-data.Rmd',
            output_file = normalizePath('reports/01-gather-data.html', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 02-1: transcript filters + isoform quants + GTF annotation ───────────
# Cohort-wide (runs once). SUPPA2 reference GTFs + filtered transcript ids +
# the rna_to_use isoform-quant subset consumed by 02-2.
rule transcript_filters:
    input:
        rmd = "02-1-transcript-filters.Rmd",
        file_table = "processed_data/file_table.csv.gz",
    output:
        filtered_tx_ids = "splicing_analysis/filtered_transcript_ids.rds",
        tpm_expr        = "splicing_analysis/suppa/tpm_expressions.tsv.gz",
        # isoform_quants subset is consumed by rnaseq_normalisation (02-2) — declare
        # it so the DAG can produce it (was an undeclared side-output before).
        isoform_quants  = "splicing_analysis/isoform_quantifications_subset.tsv.gz",
        gencode_gtfs    = expand(
            "splicing_analysis/gencode.v29.{tf}.gtf",
            tf=ALL_TRANSCRIPT_FILTERS,
        ),
    log: "logs/02-1_transcript_filters.log"
    threads: R("transcript_filters", "threads")
    resources:
        mem_mb          = R("transcript_filters", "mem_mb"),
        runtime         = R("transcript_filters", "runtime"),
        slurm_partition = _partition("transcript_filters"),
        slurm_extra     = _extra("transcript_filters"),
        qos             = _qos("transcript_filters"),
        gres            = _gres("transcript_filters"),
    shell:
        """
        Rscript -e "rmarkdown::render('02-1-transcript-filters.Rmd',
            output_file = normalizePath('reports/02-1-transcript-filters.html', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 02-2: RNA-Seq gene expression normalisation (PLAN §4.12b) ────────────
# Per-transcript_filter. Runs BEFORE 02-3 (its VST output is the input to the
# 02-3 event-expression gate) — depends only on 02-1 outputs + file_table, so
# it sits at stage 02 in the reorder. GeTMM-CPM + DESeq2 vst side by side
# (routing: getmm→splicing_ml, vst→07/09/RBP).
rule rnaseq_normalisation:
    input:
        rmd = "02-2-rnaseq-normalisation.Rmd",
        filtered_tx_ids = "splicing_analysis/filtered_transcript_ids.rds",
        file_table      = "processed_data/file_table.csv.gz",
        isoform_quants  = "splicing_analysis/isoform_quantifications_subset.tsv.gz",
    output:
        gene_expression_normalised =
            "processed_data/gene_expression_normalised_{transcript_filter}.csv.gz",
        html = "reports/02-2-rnaseq-normalisation_{transcript_filter}.html",
    log: "logs/02-2_rnaseq_normalisation_{transcript_filter}.log"
    threads: R("rnaseq_normalisation", "threads")
    resources:
        mem_mb          = R("rnaseq_normalisation", "mem_mb"),
        runtime         = R("rnaseq_normalisation", "runtime"),
        slurm_partition = _partition("rnaseq_normalisation"),
        slurm_extra     = _extra("rnaseq_normalisation"),
        qos             = _qos("rnaseq_normalisation"),
        gres            = _gres("rnaseq_normalisation"),
    shell:
        """
        TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        Rscript -e "rmarkdown::render('02-2-rnaseq-normalisation.Rmd',
            output_file = normalizePath('{output.html}', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 02-3: rMATS PSI + Procedure-2 + VST gene-expression gate ─────────────
# Per-transcript_filter (PLAN §4.12e). Depends on 02-2's per-filter normalised
# expression (the VST gate) — this is the edge that sequences normalisation
# before event filtering.
rule rmats_event_filtering:
    input:
        rmd = "02-3-rmats-event-filtering.Rmd",
        file_table = "processed_data/file_table.csv.gz",
        gene_expression_normalised =
            "processed_data/gene_expression_normalised_{transcript_filter}.csv.gz",
    output:
        psi_files = expand(
            "splicing_analysis/rmats/{{transcript_filter}}/event_{et}.psi",
            et=EVENT_TYPES,
        ),
        jc_files = expand(
            "splicing_analysis/rmats/{{transcript_filter}}/event_{et}.jc.csv.gz",
            et=EVENT_TYPES,
        ),
        vst_cutoff = "qc/vst_expression_cutoff_{transcript_filter}.csv",
    log: "logs/02-3_rmats_event_filtering_{transcript_filter}.log"
    threads: R("rmats_event_filtering", "threads")
    resources:
        mem_mb          = R("rmats_event_filtering", "mem_mb"),
        runtime         = R("rmats_event_filtering", "runtime"),
        slurm_partition = _partition("rmats_event_filtering"),
        slurm_extra     = _extra("rmats_event_filtering"),
        qos             = _qos("rmats_event_filtering"),
        gres            = _gres("rmats_event_filtering"),
    shell:
        """
        TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        Rscript -e "rmarkdown::render('02-3-rmats-event-filtering.Rmd',
            output_file = normalizePath('reports/02-3-rmats-event-filtering_{wildcards.transcript_filter}.html', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 02-4: IRFinder vs rMATS RI concordance (diagnostic only) ─────────────
# Standalone target, NOT in rule all / not on 03's critical path (03 keeps
# reading 02-3's RI .psi unchanged) - see revision/file-changes/
# 02-4-irfinder-concordance.Rmd.md. IRFinder data itself
# (/nfs/data/IHEC/RNAseq/irfinder/) is external, pre-computed, not a Snakemake
# input (thousands of per-sample files, static reference data - not tracked).
rule irfinder_concordance:
    input:
        rmd = "02-4-irfinder-concordance.Rmd",
        rmats_ri_psi = "splicing_analysis/rmats/{transcript_filter}/event_RI.psi",
        rmats_ri_jc  = "splicing_analysis/rmats/{transcript_filter}/event_RI.jc.csv.gz",
    output:
        # Runs BOTH IRFinder reference trees (annotation, TSL12) for real
        # rather than assuming one - see file-changes/02-4-irfinder-concordance.Rmd.md.
        concordance = expand(
            "qc/irfinder_concordance_{tree}_{{transcript_filter}}.csv.gz",
            tree=["annotation", "TSL12"],
        ),
        tree_comparison = "qc/irfinder_tree_comparison_{transcript_filter}.csv.gz",
        stratified      = "qc/irfinder_stratified_metrics_{transcript_filter}.csv.gz",
        # reverse direction: introns IRFinder finds with recurring signal that
        # have no rMATS RI candidate at all (see file-changes doc "reverse
        # direction" section).
        irfinder_only = expand(
            "qc/irfinder_only_introns_{tree}_{{transcript_filter}}.csv.gz",
            tree=["annotation", "TSL12"],
        ),
        hist_pdfs = expand(
            "qc/distributions/irfinder_concordance_{tree}_{{transcript_filter}}.pdf",
            tree=["annotation", "TSL12"],
        ),
        scatter_pdfs = expand(
            "qc/distributions/irfinder_scatter_{tree}_{{transcript_filter}}.pdf",
            tree=["annotation", "TSL12"],
        ),
        stratified_pdfs = expand(
            "qc/distributions/irfinder_stratified_{metric}_{{transcript_filter}}.pdf",
            metric=["mean_coverage", "intron_length", "mean_psi", "sd_psi", "n_samples"],
        ),
        html = "reports/02-4-irfinder-concordance_{transcript_filter}.html",
    log: "logs/02-4_irfinder_concordance_{transcript_filter}.log"
    threads: R("irfinder_concordance", "threads")
    resources:
        mem_mb          = R("irfinder_concordance", "mem_mb"),
        runtime         = R("irfinder_concordance", "runtime"),
        slurm_partition = _partition("irfinder_concordance"),
        slurm_extra     = _extra("irfinder_concordance"),
        qos             = _qos("irfinder_concordance"),
        gres            = _gres("irfinder_concordance"),
    shell:
        """
        TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        Rscript -e "rmarkdown::render('02-4-irfinder-concordance.Rmd',
            output_file = normalizePath('{output.html}', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 03: prepare aggregation workspace ────────────────────────────────────
rule prepare_aggregation:
    input:
        rmd = "03-prepare-aggregation.Rmd",
        psi_files = expand(
            "splicing_analysis/rmats/{{transcript_filter}}/event_{et}.psi",
            et=EVENT_TYPES,
        ),
        jc_files = expand(
            "splicing_analysis/rmats/{{transcript_filter}}/event_{et}.jc.csv.gz",
            et=EVENT_TYPES,
        ),
        filtered_tx_ids = "splicing_analysis/filtered_transcript_ids.rds",
        gencode_gtf = "splicing_analysis/gencode.v29.{transcript_filter}.gtf",
    output:
        aggregateOver_bed     = "processed_data/aggregateOver_{transcript_filter}.bed",
        keep_rows_manual      = "processed_data/keep_rows_manual_{transcript_filter}.rds",
        ijc_sjc_dt            = "processed_data/ijc_sjc_dt_{transcript_filter}.csv.gz",
        event_annotations_dt  = "processed_data/event_annotations_dt_{transcript_filter}.csv.gz",
        psi_long_dt           = "processed_data/psi_long_dt_{transcript_filter}.csv.gz",
        pangolin_events       = "processed_data/pangolin_events_{transcript_filter}.csv",
        ss5_fasta             = "processed_data/5ss_{transcript_filter}.fasta",
        ss5up_fasta           = "processed_data/5ss_up_{transcript_filter}.fasta",
        ss3_fasta             = "processed_data/3ss_{transcript_filter}.fasta",
        ss3down_fasta         = "processed_data/3ss_down_{transcript_filter}.fasta",
        sample_cols           = "processed_data/sample_cols_{transcript_filter}.rds",
        # 09-1's chromHMM-vicinity per-event feature engineering (a separate
        # local-modelling path from splicing_ml's pooled genome-wide route) —
        # these 3 used to leave 03's session only via the removed aggregating.rda
        # workspace dump; now persisted directly. Must come from the SAME run as
        # keep_rows_manual (see 03's saveRDS comment: internal index consistency).
        event_gr              = "processed_data/event_gr_{transcript_filter}.rds",
        active_chromhmm       = "processed_data/activeChromHMM_{transcript_filter}.rds",
        chromhmm_hits         = "processed_data/chromhmm_hits_{transcript_filter}.rds",
    log: "logs/03_prepare_aggregation_{transcript_filter}.log"
    threads: R("prepare_aggregation", "threads")
    resources:
        mem_mb          = R("prepare_aggregation", "mem_mb"),
        runtime         = R("prepare_aggregation", "runtime"),
        slurm_partition = _partition("prepare_aggregation"),
        slurm_extra     = _extra("prepare_aggregation"),
        qos             = _qos("prepare_aggregation"),
        gres            = _gres("prepare_aggregation"),
    shell:
        """
        TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        Rscript -e "rmarkdown::render('03-prepare-aggregation.Rmd',
            output_file = normalizePath('reports/03-prepare-aggregation_{wildcards.transcript_filter}.html', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 04-1: aggregate WGBS ─────────────────────────────────────────────────
rule aggregate_wgbs:
    input:
        script = "04-1-aggregate-WGBS-matrix.R",
        aggregateOver = "processed_data/aggregateOver_{transcript_filter}.bed",
    output:
        wgbs = "sample_dts/WGBS_agg_{transcript_filter}.csv.gz",
    log: "logs/04-1_aggregate_wgbs_{transcript_filter}.log"
    threads: R("aggregate_wgbs", "threads")
    resources:
        mem_mb          = R("aggregate_wgbs", "mem_mb"),
        runtime         = R("aggregate_wgbs", "runtime"),
        slurm_partition = _partition("aggregate_wgbs"),
        slurm_extra     = _extra("aggregate_wgbs"),
        qos             = _qos("aggregate_wgbs"),
        gres            = _gres("aggregate_wgbs"),
    shell:
        """
        TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        Rscript 04-1-aggregate-WGBS-matrix.R > {log} 2>&1
        """


# ── Step 04-2: aggregate ChIP-Seq signal ──────────────────────────────────────
# One aggregate_chip_one job per bigwig (~2295). Submitted as INDIVIDUAL SLURM
# jobs, throttled by the profile's `jobs` (concurrency) + `max-jobs-per-second`
# (submit rate). The SLURM-array grouping (`slurm-array-jobs`) was tried but
# stalled the scheduler on this cluster and is disabled — see
# profiles/slurm/config.yaml. Per-bigwig outputs stay independent, so a rerun
# only rebuilds missing tabs. CHIP_BW is discovered at parse time via
# glob_wildcards. aggregate_chip is a sentinel that fans the per-bigwig tabs
# into a single .done for 05.
rule aggregate_chip_one:
    input:
        bw  = "/nfs/data3/IHEC/ChIP-Seq/{bw}.pval.signal.bigwig",
        bed = "processed_data/aggregateOver_{transcript_filter}.bed",
    output:
        tab = "sample_dts/chip_agg_{transcript_filter}/{bw}.pval.signal.bigwig.tab.gz",
    log: "logs/chip/{transcript_filter}/{bw}.log"
    threads: R("aggregate_chip_one", "threads")
    resources:
        mem_mb          = R("aggregate_chip_one", "mem_mb"),
        runtime         = R("aggregate_chip_one", "runtime"),
        slurm_partition = _partition("aggregate_chip_one"),
        slurm_extra     = "--mail-type=NONE",
    shell:
        """
        set -euo pipefail
        out=sample_dts/chip_agg_{wildcards.transcript_filter}/{wildcards.bw}.pval.signal.bigwig.tab
        mkdir -p "$(dirname "$out")"
        # mamba run: shell.prefix sources conda/mamba but does NOT activate an
        # env, so a bare `bigWigAverageOverBed` on a compute node resolves to a
        # broken/incompatible binary (dies after "processing chromosomes"). Force
        # the ihec-as ucsc-bigwigaverageoverbed.
        # set -e above: if bigWig fails (e.g. OOM-Killed) the script aborts BEFORE
        # gzip, so no partial/empty .tab.gz is left at the output path.
        mamba run -n ihec-as bigWigAverageOverBed {input.bw} {input.bed} "$out" -minMax > {log} 2>&1
        gzip -f "$out"
        """

# Imputed ChIP tracks (impute_{bw}.pval.bw) — same aggregation, different naming.
# Output basename matches file_table's basename(file_path) so 05 finds them.
rule aggregate_chip_imp:
    input:
        bw  = "/nfs/data3/IHEC/ChIP-Seq/impute_{bw}.pval.bw",
        bed = "processed_data/aggregateOver_{transcript_filter}.bed",
    output:
        tab = "sample_dts/chip_agg_{transcript_filter}/impute_{bw}.pval.bw.tab.gz",
    log: "logs/chip/{transcript_filter}/impute_{bw}.log"
    threads: R("aggregate_chip_one", "threads")
    resources:
        mem_mb          = R("aggregate_chip_one", "mem_mb"),
        runtime         = R("aggregate_chip_one", "runtime"),
        slurm_partition = _partition("aggregate_chip_one"),
        slurm_extra     = "--mail-type=NONE",
    shell:
        """
        set -euo pipefail
        out=sample_dts/chip_agg_{wildcards.transcript_filter}/impute_{wildcards.bw}.pval.bw.tab
        mkdir -p "$(dirname "$out")"
        mamba run -n ihec-as bigWigAverageOverBed {input.bw} {input.bed} "$out" -minMax > {log} 2>&1
        gzip -f "$out"
        """

rule aggregate_chip:
    input:
        expand(
            "sample_dts/chip_agg_{{transcript_filter}}/{bw}.pval.signal.bigwig.tab.gz",
            bw=CHIP_BW,
        ),
        expand(
            "sample_dts/chip_agg_{{transcript_filter}}/impute_{bw}.pval.bw.tab.gz",
            bw=CHIP_IMP,
        ),
    output:
        done = touch("sample_dts/chip_agg_{transcript_filter}.done"),
    localrule: True

# §4.13b ChIP signal sanity — QC check, NOT an input to create_aggregated_dt,
# so it never gates the data build (in `rule all` since 2026-07-14, but as a
# leaf, not a dependency). Needs the whole per-filter ChIP aggregation together
# (streams each tab once), hence a real compute job downstream of the
# aggregate_chip sentinel — not folded into that localrule.
rule chip_signal_sanity:
    input:
        rmd = "04-2-chip-signal-sanity.Rmd",
        chip_done = "sample_dts/chip_agg_{transcript_filter}.done",
        bed = "processed_data/aggregateOver_{transcript_filter}.bed",
        event_annotations_dt = "processed_data/event_annotations_dt_{transcript_filter}.csv.gz",
        file_table = "processed_data/file_table.csv.gz",
    output:
        html = "reports/04-2-chip-signal-sanity_{transcript_filter}.html",
        per_file_summary = "qc/chip_per_file_summary_{transcript_filter}.csv.gz",
    log: "logs/04-2_chip_signal_sanity_{transcript_filter}.log"
    threads: R("chip_signal_sanity", "threads")
    resources:
        mem_mb          = R("chip_signal_sanity", "mem_mb"),
        runtime         = R("chip_signal_sanity", "runtime"),
        slurm_partition = _partition("chip_signal_sanity"),
        slurm_extra     = _extra("chip_signal_sanity"),
        qos             = _qos("chip_signal_sanity"),
        gres            = _gres("chip_signal_sanity"),
    shell:
        """
        TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        Rscript -e "rmarkdown::render('04-2-chip-signal-sanity.Rmd',
            output_file = normalizePath('{output.html}', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 04-3: MaxEntScan splice-site scores ──────────────────────────────────
rule maxentscan_scores:
    input:
        script   = "04-3-maxentscan-scores.sh",
        ss3      = "processed_data/3ss_{transcript_filter}.fasta",
        ss3down  = "processed_data/3ss_down_{transcript_filter}.fasta",
        ss5      = "processed_data/5ss_{transcript_filter}.fasta",
        ss5up    = "processed_data/5ss_up_{transcript_filter}.fasta",
    output:
        score3     = "processed_data/3scores_{transcript_filter}.txt",
        score3down = "processed_data/3down_scores_{transcript_filter}.txt",
        score5     = "processed_data/5scores_{transcript_filter}.txt",
        score5up   = "processed_data/5up_scores_{transcript_filter}.txt",
    log: "logs/04-3_maxentscan_scores_{transcript_filter}.log"
    threads: R("maxentscan_scores", "threads")
    resources:
        mem_mb          = R("maxentscan_scores", "mem_mb"),
        runtime         = R("maxentscan_scores", "runtime"),
        slurm_partition = _partition("maxentscan_scores"),
        slurm_extra     = _extra("maxentscan_scores"),
        qos             = _qos("maxentscan_scores"),
        gres            = _gres("maxentscan_scores"),
    shell:
        """
        TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        mamba run -n ihec-as bash 04-3-maxentscan-scores.sh > {log} 2>&1
        """


# ── Step 04-4: Pangolin splice-site scores ────────────────────────────────────
rule pangolin_scores:
    input:
        script     = "04-4-pangolin-scores.sh",
        events_csv = "processed_data/pangolin_events_{transcript_filter}.csv",
        fasta      = "data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
    output:
        scores_csv = "processed_data/pangolin_scores_{transcript_filter}.csv",
    log: "logs/04-4_pangolin_scores_{transcript_filter}.log"
    params:
        batch = PANGOLIN_BATCH,   # 4096 a40 / 1024 titan|CPU
        env   = PANGOLIN_ENV,     # pangolin-titan (cu121, titan) or ihec-as (cu128, a40)
    threads: R("pangolin_scores", "threads")
    resources:
        mem_mb          = R("pangolin_scores", "mem_mb"),
        runtime         = R("pangolin_scores", "runtime"),
        slurm_partition = _partition("pangolin_scores"),
        slurm_extra     = _extra("pangolin_scores"),
        qos             = _qos("pangolin_scores"),
        gres            = PANGOLIN_GRES,   # titan (gpu01), not the a40 default
    shell:
        """
        PANGOLIN_ENV={params.env} TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        bash 04-4-pangolin-scores.sh --batch-size {params.batch} > {log} 2>&1
        """


# ── Step 04-6: RBP binding-site annotation (PLAN §3.4) ────────────────────────
# Step 1 (binding-site annotation) AND Step 2 (RBP-score + core-spliceosome-factor
# expression features, on gene_expression_vst) are BOTH implemented. This rule
# takes `gene_expression_normalised_{tf}.csv.gz` (from rnaseq_normalisation) +
# `event_annotations_dt_{tf}` (from 03) as inputs and writes the Step 2 feature
# tables (`rbp_score_dt`, `splicing_factor_expression`, `rbp_wide_expression`),
# which create_aggregated_dt consumes (see its input block) — so the DAG chain is
# rnaseq_normalisation + prepare_aggregation -> rbp_binding_sites -> create_aggregated_dt.
# `gencode_gtf` is a static reference file (like `human_postar3`), not a
# pipeline-generated artifact — no upstream rule produces it.
rule rbp_binding_sites:
    input:
        rmd = "04-5-rbp-binding-sites.Rmd",
        event_annotations_dt = "processed_data/event_annotations_dt_{transcript_filter}.csv.gz",
        # Step 2 turns binding sites into per-sample expression features -> needs
        # the per-filter normalised expression (uses gene_expression_vst).
        gene_expr             = "processed_data/gene_expression_normalised_{transcript_filter}.csv.gz",
        human_postar3         = "data/human.txt.gz",
        gencode_gtf            = "splicing_analysis/gencode.v29.annotation.gtf",
    output:
        rbp_per_event    = "processed_data/rbp_per_event_{transcript_filter}.rds",
        rbp_gene_ids     = "processed_data/rbp_gene_ids_{transcript_filter}.rds",
        # Step 2 feature tables (all on gene_expression_vst):
        rbp_score        = "processed_data/rbp_score_dt_{transcript_filter}.csv.gz",
        spliceosome_expr = "processed_data/splicing_factor_expression_{transcript_filter}.csv.gz",
        rbp_wide         = "processed_data/rbp_wide_expression_{transcript_filter}.csv.gz",
        spliceosome_ids  = "processed_data/spliceosome_gene_ids_{transcript_filter}.rds",
        html             = "reports/04-5-rbp-binding-sites_{transcript_filter}.html",
    log: "logs/04-5_rbp_binding_sites_{transcript_filter}.log"
    threads: R("rbp_binding_sites", "threads")
    resources:
        mem_mb          = R("rbp_binding_sites", "mem_mb"),
        runtime         = R("rbp_binding_sites", "runtime"),
        slurm_partition = _partition("rbp_binding_sites"),
        slurm_extra     = _extra("rbp_binding_sites"),
        qos             = _qos("rbp_binding_sites"),
        gres            = _gres("rbp_binding_sites"),
    shell:
        """
        TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        Rscript -e "rmarkdown::render('04-5-rbp-binding-sites.Rmd',
            output_file = normalizePath('{output.html}', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 05: create aggregated dataset (per transcript_filter) ────────────────
rule create_aggregated_dt:
    input:
        rmd                   = "05-create-aggregated-dt.Rmd",
        keep_rows_manual      = "processed_data/keep_rows_manual_{transcript_filter}.rds",
        sample_cols           = "processed_data/sample_cols_{transcript_filter}.rds",
        file_table            = "processed_data/file_table.csv.gz",
        ijc_sjc_dt            = "processed_data/ijc_sjc_dt_{transcript_filter}.csv.gz",
        event_annotations_dt  = "processed_data/event_annotations_dt_{transcript_filter}.csv.gz",
        psi_long_dt           = "processed_data/psi_long_dt_{transcript_filter}.csv.gz",
        gene_expr             = "processed_data/gene_expression_normalised_{transcript_filter}.csv.gz",
        vst_cutoff            = "qc/vst_expression_cutoff_{transcript_filter}.csv",
        wgbs                  = "sample_dts/WGBS_agg_{transcript_filter}.csv.gz",
        chip_done             = "sample_dts/chip_agg_{transcript_filter}.done",
        qc_covariates         = "processed_data/qc_flag_covariates.csv",
        filtered_tx_ids       = "splicing_analysis/filtered_transcript_ids.rds",
        maxentscan_score3     = "processed_data/3scores_{transcript_filter}.txt",
        maxentscan_score3down = "processed_data/3down_scores_{transcript_filter}.txt",
        maxentscan_score5     = "processed_data/5scores_{transcript_filter}.txt",
        maxentscan_score5up   = "processed_data/5up_scores_{transcript_filter}.txt",
        ss3_fasta             = "processed_data/3ss_{transcript_filter}.fasta",
        ss5_fasta             = "processed_data/5ss_{transcript_filter}.fasta",
        pangolin_scores       = "processed_data/pangolin_scores_{transcript_filter}.csv",
        # RBP + core-spliceosome expression features (04-6 Step 2, §3.4/§4.12f)
        rbp_score             = "processed_data/rbp_score_dt_{transcript_filter}.csv.gz",
        spliceosome_expr      = "processed_data/splicing_factor_expression_{transcript_filter}.csv.gz",
    output:
        csv  = "processed_data/aggregated_dt_filtered_{transcript_filter}.csv.gz",
        html = "reports/05-create-aggregated-dt_{transcript_filter}.html",
    log: "logs/05_create_aggregated_dt_{transcript_filter}.log"
    threads: R("create_aggregated_dt", "threads")
    resources:
        mem_mb          = R("create_aggregated_dt", "mem_mb"),
        runtime         = R("create_aggregated_dt", "runtime"),
        slurm_partition = _partition("create_aggregated_dt"),
        slurm_extra     = _extra("create_aggregated_dt"),
        qos             = _qos("create_aggregated_dt"),
        gres            = _gres("create_aggregated_dt"),
    shell:
        """
        TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        Rscript -e "rmarkdown::render('05-create-aggregated-dt.Rmd',
            output_file = normalizePath('{output.html}', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 05b: feature matrix PCA/UMAP sanity + embedding gallery ──────────────
# Fixed to PRIMARY (not wildcarded), same pattern as correlation/ml_analysis —
# 05b is per-filter (TRANSCRIPT_FILTER env) but only the primary filter is
# wired into the DAG today. Per-(Event Type x Variability) plot files under
# plot_dir/{feature_umap_pca,psi_embedding}/ are runtime-discovered from the
# data (not statically enumerable) -- html report is the tracked output, same
# treatment as irfinder_concordance's report vs its per-tree/per-metric plots.
rule feature_pca_sanity:
    input:
        rmd            = "05b-feature-pca-sanity.Rmd",
        preprocessing  = "05zz-feature-preprocessing.R",
        aggregated_dt  = f"processed_data/aggregated_dt_filtered_{PRIMARY}.csv.gz",
    output:
        html = f"reports/05b-feature-pca-sanity_{PRIMARY}.html",
    log: f"logs/05b_feature_pca_sanity_{PRIMARY}.log"
    threads: R("feature_pca_sanity", "threads")
    resources:
        mem_mb          = R("feature_pca_sanity", "mem_mb"),
        runtime         = R("feature_pca_sanity", "runtime"),
        slurm_partition = _partition("feature_pca_sanity"),
        slurm_extra     = _extra("feature_pca_sanity"),
        qos             = _qos("feature_pca_sanity"),
        gres            = _gres("feature_pca_sanity"),
    shell:
        f"""
        TRANSCRIPT_FILTER={PRIMARY} \\
        Rscript -e "rmarkdown::render('05b-feature-pca-sanity.Rmd',
            output_file = normalizePath('{{output.html}}', mustWork = FALSE)
        )" > {{log}} 2>&1
        """


# ── Step 06: correlation analysis ─────────────────────────────────────────────
# Fixed to PRIMARY (not wildcarded), same pattern as ml_analysis — 06-correlation.Rmd
# is per-filter (TRANSCRIPT_FILTER env) but only the primary filter's correlation is
# wired into the DAG today.
rule correlation:
    input:
        rmd = "06-correlation.Rmd",
        # shared preprocessing helper source()d by 06 — declare so edits to it
        # trigger a correlation rerun (the render shell string is opaque to
        # Snakemake's code trigger). Same helper feeds 05b (not a rule; 05z deprecated).
        preprocessing = "05zz-feature-preprocessing.R",
        aggregated_dt = f"processed_data/aggregated_dt_filtered_{PRIMARY}.csv.gz",
        file_table    = "processed_data/file_table.csv.gz",
    output:
        corr_raw     = f"processed_data/correlation_intrinsic_{PRIMARY}.csv.gz",
        corr_preproc = f"processed_data/correlation_intrinsic_preproc_{PRIMARY}.csv.gz",
        html         = f"reports/06-correlation_{PRIMARY}.html",
    log: f"logs/06_correlation_{PRIMARY}.log"
    threads: R("correlation", "threads")
    resources:
        # attempt-scaled (2026-07-15): OOM-killed at flat mem_mb=48000 (5
        # oom_kill events, job 6329894) -- compute_all_cor's pbmclapply forks
        # up to 10 workers over the full 9M-row aggregated_dt, each fork's
        # copy-on-write footprint apparently exceeds the flat budget. Scale on
        # retry like the ML rules instead of guessing a single higher number.
        mem_mb          = lambda wc, attempt: R("correlation", "mem_mb") * attempt,
        runtime         = R("correlation", "runtime"),
        slurm_partition = _partition("correlation"),
        slurm_extra     = _extra("correlation"),
        qos             = _qos("correlation"),
        gres            = _gres("correlation"),
    shell:
        f"""
        TRANSCRIPT_FILTER={PRIMARY} \\
        Rscript -e "rmarkdown::render('06-correlation.Rmd',
            output_file = normalizePath('{{output.html}}', mustWork = FALSE)
        )" > {{log}} 2>&1
        """


# ── Step 07: splicing_ml Fig2B-style global comparison plot ───────────────────
# Reads splicing_ml_classification's + splicing_ml_ablation_classification's own
# JSON outputs (no aggregated_dt/large-data re-read) and recreates a comparison
# close to the preprint's Fig 2B: model=xgb balanced_accuracy/AUROC/AUPRC/MCC by
# held-out condition (Chr=seqnames, Cell=ontology -- no "Chr & Cell", that mode
# doesn't exist in this pipeline) x feature set (All/Non-Epigenetic/Epigenetic,
# from the ablation battery, seqnames only -- see ABLATION_CONFIGS above).
# Numbered 07 (the one unused slot between 06-correlation and 08-2's frozen old
# comparison) in anticipation of `splicing_ml/` eventually being renamed
# `07_splicing_ml/`. Explicit input lists (not a glob) so Snakemake's DAG
# correctly requires exactly the same 12+12 targets `rule all` already commits
# to for splicing_ml_classification/splicing_ml_ablation_classification --
# nothing extra is forced to build just for this plot.
_ML_GLOBAL_COMPARISON_MAIN_JSONS = [
    f"splicing_ml/output/{et}_{PRIMARY}_{var}_{gc}/splicing_ml_results_classification.json.gz"
    for et in EVENT_TYPES for var in VARIABILITIES for gc in GROUP_COLS
]
_ML_GLOBAL_COMPARISON_ABLATION_JSONS = [
    f"splicing_ml/output_ablation/{fg}/{et}_{PRIMARY}_{var}_seqnames/splicing_ml_results_classification.json.gz"
    for fg, et, var in ABLATION_CONFIGS
]

rule ml_global_comparison:
    input:
        rmd      = "07-2-ml-global-comparison.Rmd",
        wandb_py = "scripts/upload_comparison_to_wandb.py",
        main_jsons     = _ML_GLOBAL_COMPARISON_MAIN_JSONS,
        ablation_jsons = _ML_GLOBAL_COMPARISON_ABLATION_JSONS,
    output:
        html = f"reports/07-2-ml-global-comparison_{PRIMARY}.html",
        pdf  = "images/Rplots/07-2_ml_global_comparison.pdf",
        csv  = "images/Rplots/07-2_ml_global_comparison_data.csv",
    log: f"logs/07_2_ml_global_comparison_{PRIMARY}.log"
    threads: R("ml_global_comparison", "threads")
    resources:
        mem_mb          = R("ml_global_comparison", "mem_mb"),
        runtime         = R("ml_global_comparison", "runtime"),
        slurm_partition = _partition("ml_global_comparison"),
        slurm_extra     = _extra("ml_global_comparison"),
        qos             = _qos("ml_global_comparison"),
        gres            = _gres("ml_global_comparison"),
    params:
        wandb_login = wandb_login_cmd(),
    shell:
        f"""
        {{params.wandb_login}}Rscript -e "rmarkdown::render('07-2-ml-global-comparison.Rmd',
            output_file = normalizePath('{{output.html}}', mustWork = FALSE)
        )" > {{log}} 2>&1
        """


# ── Step 08: splicing_ml global models ────────────────────────────────────────
# Split into classification/regression rules (2026-07-14) so the two tasks run as
# independent, parallelizable SLURM jobs instead of serially inside one job --
# they use the same input data but train unrelated models, so there's no reason
# to pay for both sequentially in the same allocation. The only duplicated cost
# is the initial load_dataset() call (~seconds), shared inside a single
# run_splicing_ml.py invocation but now paid twice across the two jobs --
# negligible next to per-task training time (minutes to hours).
rule splicing_ml_classification:
    input:
        src  = SPLICING_ML_SRC,
        data = "processed_data/aggregated_dt_filtered_{transcript_filter}.csv.gz",
    output:
        pkl  = "splicing_ml/output/{event_type}_{transcript_filter}_{variability}_{group_col}/splicing_ml_results_classification.pkl.gz",
        json = "splicing_ml/output/{event_type}_{transcript_filter}_{variability}_{group_col}/splicing_ml_results_classification.json.gz",
    log: "logs/splicing_ml_classification_{event_type}_{transcript_filter}_{variability}_{group_col}.log"
    threads: lambda wc: ml_resource(wc.event_type, wc.variability, "threads")
    resources:
        mem_mb          = lambda wc, attempt: ml_resource(wc.event_type, wc.variability, "mem_mb") * attempt,
        runtime         = lambda wc: ml_resource(wc.event_type, wc.variability, "runtime"),
        slurm_partition = _ml_partition,
        slurm_extra     = _ml_extra,
        qos             = _ml_qos,
        gres            = _ml_gres,
    params:
        wandb_login = wandb_login_cmd(),
        wandb_args  = wandb_cli_args(),
    shell:
        """
        {params.wandb_login}PYTHONUNBUFFERED=1 mamba run --no-capture-output -n ihec-as python run_splicing_ml.py \
            --debug \
            --optuna \
            --data-path {input.data} \
            --output-dir splicing_ml/output/{wildcards.event_type}_{wildcards.transcript_filter}_{wildcards.variability}_{wildcards.group_col} \
            --only-event-type {wildcards.event_type} \
            --only-transcript-filter {wildcards.transcript_filter} \
            --only-variability {wildcards.variability} \
            --only-group {wildcards.group_col} \
            --skip-regression \
            --max-cores {threads} \
            {params.wandb_args} \
        > {log} 2>&1
        """


rule splicing_ml_regression:
    input:
        src  = SPLICING_ML_SRC,
        data = "processed_data/aggregated_dt_filtered_{transcript_filter}.csv.gz",
    output:
        pkl  = "splicing_ml/output/{event_type}_{transcript_filter}_{variability}_{group_col}/splicing_ml_results_regression.pkl.gz",
        json = "splicing_ml/output/{event_type}_{transcript_filter}_{variability}_{group_col}/splicing_ml_results_regression.json.gz",
    log: "logs/splicing_ml_regression_{event_type}_{transcript_filter}_{variability}_{group_col}.log"
    threads: lambda wc: ml_resource(wc.event_type, wc.variability, "threads")
    resources:
        mem_mb          = lambda wc, attempt: ml_resource(wc.event_type, wc.variability, "mem_mb") * attempt,
        runtime         = lambda wc: ml_resource(wc.event_type, wc.variability, "runtime"),
        slurm_partition = _ml_partition,
        slurm_extra     = _ml_extra,
        qos             = _ml_qos,
        gres            = _ml_gres,
    params:
        wandb_login = wandb_login_cmd(),
        wandb_args  = wandb_cli_args(),
    shell:
        """
        {params.wandb_login}PYTHONUNBUFFERED=1 mamba run --no-capture-output -n ihec-as python run_splicing_ml.py \
            --debug \
            --optuna \
            --data-path {input.data} \
            --output-dir splicing_ml/output/{wildcards.event_type}_{wildcards.transcript_filter}_{wildcards.variability}_{wildcards.group_col} \
            --only-event-type {wildcards.event_type} \
            --only-transcript-filter {wildcards.transcript_filter} \
            --only-variability {wildcards.variability} \
            --only-group {wildcards.group_col} \
            --skip-classification \
            --max-cores {threads} \
            {params.wandb_args} \
        > {log} 2>&1
        """


# ── Step 08b: splicing_ml feature-group ablations (§4.15) ─────────────────────
# In `rule all` since 2026-07-14 (ABLATION_CONFIGS: sequence + histone+dnam x
# both event types x all 3 variability strata, PRIMARY filter, seqnames only
# -- 12 targets, "all should be all" per user direction). The full 8-group
# FEATURE_GROUPS space (per-group-only: rbp/spliceosome/gene_expression/dnam/
# histone alone) is NOT in ABLATION_CONFIGS -- only the two combos actually
# requested so far. Extend ABLATION_FEATURE_GROUPS above to add more, or
# request other combos directly, e.g.:
#   snakemake --profile profiles/slurm \
#     splicing_ml/output_ablation/histone/RI_biotype_filtered_Low_seqnames/splicing_ml_results_classification.pkl.gz
# {feature_groups} is a "+"-joined FEATURE_GROUPS combo (matches
# wandb_tracker.py's own run-name convention for the same concept, e.g.
# "sequence" or "histone+dnam") -- translated back to `--feature-groups a b`
# (space-separated) for the CLI.
# group_col restricted to seqnames ONLY (see wildcard_constraints below):
# ontology-grouped CV leaks event identity on low-variability subsets (99.7%
# of RI/Low events span multiple ontologies -- ontology-grouping doesn't
# isolate whole events the way seqnames-grouping does, see memory
# project_splicing_ml_ontology_cv_event_leak.md) -- ablation conclusions on
# ontology-split numbers would be unreliable for exactly the subsets battery
# A cares most about. User-directed 2026-07-14: ablations run on seqnames only.
rule splicing_ml_ablation_classification:
    wildcard_constraints:
        group_col = "seqnames",
    input:
        src  = SPLICING_ML_SRC,
        data = "processed_data/aggregated_dt_filtered_{transcript_filter}.csv.gz",
    output:
        pkl  = "splicing_ml/output_ablation/{feature_groups}/{event_type}_{transcript_filter}_{variability}_{group_col}/splicing_ml_results_classification.pkl.gz",
        json = "splicing_ml/output_ablation/{feature_groups}/{event_type}_{transcript_filter}_{variability}_{group_col}/splicing_ml_results_classification.json.gz",
    log: "logs/splicing_ml_ablation_classification_{feature_groups}_{event_type}_{transcript_filter}_{variability}_{group_col}.log"
    threads: lambda wc: ml_resource(wc.event_type, wc.variability, "threads")
    resources:
        mem_mb          = lambda wc, attempt: ml_resource(wc.event_type, wc.variability, "mem_mb") * attempt,
        runtime         = lambda wc: ml_resource(wc.event_type, wc.variability, "runtime"),
        slurm_partition = _ml_partition,
        slurm_extra     = _ml_extra,
        qos             = _ml_qos,
        gres            = _ml_gres,
    params:
        wandb_login  = wandb_login_cmd(),
        wandb_args   = wandb_cli_args(),
        feature_args = lambda wc: wc.feature_groups.replace("+", " "),
    shell:
        """
        {params.wandb_login}PYTHONUNBUFFERED=1 mamba run --no-capture-output -n ihec-as python run_splicing_ml.py \
            --debug \
            --optuna \
            --data-path {input.data} \
            --output-dir splicing_ml/output_ablation/{wildcards.feature_groups}/{wildcards.event_type}_{wildcards.transcript_filter}_{wildcards.variability}_{wildcards.group_col} \
            --only-event-type {wildcards.event_type} \
            --only-transcript-filter {wildcards.transcript_filter} \
            --only-variability {wildcards.variability} \
            --only-group {wildcards.group_col} \
            --feature-groups {params.feature_args} \
            --skip-regression \
            --max-cores {threads} \
            {params.wandb_args} \
        > {log} 2>&1
        """


rule splicing_ml_ablation_regression:
    wildcard_constraints:
        group_col = "seqnames",
    input:
        src  = SPLICING_ML_SRC,
        data = "processed_data/aggregated_dt_filtered_{transcript_filter}.csv.gz",
    output:
        pkl  = "splicing_ml/output_ablation/{feature_groups}/{event_type}_{transcript_filter}_{variability}_{group_col}/splicing_ml_results_regression.pkl.gz",
        json = "splicing_ml/output_ablation/{feature_groups}/{event_type}_{transcript_filter}_{variability}_{group_col}/splicing_ml_results_regression.json.gz",
    log: "logs/splicing_ml_ablation_regression_{feature_groups}_{event_type}_{transcript_filter}_{variability}_{group_col}.log"
    threads: lambda wc: ml_resource(wc.event_type, wc.variability, "threads")
    resources:
        mem_mb          = lambda wc, attempt: ml_resource(wc.event_type, wc.variability, "mem_mb") * attempt,
        runtime         = lambda wc: ml_resource(wc.event_type, wc.variability, "runtime"),
        slurm_partition = _ml_partition,
        slurm_extra     = _ml_extra,
        qos             = _ml_qos,
        gres            = _ml_gres,
    params:
        wandb_login  = wandb_login_cmd(),
        wandb_args   = wandb_cli_args(),
        feature_args = lambda wc: wc.feature_groups.replace("+", " "),
    shell:
        """
        {params.wandb_login}PYTHONUNBUFFERED=1 mamba run --no-capture-output -n ihec-as python run_splicing_ml.py \
            --debug \
            --optuna \
            --data-path {input.data} \
            --output-dir splicing_ml/output_ablation/{wildcards.feature_groups}/{wildcards.event_type}_{wildcards.transcript_filter}_{wildcards.variability}_{wildcards.group_col} \
            --only-event-type {wildcards.event_type} \
            --only-transcript-filter {wildcards.transcript_filter} \
            --only-variability {wildcards.variability} \
            --only-group {wildcards.group_col} \
            --feature-groups {params.feature_args} \
            --skip-classification \
            --max-cores {threads} \
            {params.wandb_args} \
        > {log} 2>&1
        """


# ── Event-specific models: build → Tier-1 ridge screen → aggregate → Tier-2 EN ─
# Two-tier redesign (2026-07-22). Tier-1 (09s-*) is a fit-free closed-form ridge
# screen that decides significance for ALL events cheaply (feature-rotation null
# → qvalue FDR). Tier-2 (09-1 → 09zz elastic-net) runs ONLY on screen hits, real
# PSI, for feature interpretation. DAG:
#   build_feature_tables → event_screen → screen_aggregate → event_models
# NOTE on async SLURM arrays (same fire-and-forget convention as the rest of this
# pipeline): `event_screen` and `event_models` SUBMIT arrays and return; the
# per-event work runs afterwards. Run `screen_aggregate` only after the screen
# array has finished (it warns if per-event files are incomplete), and re-run
# `ml_analysis` only after the elastic-net array has finished.

# ── Step 09-1a: build per-event feature tables (Phase 1, build-only) ───────────
rule build_feature_tables:
    input:
        script           = "09-1-ml-local.R",
        local_glmnet     = "09zz-ml-event-glmnet-tidymodels.R",  # source()d by 09-1
        shared           = "09-ml-shared.R",
        aggregated_dt    = "processed_data/aggregated_dt_filtered_{transcript_filter}.csv.gz",
        keep_rows_manual = "processed_data/keep_rows_manual_{transcript_filter}.rds",
        sample_cols      = "processed_data/sample_cols_{transcript_filter}.rds",
        event_annotations_dt = "processed_data/event_annotations_dt_{transcript_filter}.csv.gz",
        psi_long_dt      = "processed_data/psi_long_dt_{transcript_filter}.csv.gz",
        event_gr         = "processed_data/event_gr_{transcript_filter}.rds",
        active_chromhmm  = "processed_data/activeChromHMM_{transcript_filter}.rds",
        chromhmm_hits    = "processed_data/chromhmm_hits_{transcript_filter}.rds",
        wgbs             = "sample_dts/WGBS_agg_{transcript_filter}.csv.gz",
        file_table       = "processed_data/file_table.csv.gz",
        rbp_wide         = "processed_data/rbp_wide_expression_{transcript_filter}.csv.gz",
        rbp_per_event    = "processed_data/rbp_per_event_{transcript_filter}.rds",
    output:
        session = "processed_data/session_09_1_ml_local_{transcript_filter}.rds",
        all_ids = "processed_data/event_glmnet_all_ids_{transcript_filter}.txt",
        cfg     = "processed_data/event_glmnet_cfg_{transcript_filter}.rds",
    log: "logs/09_1a_build_feature_tables_{transcript_filter}.log"
    threads: R("event_models", "threads")
    resources:
        mem_mb          = R("event_models", "mem_mb"),
        runtime         = R("event_models", "runtime"),
        slurm_partition = _partition("event_models"),
        slurm_extra     = _extra("event_models"),
        qos             = _qos("event_models"),
        gres            = _gres("event_models"),
    shell:
        """
        EPIATLAS_AS_ML_PHASE=build \
        TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        Rscript 09-1-ml-local.R > {log} 2>&1
        """

# ── Step 09s: Tier-1 fit-free ridge screen (fires SLURM array over all events) ─
rule event_screen:
    input:
        dispatch = "09s-dispatch.R",
        screen   = "09s-ridge-screen.R",
        shared   = "09-ml-shared.R",
        array_sh = "09s-ridge-screen-array.sh",
        cfg      = "processed_data/event_glmnet_cfg_{transcript_filter}.rds",
        all_ids  = "processed_data/event_glmnet_all_ids_{transcript_filter}.txt",
        session  = "processed_data/session_09_1_ml_local_{transcript_filter}.rds",
    output:
        dispatched = touch("processed_data/event_models/{transcript_filter}/.screen_dispatched"),
    log: "logs/09s_event_screen_{transcript_filter}.log"
    threads: R("analysis", "threads")
    resources:
        mem_mb          = R("analysis", "mem_mb"),
        runtime         = R("analysis", "runtime"),
        slurm_partition = _partition("analysis"),
        slurm_extra     = _extra("analysis"),
        qos             = _qos("analysis"),
        gres            = _gres("analysis"),
    shell:
        """
        Rscript 09s-dispatch.R {input.cfg} > {log} 2>&1
        """

# ── Step 09s: aggregate the screen → qvalues + Tier-1 hit list ────────────────
# Run only after the event_screen SLURM array has FINISHED (guarded: warns on
# incomplete per-event files). Produces the hit list the elastic-net runs on.
rule screen_aggregate:
    input:
        aggregate  = "09s-aggregate.R",
        cfg        = "processed_data/event_glmnet_cfg_{transcript_filter}.rds",
        dispatched = "processed_data/event_models/{transcript_filter}/.screen_dispatched",
    output:
        results = "processed_data/event_models/{transcript_filter}/screen_results.csv.gz",
        hits    = "processed_data/event_models/{transcript_filter}/tier1_hits_{transcript_filter}.txt",
    log: "logs/09s_screen_aggregate_{transcript_filter}.log"
    threads: R("analysis", "threads")
    resources:
        mem_mb          = R("analysis", "mem_mb"),
        runtime         = R("analysis", "runtime"),
        slurm_partition = _partition("analysis"),
        slurm_extra     = _extra("analysis"),
        qos             = _qos("analysis"),
        gres            = _gres("analysis"),
    shell:
        """
        Rscript 09s-aggregate.R {input.cfg} > {log} 2>&1
        """

# ── Step 09-1: Tier-2 elastic-net on screen hits (fires SLURM array) ──────────
rule event_models:
    input:
        script           = "09-1-ml-local.R",
        local_glmnet     = "09zz-ml-event-glmnet-tidymodels.R",  # source()d by 09-1
        shared           = "09-ml-shared.R",
        hits             = "processed_data/event_models/{transcript_filter}/tier1_hits_{transcript_filter}.txt",
        session          = "processed_data/session_09_1_ml_local_{transcript_filter}.rds",
        cfg              = "processed_data/event_glmnet_cfg_{transcript_filter}.rds",
        aggregated_dt    = "processed_data/aggregated_dt_filtered_{transcript_filter}.csv.gz",
        keep_rows_manual = "processed_data/keep_rows_manual_{transcript_filter}.rds",
        sample_cols      = "processed_data/sample_cols_{transcript_filter}.rds",
        event_annotations_dt = "processed_data/event_annotations_dt_{transcript_filter}.csv.gz",
        psi_long_dt      = "processed_data/psi_long_dt_{transcript_filter}.csv.gz",
        event_gr         = "processed_data/event_gr_{transcript_filter}.rds",
        active_chromhmm  = "processed_data/activeChromHMM_{transcript_filter}.rds",
        chromhmm_hits    = "processed_data/chromhmm_hits_{transcript_filter}.rds",
        wgbs             = "sample_dts/WGBS_agg_{transcript_filter}.csv.gz",
        file_table       = "processed_data/file_table.csv.gz",
        rbp_wide         = "processed_data/rbp_wide_expression_{transcript_filter}.csv.gz",
        rbp_per_event    = "processed_data/rbp_per_event_{transcript_filter}.rds",
    output:
        done    = touch("processed_data/event_models/{transcript_filter}/.done"),
    log: "logs/09_1_event_models_{transcript_filter}.log"
    threads: R("event_models", "threads")
    resources:
        mem_mb          = R("event_models", "mem_mb"),
        runtime         = R("event_models", "runtime"),
        slurm_partition = _partition("event_models"),
        slurm_extra     = _extra("event_models"),
        qos             = _qos("event_models"),
        gres            = _gres("event_models"),
    shell:
        """
        TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        Rscript 09-1-ml-local.R > {log} 2>&1
        """


# ── Step 09-2: ML analysis report ─────────────────────────────────────────────
rule ml_analysis:
    input:
        aggregated_dt   = f"processed_data/aggregated_dt_filtered_{PRIMARY}.csv.gz",
        session         = f"processed_data/session_09_1_ml_local_{PRIMARY}.rds",
        active_chromhmm = f"processed_data/activeChromHMM_{PRIMARY}.rds",
        event_models    = f"processed_data/event_models/{PRIMARY}/.done",
        screen_results  = f"processed_data/event_models/{PRIMARY}/screen_results.csv.gz",
        rmd             = "09-2-ml-local-new.Rmd",
        # global model results (all configs for primary filter)
        splicing_ml = [
            f"splicing_ml/output/{et}_{PRIMARY}_{var}_{gc}/splicing_ml_results_classification.pkl.gz"
            for et, tf, var, gc in ML_CONFIGS if tf == PRIMARY
        ],
    output:
        html = f"reports/09-2-ml-local-new_{PRIMARY}.html",
    log: f"logs/09_2_ml_analysis_{PRIMARY}.log"
    threads: R("analysis", "threads")
    resources:
        mem_mb          = R("analysis", "mem_mb"),
        runtime         = R("analysis", "runtime"),
        slurm_partition = _partition("analysis"),
        slurm_extra     = _extra("analysis"),
        qos             = _qos("analysis"),
        gres            = _gres("analysis"),
    shell:
        f"""
        TRANSCRIPT_FILTER={PRIMARY} \\
        Rscript -e "rmarkdown::render('09-2-ml-local-new.Rmd',
            output_file = normalizePath('{{output.html}}', mustWork = FALSE)
        )" > {{log}} 2>&1
        """


# ── Step 10: experimental events ──────────────────────────────────────────────
rule experimental_events:
    input:
        rmd = "10-experimental-events.Rmd",
        aggregated_dt         = f"processed_data/aggregated_dt_filtered_{PRIMARY}.csv.gz",
        event_annotations_dt  = f"processed_data/event_annotations_dt_{PRIMARY}.csv.gz",
        active_chromhmm       = f"processed_data/activeChromHMM_{PRIMARY}.rds",
        event_gr              = f"processed_data/event_gr_{PRIMARY}.rds",
        chromhmm_hits         = f"processed_data/chromhmm_hits_{PRIMARY}.rds",
        keep_rows_manual      = f"processed_data/keep_rows_manual_{PRIMARY}.rds",
        event_models          = f"processed_data/event_models/{PRIMARY}/.done",
    output:
        html = f"reports/10-experimental-events_{PRIMARY}.html",
    log: f"logs/10_experimental_events_{PRIMARY}.log"
    threads: R("analysis", "threads")
    resources:
        mem_mb          = R("analysis", "mem_mb"),
        runtime         = R("analysis", "runtime"),
        slurm_partition = _partition("analysis"),
        slurm_extra     = _extra("analysis"),
        qos             = _qos("analysis"),
        gres            = _gres("analysis"),
    shell:
        f"""
        TRANSCRIPT_FILTER={PRIMARY} \\
        Rscript -e "rmarkdown::render('10-experimental-events.Rmd',
            output_file = normalizePath('{{output.html}}', mustWork = FALSE)
        )" > {{log}} 2>&1
        """


# ── Utility: clean generated outputs ──────────────────────────────────────────
rule clean:
    shell:
        """
        rm -f processed_data/aggregated_dt_filtered_*.csv.gz
        rm -f processed_data/correlation_intrinsic_*.csv.gz
        rm -f processed_data/correlation_intrinsic_preproc_*.csv.gz
        rm -rf processed_data/event_models/*/
        rm -rf processed_data/session_09_1_ml_local_*.rds
        rm -f processed_data/event_glmnet_all_ids_*.txt
        rm -f processed_data/event_glmnet_cfg_*.rds
        rm -f processed_data/event_glmnet_ids_*.txt
        rm -rf reports/
        echo "Cleaned outputs. Upstream steps (01-04) preserved."
        """
