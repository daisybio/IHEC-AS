"""
IHEC-AS Snakemake workflow.

Run (dry-run):
    mamba run -n ihec-as snakemake --profile profiles/slurm -n

Run (SLURM):
    mamba run -n ihec-as snakemake --profile profiles/slurm

Run single filter (e.g. biotype_filtered only):
    mamba run -n ihec-as snakemake --profile profiles/slurm \
        --config transcript_filters='["biotype_filtered"]'

Steps:
    01  gather-data          →  processed_data/file_table.csv.gz
    02  SUPPA2-analysis      →  splicing_analysis/filtered_transcript_ids.rds
    02z rmats-analysis       →  splicing_analysis/rmats/{tf}/event_{et}.{psi,jc.csv.gz}
    03  prepare-aggregation  →  processed_data/{aggregateOver,file_metadata,sample_cols,keep_rows_manual}.rds, processed_data/ijc_sjc_dt.csv.gz, event_annotations_dt.csv.gz, psi_long_dt.csv.gz
    04  aggregate-WGBS       →  processed_data/sample_dts/WGBS_agg.csv.gz
    05  create-aggregated-dt →  processed_data/aggregated_dt_filtered.csv.gz  (all transcript_filters as rows)
    06  correlation          →  processed_data/correlation_intrinsic.csv.gz
    08  splicing_ml          →  splicing_ml/output/{et}_{tf}_{var}_{gc}/
    09-1 event-models        →  processed_data/event_models/{tf}/.done
    09-2 ml-analysis         →  reports/09-2-ml-local-new_{primary}.html
    10  experimental-events  →  reports/10-experimental-events.html
"""

from itertools import product

configfile: "config/snakemake_config.yaml"

localrules: all, clean

# ── Wildcard values ────────────────────────────────────────────────────────────
TRANSCRIPT_FILTERS = config["transcript_filters"]
EVENT_TYPES        = config["event_types"]
VARIABILITIES      = config["variabilities"]
GROUP_COLS         = config["group_cols"]
PRIMARY            = config["primary_filter"]

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

wildcard_constraints:
    transcript_filter = "|".join(TRANSCRIPT_FILTERS),
    event_type        = "|".join(EVENT_TYPES),
    variability       = "|".join(VARIABILITIES),
    group_col         = "|".join(GROUP_COLS),


# ── Helper: resource lookup ────────────────────────────────────────────────────
def R(rule_key, field):
    return config["resources"][rule_key][field]

def ml_resource(et, var, field):
    """Map (event_type, variability) to resource tier."""
    if et == "SE" and var == "both":
        tier = "splicing_ml_L"
    elif et == "RI" and var in ("High", "Low"):
        tier = "splicing_ml_S"
    else:
        tier = "splicing_ml_M"
    return config["resources"][tier][field]

def _partition(rule_key):
    """Return SLURM partition for a rule based on its gpu flag."""
    return (config["slurm_gpu_partition"] if config["resources"][rule_key]["gpu"]
            else config["slurm_cpu_partition"])

def _extra(rule_key):
    """Return extra sbatch args including GPU gres if needed."""
    base = "--mail-type=FAIL --mail-user=quirin.manz@tum.de"
    if config["resources"][rule_key]["gpu"]:
        g = config
        return (f"--gres={g['slurm_gpu_gres']} --qos={g['slurm_gpu_qos']} "
                f"--exclude={g['slurm_gpu_exclude']} {base}")
    return base

def _ml_partition(wc):
    gpu = ml_resource(wc.event_type, wc.variability, "gpu")
    return config["slurm_gpu_partition"] if gpu else config["slurm_cpu_partition"]

def _ml_extra(wc):
    base = "--mail-type=FAIL --mail-user=quirin.manz@tum.de"
    if ml_resource(wc.event_type, wc.variability, "gpu"):
        g = config
        return (f"--gres={g['slurm_gpu_gres']} --qos={g['slurm_gpu_qos']} "
                f"--exclude={g['slurm_gpu_exclude']} {base}")
    return base


# ── ChIP BigWig file list (for aggregate_chip_one wildcard) ───────────────────
CHIP_BW, = glob_wildcards("/nfs/data3/IHEC/ChIP-Seq/{bw}.pval.signal.bigwig")


# ── Rule all ──────────────────────────────────────────────────────────────────
rule all:
    input:
        # Combined aggregated data and correlations (all transcript_filters as rows)
        "processed_data/aggregated_dt_filtered.csv.gz",
        "processed_data/correlation_intrinsic.csv.gz",
        # ML global models (all configs)
        [f"splicing_ml/output/{et}_{tf}_{var}_{gc}/splicing_ml_results_classification.pkl.gz"
         for et, tf, var, gc in ML_CONFIGS],
        # Event-specific models (primary filter only — extend if needed)
        f"processed_data/event_models/{PRIMARY}/.done",
        # Final analysis reports
        f"reports/09-2-ml-local-new_{PRIMARY}.html",
        "reports/10-experimental-events.html",


# ── Step 01: gather data ───────────────────────────────────────────────────────
rule gather_data:
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
    shell:
        """
        Rscript -e "rmarkdown::render('01-gather-data.Rmd',
            output_file = normalizePath('reports/01-gather-data.html', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 02: SUPPA2 isoform quants + GTF annotation filter ────────────────────
rule suppa_analysis:
    input:
        file_table = "processed_data/file_table.csv.gz",
    output:
        filtered_tx_ids = "splicing_analysis/filtered_transcript_ids.rds",
        tpm_expr        = "splicing_analysis/suppa/tpm_expressions.tsv.gz",
        gencode_gtfs    = expand(
            "splicing_analysis/gencode.v29.{tf}.gtf",
            tf=TRANSCRIPT_FILTERS,
        ),
    log: "logs/02_suppa_analysis.log"
    threads: R("suppa_analysis", "threads")
    resources:
        mem_mb          = R("suppa_analysis", "mem_mb"),
        runtime         = R("suppa_analysis", "runtime"),
        slurm_partition = _partition("suppa_analysis"),
        slurm_extra     = _extra("suppa_analysis"),
    shell:
        """
        Rscript -e "rmarkdown::render('02-SUPPA2-analysis.Rmd',
            output_file = normalizePath('reports/02-SUPPA2-analysis.html', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 02z: rMATS PSI + junction-count filter (PLAN §4.12d, §4.8) ──────────
rule rmats_analysis:
    input:
        file_table = "processed_data/file_table.csv.gz",
    output:
        psi_files = expand(
            "splicing_analysis/rmats/{tf}/event_{et}.psi",
            tf=TRANSCRIPT_FILTERS, et=EVENT_TYPES,
        ),
        jc_files = expand(
            "splicing_analysis/rmats/{tf}/event_{et}.jc.csv.gz",
            tf=TRANSCRIPT_FILTERS, et=EVENT_TYPES,
        ),
    log: "logs/02z_rmats_analysis.log"
    threads: R("rmats_analysis", "threads")
    resources:
        mem_mb          = R("rmats_analysis", "mem_mb"),
        runtime         = R("rmats_analysis", "runtime"),
        slurm_partition = _partition("rmats_analysis"),
        slurm_extra     = _extra("rmats_analysis"),
    shell:
        """
        Rscript -e "rmarkdown::render('02z-rmats-analysis.Rmd',
            output_file = normalizePath('reports/02z-rmats-analysis.html', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 03: prepare aggregation workspace ────────────────────────────────────
rule prepare_aggregation:
    input:
        psi_files = expand(
            "splicing_analysis/rmats/{tf}/event_{et}.psi",
            tf=TRANSCRIPT_FILTERS, et=EVENT_TYPES,
        ),
        jc_files = expand(
            "splicing_analysis/rmats/{tf}/event_{et}.jc.csv.gz",
            tf=TRANSCRIPT_FILTERS, et=EVENT_TYPES,
        ),
        filtered_tx_ids = "splicing_analysis/filtered_transcript_ids.rds",
        gencode_gtfs = expand(
            "splicing_analysis/gencode.v29.{tf}.gtf",
            tf=TRANSCRIPT_FILTERS,
        ),
    output:
        aggregateOver_bed = "processed_data/aggregateOver.bed",
        sample_cols      = "processed_data/sample_cols.rds",
        keep_rows_manual = "processed_data/keep_rows_manual.rds",
        ijc_sjc_dt            = "processed_data/ijc_sjc_dt.csv.gz",
        event_annotations_dt  = "processed_data/event_annotations_dt.csv.gz",
        psi_long_dt           = "processed_data/psi_long_dt.csv.gz",
        pangolin_events       = "processed_data/pangolin_events.csv",
        ss5_fasta             = "processed_data/5ss.fasta",
        ss5up_fasta           = "processed_data/5ss_up.fasta",
        ss3_fasta             = "processed_data/3ss.fasta",
        ss3down_fasta         = "processed_data/3ss_down.fasta",
    log: "logs/03_prepare_aggregation.log"
    threads: R("prepare_aggregation", "threads")
    resources:
        mem_mb          = R("prepare_aggregation", "mem_mb"),
        runtime         = R("prepare_aggregation", "runtime"),
        slurm_partition = _partition("prepare_aggregation"),
        slurm_extra     = _extra("prepare_aggregation"),
    shell:
        """
        Rscript -e "rmarkdown::render('03-prepare-aggregation.Rmd',
            output_file = normalizePath('reports/03-prepare-aggregation.html', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 04-3: MaxEntScan splice-site scores ─────────────────────────────────
rule maxentscan_scores:
    input:
        ss3      = "processed_data/3ss.fasta",
        ss3down  = "processed_data/3ss_down.fasta",
        ss5      = "processed_data/5ss.fasta",
        ss5up    = "processed_data/5ss_up.fasta",
    output:
        score3     = "processed_data/3scores.txt",
        score3down = "processed_data/3down_scores.txt",
        score5     = "processed_data/5scores.txt",
        score5up   = "processed_data/5up_scores.txt",
    log: "logs/04-3_maxentscan_scores.log"
    threads: R("maxentscan_scores", "threads")
    resources:
        mem_mb          = R("maxentscan_scores", "mem_mb"),
        runtime         = R("maxentscan_scores", "runtime"),
        slurm_partition = _partition("maxentscan_scores"),
        slurm_extra     = _extra("maxentscan_scores"),
    shell:
        """
        mamba run -n ihec-as bash 04-3-maxentscan-scores.sh > {log} 2>&1
        """


# ── Step 04-4: Pangolin splice-site scores ───────────────────────────────────
rule pangolin_scores:
    input:
        events_csv = "processed_data/pangolin_events.csv",
        fasta      = "data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
    output:
        scores_csv = "processed_data/pangolin_scores.csv",
    log: "logs/04-4_pangolin_scores.log"
    threads: R("pangolin_scores", "threads")
    resources:
        mem_mb          = R("pangolin_scores", "mem_mb"),
        runtime         = R("pangolin_scores", "runtime"),
        slurm_partition = _partition("pangolin_scores"),
        slurm_extra     = _extra("pangolin_scores"),
    shell:
        """
        bash 04-4-pangolin-scores.sh > {log} 2>&1
        """


# ── Step 04b: aggregate ChIP-Seq signal ──────────────────────────────────────
# Submitted as SLURM arrays via snakemake-executor-plugin-slurm:
#   profile groups aggregate_chip_one → chip_array (100 tasks per array).
# CHIP_BW discovered at parse time via glob_wildcards.

rule aggregate_chip_one:
    input:
        bw  = "/nfs/data3/IHEC/ChIP-Seq/{bw}.pval.signal.bigwig",
        bed = "processed_data/aggregateOver.bed",
    output:
        tab = "sample_dts/{bw}.pval.signal.bigwig.tab.gz",
    log: "logs/chip/{bw}.log"
    threads: R("aggregate_chip_one", "threads")
    resources:
        mem_mb          = R("aggregate_chip_one", "mem_mb"),
        runtime         = R("aggregate_chip_one", "runtime"),
        slurm_partition = _partition("aggregate_chip_one"),
        slurm_extra     = "--mail-type=NONE",
    shell:
        """
        bigWigAverageOverBed {input.bw} {input.bed} \
            sample_dts/{wildcards.bw}.pval.signal.bigwig.tab -minMax \
            > {log} 2>&1
        gzip -f sample_dts/{wildcards.bw}.pval.signal.bigwig.tab
        """

rule aggregate_chip:
    input:
        expand("sample_dts/{bw}.pval.signal.bigwig.tab.gz", bw=CHIP_BW),
    output:
        done = touch("sample_dts/chip_agg.done"),
    localrule: True


# ── Step 04: aggregate WGBS ───────────────────────────────────────────────────
rule aggregate_wgbs:
    input:
        aggregateOver = "processed_data/aggregateOver.bed",
    output:
        wgbs = "sample_dts/WGBS_agg.csv.gz",
    log: "logs/04_aggregate_wgbs.log"
    threads: R("aggregate_wgbs", "threads")
    resources:
        mem_mb          = R("aggregate_wgbs", "mem_mb"),
        runtime         = R("aggregate_wgbs", "runtime"),
        slurm_partition = _partition("aggregate_wgbs"),
        slurm_extra     = _extra("aggregate_wgbs"),
    shell:
        "Rscript 04-1-aggregate-WGBS-matrix.R > {log} 2>&1"


# ── Step 05: create aggregated dataset (all transcript_filters as rows) ────────
rule create_aggregated_dt:
    input:
        keep_rows_manual      = "processed_data/keep_rows_manual.rds",
        file_table            = "processed_data/file_table.csv.gz",
        ijc_sjc_dt            = "processed_data/ijc_sjc_dt.csv.gz",
        event_annotations_dt  = "processed_data/event_annotations_dt.csv.gz",
        psi_long_dt           = "processed_data/psi_long_dt.csv.gz",
        wgbs                  = "sample_dts/WGBS_agg.csv.gz",
        chip_done             = "sample_dts/chip_agg.done",
        qc_covariates         = "processed_data/qc_flag_covariates.csv",
        filtered_tx_ids       = "splicing_analysis/filtered_transcript_ids.rds",
        maxentscan_score3     = "processed_data/3scores.txt",
        maxentscan_score3down = "processed_data/3down_scores.txt",
        maxentscan_score5     = "processed_data/5scores.txt",
        maxentscan_score5up   = "processed_data/5up_scores.txt",
        ss3_fasta             = "processed_data/3ss.fasta",
        ss5_fasta             = "processed_data/5ss.fasta",
        pangolin_scores       = "processed_data/pangolin_scores.csv",
    output:
        csv  = "processed_data/aggregated_dt_filtered.csv.gz",
        html = "reports/05-create-aggregated-dt.html",
    log: "logs/05_create_aggregated_dt.log"
    threads: R("create_aggregated_dt", "threads")
    resources:
        mem_mb          = R("create_aggregated_dt", "mem_mb"),
        runtime         = R("create_aggregated_dt", "runtime"),
        slurm_partition = _partition("create_aggregated_dt"),
        slurm_extra     = _extra("create_aggregated_dt"),
    shell:
        """
        Rscript -e "rmarkdown::render('05-create-aggregated-dt.Rmd',
            output_file = normalizePath('{output.html}', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 06: correlation analysis ────────────────────────────────────────────
rule correlation:
    input:
        aggregated_dt = "processed_data/aggregated_dt_filtered.csv.gz",
        file_table    = "processed_data/file_table.csv.gz",
    output:
        corr_raw     = "processed_data/correlation_intrinsic.csv.gz",
        corr_preproc = "processed_data/correlation_intrinsic_preproc.csv.gz",
        html         = "reports/06-correlation.html",
    log: "logs/06_correlation.log"
    threads: R("correlation", "threads")
    resources:
        mem_mb          = R("correlation", "mem_mb"),
        runtime         = R("correlation", "runtime"),
        slurm_partition = _partition("correlation"),
        slurm_extra     = _extra("correlation"),
    shell:
        """
        Rscript -e "rmarkdown::render('06-correlation.Rmd',
            output_file = normalizePath('{output.html}', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 08: splicing_ml global models ────────────────────────────────────────
rule splicing_ml:
    input:
        data = "processed_data/aggregated_dt_filtered.csv.gz",
    output:
        pkl  = "splicing_ml/output/{event_type}_{transcript_filter}_{variability}_{group_col}/splicing_ml_results_classification.pkl.gz",
        json = "splicing_ml/output/{event_type}_{transcript_filter}_{variability}_{group_col}/splicing_ml_results_classification.json.gz",
    log: "logs/splicing_ml_{event_type}_{transcript_filter}_{variability}_{group_col}.log"
    threads: lambda wc: ml_resource(wc.event_type, wc.variability, "threads")
    resources:
        mem_mb          = lambda wc, attempt: ml_resource(wc.event_type, wc.variability, "mem_mb") * attempt,
        runtime         = lambda wc: ml_resource(wc.event_type, wc.variability, "runtime"),
        slurm_partition = _ml_partition,
        slurm_extra     = _ml_extra,
    shell:
        """
        mamba run -n ihec-as python run_splicing_ml.py \
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
        > {log} 2>&1
        """


# ── Step 09-1: event-specific glmnet models ───────────────────────────────────
rule event_models:
    input:
        aggregated_dt    = "processed_data/aggregated_dt_filtered.csv.gz",
        keep_rows_manual = "processed_data/keep_rows_manual.rds",
        sample_cols      = "processed_data/sample_cols.rds",
        file_table       = "processed_data/file_table.csv.gz",
    output:
        session = "processed_data/session_09_1_ml_local_{transcript_filter}.rds",
        done    = touch("processed_data/event_models/{transcript_filter}/.done"),
    log: "logs/09_1_event_models_{transcript_filter}.log"
    threads: R("event_models", "threads")
    resources:
        mem_mb          = R("event_models", "mem_mb"),
        runtime         = R("event_models", "runtime"),
        slurm_partition = _partition("event_models"),
        slurm_extra     = _extra("event_models"),
    shell:
        """
        TRANSCRIPT_FILTER={wildcards.transcript_filter} \
        Rscript 09-1-ml-local.R > {log} 2>&1
        """


# ── Step 09-2: ML analysis report ─────────────────────────────────────────────
rule ml_analysis:
    input:
        aggregated_dt = "processed_data/aggregated_dt_filtered.csv.gz",
        session       = f"processed_data/session_09_1_ml_local_{PRIMARY}.rds",
        event_models  = f"processed_data/event_models/{PRIMARY}/.done",
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
    shell:
        """
        Rscript -e "rmarkdown::render('09-2-ml-local-new.Rmd',
            output_file = normalizePath('{output.html}', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Step 10: experimental events ──────────────────────────────────────────────
rule experimental_events:
    input:
        aggregated_dt = "processed_data/aggregated_dt_filtered.csv.gz",
        sample_cols   = "processed_data/sample_cols.rds",
        session       = f"processed_data/session_09_1_ml_local_{PRIMARY}.rds",
        event_models  = f"processed_data/event_models/{PRIMARY}/.done",
    output:
        html = "reports/10-experimental-events.html",
    log: "logs/10_experimental_events.log"
    threads: R("analysis", "threads")
    resources:
        mem_mb          = R("analysis", "mem_mb"),
        runtime         = R("analysis", "runtime"),
        slurm_partition = _partition("analysis"),
        slurm_extra     = _extra("analysis"),
    shell:
        """
        Rscript -e "rmarkdown::render('10-experimental-events.Rmd',
            output_file = normalizePath('{output.html}', mustWork = FALSE)
        )" > {log} 2>&1
        """


# ── Utility: clean generated outputs ──────────────────────────────────────────
rule clean:
    shell:
        """
        rm -f processed_data/aggregated_dt_filtered.csv.gz
        rm -f processed_data/correlation_intrinsic.csv.gz
        rm -f processed_data/correlation_intrinsic_preproc.csv.gz
        rm -rf processed_data/event_models/*/
        rm -rf processed_data/session_09_1_ml_local_*.rds
        rm -rf reports/
        echo "Cleaned outputs. Upstream steps (01-04) preserved."
        """
