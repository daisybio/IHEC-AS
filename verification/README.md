# Verification scripts

Checks that re-derive claims the pipeline's own comments assert. They are not part of any Snakemake
rule and nothing depends on their output — run them by hand after the change each one guards.

All read production output **read-only**. None writes into `processed_data/`, `splicing_analysis/` or
any event `screen/` directory.

| script | checks | run it after |
| --- | --- | --- |
| `verify_psi_independent_controls.R` | the PSI-independent control features: that the complete grid (`aggregated_dt_{tf}.csv.gz`) agrees with the PSI-masked table on every shared column, and that `09-ml-shared.R`'s `build_full_event_rows()` assembles the full-cohort rows it is supposed to. Takes an event count, e.g. `12`. ~10 min, mostly I/O. | any re-run of stage `05`, and any change to the feature-table build. `FEATURE_TABLE_VERSION` 2 and RI's full 818-control pool both depend on this holding. |
| `verify_09_2_reported_numbers.R` | 54 figures that `09-2-ml-local-new.Rmd`'s §1 prose states as literals, against the real `screen_results.csv.gz` and `09-1`'s session rds. Fit-free, ~5 min. | any screen re-run, and any edit to that prose. The report hardcodes measured values, so they go stale silently. |
| `check_pooled_z_strata.R` | the pooled-z stratification in `09s-aggregate.R` — that the reference is standardised within `(Event Type x feature_set)` and not pooled across them. Pooling reads as calibrated at ratio 1.01 while being 1.43x liberal for SE and 0.41x conservative for RI, because SE outnumbers RI ~18:1. | any change to `09s-aggregate.R`'s z construction. `p_pooled_z` is the floor-free statistic the candidate set rests on. |
| `run_screen_examples.sh` | one event end-to-end through the **real** Tier-1 worker in a sandbox, e.g. `bash run_screen_examples.sh reuse 20541`. Production `screen/` is untouched — the worker hard-`stop()`s on a `normalizePath` comparison if an output path resolves inside it. | any edit to `09s-ridge-screen.R` or `09-ml-shared.R`, before committing to a full re-screen. |

## If you followed a `.claude/scratch/<name>` path here

These four moved from `.claude/scratch/` to `verification/` on 2026-09-03 (`a0a982b`), so they could
be tracked and satisfy the code-at-submission requirement. Dated records under `revision/` still name
the old path — correctly, as a record of where the script was when that check was run — and are
deliberately not being rewritten. Known instances: `revision-status-code-verified.md`,
`file-changes/09-confound-projection-instability.md`, `file-changes/09-psi-independent-controls.md`,
`file-changes/11-paper-figures.Rmd.md`, `file-changes/09-2-restructure-at-zero-hits.md`. The scripts
are here; the file names did not change.

Scripts still under `.claude/scratch/` are a different case: they are untracked, do not ship, and are
development scaffolding (`test_*` fixtures, one-off investigations). A path that resolves there today
is not stale — it is simply not part of the release.

## A note on paths in the pipeline's comments

Three files still refer to these scripts by their former location, `.claude/scratch/<name>`:

* `05-create-aggregated-dt.Rmd` → `verify_psi_independent_controls.R`
* `09-ml-shared.R` → `run_screen_examples.sh`
* `09s-aggregate.R` → `check_pooled_z_strata.R`

The paths are stale; the scripts are here. They were **deliberately not corrected**, because all three
files are declared Snakemake inputs and a comment-only edit to any of them invalidates far more than
the comment is worth: `05-create-aggregated-dt.Rmd` re-runs `create_aggregated_dt` and everything
downstream of it (the 117 GB feature-table build, then all 34,146 screen jobs); `09-ml-shared.R` is an
input of **five** rules including that same build; `09s-aggregate.R` re-runs the `screen_aggregate`
checkpoint and, through it, the 4,768-job floor sweep.

They will be corrected the next time each of those stages re-runs for a substantive reason. The same
discipline applies to a handful of other cosmetic fixes in this pipeline, for the same cost reason.
`09-2-ml-local-new.Rmd`'s reference **was** updated, because `ml_analysis` is a leaf and re-rendering
it costs one report.
