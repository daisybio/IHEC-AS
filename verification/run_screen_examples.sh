#!/bin/bash
# Run the Tier-1 ridge screen on a few example events in a SANDBOX, leaving the
# production screen/ directory untouched.
#
#   bash verification/run_screen_examples.sh [reuse|fresh] [id ...]
#
#   reuse  (default) copy each event's existing per-event screen output into the
#          sandbox first, so the resume gate keeps `long` verbatim and only
#          short/local are computed  (~2 feature sets, faster)
#   fresh  compute all three feature sets from scratch; `long` should then
#          reproduce the production numbers exactly — a useful self-check
#
# Fit-free (closed-form ridge, no tuning), so this is fine on the login node.
set -euo pipefail

# Repo root is derived from this script's own location rather than hardcoded, so the check runs
# from any clone. PROJ/TF/CORES stay overridable from the environment.
PROJ="${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
TF="${TRANSCRIPT_FILTER:-biotype_filtered}"
CORES="${CORES:-4}"
# A fresh temp directory per run. The sandbox is wiped and rebuilt on every invocation anyway, so
# nothing is gained by a fixed path -- and `rm -rf` on a variable that could resolve somewhere
# unintended is not worth the convenience. Override with SCREEN_SANDBOX to keep the output around.
SANDBOX="${SCREEN_SANDBOX:-$(mktemp -d -t screen_examples.XXXXXXXX)}"

MODE="${1:-reuse}"; shift || true
IDS=("$@")
if [ ${#IDS[@]} -eq 0 ]; then
  # 2 top-effect (both tiny-n, expect wild nulls), 2 top-raw-R2 (fingerprint
  # suspects), 1 well-powered event for contrast
  IDS=(53808 16754 38458 7796 10000)
fi

cd "$PROJ"
module load r/4.2.1 2>/dev/null || true

PROD_SCREEN="processed_data/event_models/$TF/screen"
case "$SANDBOX" in
  */screen_examples.*|*/screen_examples) ;;
  *) echo "refusing to rm -rf a sandbox path that is not a screen_examples dir: $SANDBOX" >&2; exit 1 ;;
esac
rm -rf "$SANDBOX"; mkdir -p "$SANDBOX/screen"

if [ "$MODE" = "reuse" ]; then
  for id in "${IDS[@]}"; do
    for suf in _screen.csv.gz _screen_null.csv.gz; do
      [ -f "$PROD_SCREEN/${id}${suf}" ] && cp "$PROD_SCREEN/${id}${suf}" "$SANDBOX/screen/"
    done
  done
  echo "mode=reuse: copied existing long results for ${#IDS[@]} event(s) into the sandbox"
else
  echo "mode=fresh: computing all feature sets from scratch"
fi

# Clone the production cfg, pointing event_dir at the sandbox. Everything else
# (feature tables, session rds) is read-only, so it is shared.
# NB: export BEFORE the Rscript that reads these — otherwise Sys.getenv() returns ""
# and `set -e` kills the script on the resulting bad readRDS path.
export TF SANDBOX
Rscript -e '
cfg <- readRDS(sprintf("processed_data/event_glmnet_cfg_%s.rds", Sys.getenv("TF")))
cfg$event_dir <- Sys.getenv("SANDBOX")
saveRDS(cfg, file.path(Sys.getenv("SANDBOX"), "cfg.rds"))
' >/dev/null 2>&1

for id in "${IDS[@]}"; do
  echo "=== event $id ==="
  /usr/bin/time -f "  wall %es  maxRSS %MkB" \
    Rscript 09s-ridge-screen.R "$SANDBOX/cfg.rds" "$id" "$CORES" 2>&1 \
    | grep -E "Resuming|Done screen|Already computed|Error|wall " || true
done

echo
echo "=== results ==="
Rscript -e '
suppressPackageStartupMessages(library(data.table))
f <- list.files(file.path(Sys.getenv("SANDBOX"), "screen"),
                pattern = "_screen\\.csv\\.gz$", full.names = TRUE)
f <- f[!grepl("_screen_null", f)]
d <- rbindlist(lapply(f, fread), fill = TRUE)
setorder(d, ID, feature_set)
print(d[, .(ID, feature_set, n_samples, n_features, R_used,
            screen_R2 = round(screen_R2, 4),
            null_R2_mean = round(null_R2_mean, 4),
            effect = round(effect, 4),
            p_emp = round(p_emp, 4), note)])
cat("\nNOTE: p_emp here is RAW. q-values/hits require FDR over the FULL event set\n",
    "(09s-aggregate.R, grouped by transcript_filter x Event Type x feature_set) —\n",
    "a handful of events cannot give you significance.\n")
' 2>&1 | grep -vE "^Loading|^Attaching|masked|^$|R\.(oo|utils|methodsS3)|expand.grid|shift|trim|set mc.cores|^ +[A-Za-z]+,"
echo
echo "sandbox: $SANDBOX  (production screen/ untouched)"
