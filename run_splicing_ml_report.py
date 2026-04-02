#!/usr/bin/env python
"""Generate plotting/evaluation HTML from pretrained run artifacts.

This entrypoint does not retrain models. It reads existing per-task pipeline
results (JSON/JSON.GZ) and regenerates reports.

Usage examples
--------------
# Both tasks (default) — auto-detected from --output-dir:
python run_splicing_ml_report.py --output-dir splicing_ml/output/SE_transcripts_both_seqnames

# One task only:
python run_splicing_ml_report.py --output-dir ... --task classification

# Explicit file path (legacy combined artifact or any single file):
python run_splicing_ml_report.py --results-file path/to/splicing_ml_results_regression.json.gz
"""

from __future__ import annotations

import argparse

from splicing_ml.reporting import (
    generate_html_reports_from_results_file,
    generate_html_reports_from_task_dir,
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description="Generate reports from pretrained run outputs"
    )

    source = p.add_mutually_exclusive_group()
    source.add_argument(
        "--results-file",
        default=None,
        help=(
            "Path to a specific results artifact (.json.gz or .json). "
            "When given, --task is ignored."
        ),
    )
    source.add_argument(
        "--task",
        choices=["classification", "regression", "both"],
        default="both",
        help=(
            "Which task results to load from --output-dir. "
            "Looks for splicing_ml_results_{task}.json.gz. "
            "Default: both (loads and merges classification + regression)."
        ),
    )

    p.add_argument(
        "--output-dir",
        default="splicing_ml/output/ml_splicing_outputs",
        help="Directory containing result artifacts and where HTML reports are written.",
    )
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def main() -> None:
    """Run the command-line entrypoint."""
    args = parse_args()

    if args.results_file is not None:
        generate_html_reports_from_results_file(
            results_file=args.results_file,
            output_dir=args.output_dir,
            verbose=args.verbose,
        )
    else:
        generate_html_reports_from_task_dir(
            output_dir=args.output_dir,
            task=args.task,
            verbose=args.verbose,
        )

    print(f"Report generation complete: {args.output_dir}")


if __name__ == "__main__":
    main()
