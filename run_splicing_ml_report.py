"""Generate plotting/evaluation HTML from pretrained run artifacts.

This entrypoint does not retrain models. It reads existing pipeline results
(JSON/JSON.GZ) and regenerates reports.
"""

from __future__ import annotations

import argparse

from splicing_ml.reporting import generate_html_reports_from_results_file


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate reports from pretrained run outputs"
    )
    p.add_argument(
        "--results-file",
        default="processed_data/ml_splicing_outputs/splicing_ml_results.json.gz",
        help="Path to existing results artifact (.json.gz or .json)",
    )
    p.add_argument(
        "--output-dir",
        default="processed_data/ml_splicing_outputs_pretrained_report",
        help="Directory where HTML reports will be written",
    )
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    generate_html_reports_from_results_file(
        results_file=args.results_file,
        output_dir=args.output_dir,
        verbose=args.verbose,
    )
    print(f"Report generation complete: {args.output_dir}")


if __name__ == "__main__":
    main()
