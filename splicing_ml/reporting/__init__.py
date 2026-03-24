"""Reporting package entry points and compatibility exports."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

from .html import (
    config_key,
    generate_html_reports,
    slugify_config_key,
    write_subset_html_report,
)
from .payload import build_task_plot_payload, important_params_table_rows

__all__ = [
    "config_key",
    "slugify_config_key",
    "build_task_plot_payload",
    "important_params_table_rows",
    "write_subset_html_report",
    "generate_html_reports",
    "generate_html_reports_from_results_file",
]


def generate_html_reports_from_results_file(
    results_file: str,
    output_dir: str,
    verbose: bool = False,
) -> None:
    """Generate reports from an existing results artifact (no retraining)."""
    p = Path(results_file)
    if p.suffixes[-2:] == [".json", ".gz"]:
        with gzip.open(p, "rt", encoding="utf-8") as f:
            payload = json.load(f)
    elif p.suffix == ".json":
        payload = json.loads(p.read_text(encoding="utf-8"))
    else:
        raise ValueError("Unsupported results file type; expected .json or .json.gz")

    all_results = payload.get("results", [])
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    generate_html_reports(all_results, out_dir, verbose=verbose)
