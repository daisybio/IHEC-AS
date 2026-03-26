"""Reporting package entry points and compatibility exports."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from .html import (
    config_key,
    generate_html_reports,
    slugify_config_key,
    write_subset_html_report,
)
from .payload import build_task_plot_payload, important_params_table_rows
from ..io_utils import load_json_gz

__all__ = [
    "config_key",
    "slugify_config_key",
    "build_task_plot_payload",
    "important_params_table_rows",
    "write_subset_html_report",
    "generate_html_reports",
    "generate_html_reports_from_results_file",
    "generate_html_reports_from_task_dir",
]


def generate_html_reports_from_results_file(
    results_file: str,
    output_dir: str,
    verbose: bool = False,
) -> None:
    """Generate reports from a single results artifact (no retraining).

    Accepts .json.gz or .json files — either the legacy combined artifact
    or the per-task artifacts produced by the current pipeline.
    """
    p = Path(results_file)
    if p.suffix not in {".gz", ".json"}:
        raise ValueError("Unsupported results file type; expected .json or .json.gz")
    payload = load_json_gz(p)
    all_results = payload.get("results", [])
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    generate_html_reports(all_results, out_dir, verbose=verbose)


def generate_html_reports_from_task_dir(
    output_dir: str,
    task: Literal["classification", "regression", "both"] = "both",
    verbose: bool = False,
) -> None:
    """Generate reports by auto-loading per-task result files from *output_dir*.

    Looks for ``splicing_ml_results_{task}.json.gz`` files written by the
    pipeline.  When *task* is ``"both"``, both files are loaded and merged so
    that classification and regression reports are generated together.

    Parameters
    ----------
    output_dir:
        Directory that contains the ``splicing_ml_results_*.json.gz`` files.
    task:
        Which task(s) to include: ``"classification"``, ``"regression"``, or
        ``"both"`` (default).
    verbose:
        Emit progress messages.
    """
    tasks = ["classification", "regression"] if task == "both" else [task]
    out_dir = Path(output_dir)
    all_results: list = []
    for t in tasks:
        p = out_dir / f"splicing_ml_results_{t}.json.gz"
        if not p.exists():
            if verbose:
                print(f"[reporting] {p.name} not found — skipping {t}")
            continue
        payload = load_json_gz(p)
        all_results.extend(payload.get("results", []))
        if verbose:
            print(f"[reporting] Loaded {len(payload.get('results', []))} {t} results from {p.name}")

    if not all_results:
        raise FileNotFoundError(
            f"No result files found in {out_dir} for task={task!r}. "
            "Expected splicing_ml_results_classification.json.gz and/or "
            "splicing_ml_results_regression.json.gz"
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    generate_html_reports(all_results, out_dir, verbose=verbose)
