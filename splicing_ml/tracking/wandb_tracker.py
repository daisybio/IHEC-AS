from __future__ import annotations

"""W&B tracking adapter with a NullTracker fallback.

This module keeps experiment tracking optional and non-fatal. Call sites can
use a single tracker object without guard conditionals.
"""

import contextlib
import json
import os
import warnings
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np


__all__ = ["NullTracker", "WandbTracker", "make_tracker"]


class NullTracker:
    """No-op tracker matching the W&B tracker interface."""

    def start_run(self, *args: Any, **kwargs: Any) -> None:
        """Start run."""
        return None

    def start_config_task(self, *args: Any, **kwargs: Any) -> None:
        """Start config task."""
        return None

    def log_fold_result(self, *args: Any, **kwargs: Any) -> None:
        """Log fold result."""
        return None

    def log_config_task_summary(self, *args: Any, **kwargs: Any) -> None:
        """Log config task summary."""
        return None

    def log_artifact(self, *args: Any, **kwargs: Any) -> None:
        """Log artifact."""
        return None

    def finish_config_task(self, *args: Any, **kwargs: Any) -> None:
        """Finish config task."""
        return None

    def finish_run(self, *args: Any, **kwargs: Any) -> None:
        """Finish run."""
        return None


class WandbTracker:
    """W&B-backed tracker.

    Creates one parent orchestrator run and one child run per
    (subset configuration x task) workload.
    """

    def __init__(self, wandb_module: Any, run_cfg: Any):
        """Initialize a WandbTracker instance."""
        self._wandb = wandb_module
        self._project = str(getattr(run_cfg, "wandb_project", "splicing-ml"))
        self._entity = getattr(run_cfg, "wandb_entity", None)
        self._run_cfg = run_cfg
        self._parent_run = None
        self._child_run = None
        self._child_context: dict[str, Any] = {}
        self._fold_rows: list[dict[str, Any]] = []
        self._cv_results_by_model: dict[str, list[dict[str, Any]]] = {}

    @staticmethod
    def _safe_float(value: Any) -> float | None:
        """Internal helper for safe float."""
        if not isinstance(value, (int, float)):
            return None
        fval = float(value)
        if not np.isfinite(fval):
            return None
        return fval

    def _extract_tuning_summary(self, tuning_info: dict[str, Any]) -> dict[str, Any]:
        """Internal helper for extract tuning summary."""
        best_score = self._safe_float(tuning_info.get("best_score"))
        best_params = tuning_info.get("best_params", {})
        if not isinstance(best_params, dict):
            best_params = {}
        param_keys = sorted(str(k) for k in best_params.keys())
        summary: dict[str, Any] = {
            "best_score": best_score,
            "best_param_count": int(len(best_params)),
            "best_param_keys": ",".join(param_keys),
            "best_params_json": json.dumps(best_params, sort_keys=True, default=str),
        }
        if "tree_es_applied" in tuning_info:
            summary["tree_es_applied"] = bool(tuning_info.get("tree_es_applied"))
        if "tree_es_n_estimators" in tuning_info:
            n_estimators = tuning_info.get("tree_es_n_estimators")
            if isinstance(n_estimators, (int, float)):
                summary["tree_es_n_estimators"] = int(n_estimators)
        if "beta_convergence_risk" in tuning_info:
            summary["beta_convergence_risk"] = json.dumps(
                tuning_info.get("beta_convergence_risk"),
                sort_keys=True,
                default=str,
            )
        return summary

    def _log_fold_subrun(
        self, payload: dict[str, Any], fold_row: dict[str, Any]
    ) -> None:
        """Internal helper for log fold subrun."""
        if not bool(getattr(self._run_cfg, "wandb_fold_subruns", False)):
            return
        if self._child_run is None:
            return
        child_name = str(getattr(self._child_run, "name", "subset-task"))
        run_name = (
            f"{child_name}-fold{int(fold_row['fold_id']):02d}-{fold_row['model_name']}"
        )
        subrun = None
        with contextlib.suppress(Exception):
            subrun = self._wandb.init(
                project=self._project,
                entity=self._entity,
                group=str(getattr(self._child_run, "group", None)),
                job_type="subset-task-fold",
                name=run_name,
                tags=["fold", str(fold_row.get("task", ""))],
                config={
                    "fold_id": int(fold_row["fold_id"]),
                    "model_name": str(fold_row["model_name"]),
                    "task": str(fold_row["task"]),
                    "subset_context": dict(self._child_context),
                },
                reinit=True,
            )
            subrun.log(payload)
            subrun.summary.update(
                {
                    "primary_metric": fold_row.get("primary_metric"),
                    "primary_metric_value": fold_row.get("primary_metric_value"),
                    "n_test_samples": fold_row.get("n_test_samples"),
                }
            )
            subrun.finish()

    def _log_fold_table(self) -> None:
        """Internal helper for log fold table."""
        if self._child_run is None:
            return
        if not bool(getattr(self._run_cfg, "wandb_log_fold_table", True)):
            return
        if not self._fold_rows:
            return

        columns = [
            "fold_id",
            "model_name",
            "task",
            "primary_metric",
            "primary_metric_value",
            "train_primary_metric_value",
            "generalization_gap",
            "n_test_samples",
            "n_train_samples",
            "best_inner_score",
            "best_param_count",
            "best_param_keys",
            "best_params_json",
            "baseline_primary_metric_value",
            "delta_vs_baseline_primary",
        ]
        table = self._wandb.Table(columns=columns)
        for row in self._fold_rows:
            table.add_data(
                row.get("fold_id"),
                row.get("model_name"),
                row.get("task"),
                row.get("primary_metric"),
                row.get("primary_metric_value"),
                row.get("train_primary_metric_value"),
                row.get("generalization_gap"),
                row.get("n_test_samples"),
                row.get("n_train_samples"),
                row.get("best_inner_score"),
                row.get("best_param_count"),
                row.get("best_param_keys"),
                row.get("best_params_json"),
                row.get("baseline_primary_metric_value"),
                row.get("delta_vs_baseline_primary"),
            )
        with contextlib.suppress(Exception):
            self._child_run.log({"folds/by_model": table})

    def start_run(self, run_cfg: Any) -> None:
        """Start run."""
        tags = ["splicing-ml", "orchestrator"]
        if bool(getattr(run_cfg, "smoke_mode", False)):
            tags.append("smoke")
        tags.append(f"log_level:{getattr(run_cfg, 'log_level', 'info')}")
        self._parent_run = self._wandb.init(
            project=self._project,
            entity=self._entity,
            job_type="orchestrator",
            config=asdict(run_cfg),
            tags=tags,
            reinit=True,
        )

    def start_config_task(
        self,
        cfg: Any,
        task: str,
        run_cfg: Any,
        n_samples: int,
        prep_details: dict[str, Any],
    ) -> None:
        """Start config task."""
        group = getattr(self._parent_run, "name", None)
        run_name = (
            f"{cfg.event_type}-{cfg.transcript_filter}-{cfg.variability}-"
            f"{cfg.group_col}-{task}"
        )
        tags = [
            f"event:{cfg.event_type}",
            f"tx:{cfg.transcript_filter}",
            f"var:{cfg.variability}",
            f"group:{cfg.group_col}",
            f"task:{task}",
            "subset-task",
        ]
        if bool(getattr(run_cfg, "smoke_mode", False)):
            tags.append("smoke")
        if bool(getattr(run_cfg, "wandb_fold_subruns", False)):
            tags.append("fold-subruns")
        self._child_run = self._wandb.init(
            project=self._project,
            entity=self._entity,
            group=group,
            job_type="subset-task",
            name=run_name,
            tags=tags,
            config={
                "subset_config": asdict(cfg),
                "task": task,
                "run_config": asdict(run_cfg),
                "n_samples": int(n_samples),
            },
            reinit=True,
        )
        self._child_context = {
            "event_type": str(cfg.event_type),
            "transcript_filter": str(cfg.transcript_filter),
            "variability": str(cfg.variability),
            "group_col": str(cfg.group_col),
            "task": str(task),
        }
        self._fold_rows = []
        self._cv_results_by_model = {}
        # Without this, every logged metric is plotted against W&B's internal
        # auto-incrementing _step (i.e. call order), not the fold number — since
        # each fold triggers 2-3 separate .log() calls (scores, tuning_summary,
        # sklearn diagnostics) across multiple models, the default charts end up
        # as scrambled/sawtooth lines that don't compare across models at all.
        # Declaring fold_id as the step metric makes every chart use it as the
        # x-axis instead, regardless of call order.
        with contextlib.suppress(Exception):
            self._child_run.define_metric("fold_id")
            self._child_run.define_metric("*", step_metric="fold_id")
        with contextlib.suppress(Exception):
            prep_table = self._wandb.Table(columns=["key", "value"])
            for key, value in prep_details.items():
                prep_table.add_data(str(key), str(value))
            self._child_run.log({"preprocessing/details": prep_table})

    def log_fold_result(
        self,
        fold_id: int,
        model_name: str,
        scores: dict[str, float],
        tuning_info: dict[str, Any],
        y_true: np.ndarray | list[float],
        y_pred: np.ndarray | list[float],
        estimator: Any,
        task: str,
        x_train: Any,
        train_scores: dict[str, float] | None = None,
        primary_metric: str | None = None,
        baseline_scores: dict[str, float] | None = None,
        delta_vs_baseline: dict[str, float] | None = None,
    ) -> None:
        """Log fold result."""
        if self._child_run is None:
            return

        # Prefix model metrics so multiple model families can share one run.
        payload = {
            f"{model_name}/{k}": float(v)
            for k, v in scores.items()
            if isinstance(v, (int, float)) and np.isfinite(float(v))
        }
        payload[f"{model_name}/fold_id"] = int(fold_id)
        # Unprefixed fold_id (paired with the define_metric step_metric wiring
        # in start_config_task) so every chart in this run uses fold number as
        # its x-axis instead of call order.
        payload["fold_id"] = int(fold_id)
        # Log training scores (prefixed train/) to expose overfitting gap.
        if train_scores:
            for k, v in train_scores.items():
                if isinstance(v, (int, float)) and np.isfinite(float(v)):
                    payload[f"{model_name}/train_{k}"] = float(v)
        if bool(getattr(self._run_cfg, "wandb_log_baseline_metrics", True)):
            for k, v in (baseline_scores or {}).items():
                if isinstance(v, (int, float)) and np.isfinite(float(v)):
                    payload[f"{model_name}/baseline_{k}"] = float(v)
            for k, v in (delta_vs_baseline or {}).items():
                if isinstance(v, (int, float)) and np.isfinite(float(v)):
                    payload[f"{model_name}/delta_vs_baseline_{k}"] = float(v)
        self._child_run.log(payload)

        if bool(getattr(self._run_cfg, "wandb_log_tuning_details", True)):
            # Previously also logged the raw tuning_info dict wholesale
            # (f"{model_name}/tuning") — nested dicts (best_params objects,
            # full cv_results arrays) don't render as a chart, just clutter as
            # unreadable panels. The flat summary below is the useful part;
            # _log_cv_result_tables() already exposes the full cv_results as a
            # proper Table.
            with contextlib.suppress(Exception):
                tuning_summary = self._extract_tuning_summary(tuning_info)
                self._child_run.log({f"{model_name}/tuning_summary": tuning_summary})

        primary_key = str(primary_metric or "")
        primary_value = (
            self._safe_float(scores.get(primary_key)) if primary_key else None
        )
        train_primary = (
            self._safe_float((train_scores or {}).get(primary_key))
            if primary_key
            else None
        )
        gap = (
            float(train_primary - primary_value)
            if train_primary is not None and primary_value is not None
            else None
        )
        tuning_summary_for_row = self._extract_tuning_summary(tuning_info)
        baseline_primary = (
            self._safe_float((baseline_scores or {}).get(primary_key))
            if primary_key
            else None
        )
        delta_primary = (
            self._safe_float((delta_vs_baseline or {}).get(primary_key))
            if primary_key
            else None
        )
        fold_row = {
            "fold_id": int(fold_id),
            "model_name": str(model_name),
            "task": str(task),
            "primary_metric": primary_key,
            "primary_metric_value": primary_value,
            "train_primary_metric_value": train_primary,
            "generalization_gap": gap,
            "n_test_samples": int(len(y_true)) if y_true is not None else None,
            "n_train_samples": int(getattr(x_train, "shape", [len(x_train)])[0]),
            "best_inner_score": tuning_summary_for_row.get("best_score"),
            "best_param_count": tuning_summary_for_row.get("best_param_count"),
            "best_param_keys": tuning_summary_for_row.get("best_param_keys"),
            "best_params_json": tuning_summary_for_row.get("best_params_json"),
            "baseline_primary_metric_value": baseline_primary,
            "delta_vs_baseline_primary": delta_primary,
        }
        self._fold_rows.append(fold_row)
        cv_results = (
            tuning_info.get("cv_results") if isinstance(tuning_info, dict) else None
        )
        if cv_results and isinstance(cv_results, dict):
            entry = {"outer_fold": int(fold_id), **cv_results}
            self._cv_results_by_model.setdefault(model_name, []).append(entry)
        self._log_fold_subrun(payload, fold_row)

        # Optional sklearn diagnostics; never fail pipeline execution.
        if task == "classification":
            yt = np.asarray(y_true, dtype=int)
            yp = np.asarray(y_pred, dtype=float)
            # wandb.sklearn functions expect full predict_proba output (n, 2).
            y_probas = np.column_stack([1.0 - yp, yp])
            labels = ["0", "1"]

            with contextlib.suppress(Exception):
                hard = (yp >= 0.5).astype(int)
                self._child_run.log(
                    {
                        f"{model_name}/confusion_matrix": self._wandb.plot.confusion_matrix(
                            y_true=yt.tolist(),
                            preds=hard.tolist(),
                            class_names=labels,
                        )
                    }
                )

            with contextlib.suppress(Exception):
                if len(np.unique(yt)) >= 2:
                    self._child_run.log(
                        {
                            f"{model_name}/roc": self._wandb.plot.roc_curve(
                                yt, y_probas, labels
                            )
                        }
                    )

            with contextlib.suppress(Exception):
                if len(np.unique(yt)) >= 2:
                    self._child_run.log(
                        {
                            f"{model_name}/pr_curve": self._wandb.plot.pr_curve(
                                yt, y_probas, labels
                            )
                        }
                    )

        with contextlib.suppress(Exception):
            if not bool(getattr(self._run_cfg, "wandb_fold_subruns", False)):
                model = estimator.named_steps.get("model", estimator)
                prep = estimator.named_steps.get("prep")
                feature_names = None
                if prep is not None:
                    feature_names = prep.get_feature_names_out()
                self._wandb.sklearn.plot_feature_importances(model, feature_names)

    def _log_model_comparison(self, primary_metric: str) -> None:
        """Log outer-fold model comparison: bar chart (mean) + fold-detail table."""
        if self._child_run is None or not self._fold_rows:
            return
        from collections import defaultdict

        scores_by_model: dict[str, list[float]] = defaultdict(list)
        for row in self._fold_rows:
            v = row.get("primary_metric_value")
            if v is not None:
                scores_by_model[str(row["model_name"])].append(float(v))
        if not scores_by_model:
            return

        summary_table = self._wandb.Table(columns=["model", "mean", "std", "n_folds"])
        for model_name, vals in sorted(scores_by_model.items()):
            summary_table.add_data(
                model_name,
                float(np.mean(vals)),
                float(np.std(vals)),
                len(vals),
            )

        detail_table = self._wandb.Table(columns=["model", "fold_id", "score"])
        for row in self._fold_rows:
            v = row.get("primary_metric_value")
            if v is not None:
                detail_table.add_data(
                    str(row["model_name"]),
                    int(row["fold_id"]),
                    float(v),
                )

        # One multi-line chart (all models, x=fold_id) instead of relying on
        # the per-model auto-generated single-line panels from log_fold_result
        # — those land at different points on W&B's default x-axis whenever
        # models are logged in different .log() calls, so they were never
        # visually comparable without this. Sort each model's (fold_id, score)
        # pairs by fold_id first so the line doesn't zigzag if folds completed
        # out of order (parallel fitting).
        fold_score_pairs_by_model: dict[str, list[tuple[int, float]]] = defaultdict(
            list
        )
        for row in self._fold_rows:
            v = row.get("primary_metric_value")
            if v is not None:
                fold_score_pairs_by_model[str(row["model_name"])].append(
                    (int(row["fold_id"]), float(v))
                )
        model_names_sorted = sorted(fold_score_pairs_by_model)
        line_chart = None
        with contextlib.suppress(Exception):
            xs, ys = [], []
            for m in model_names_sorted:
                pairs = sorted(fold_score_pairs_by_model[m])
                xs.append([p[0] for p in pairs])
                ys.append([p[1] for p in pairs])
            line_chart = self._wandb.plot.line_series(
                xs=xs,
                ys=ys,
                keys=model_names_sorted,
                title=f"{primary_metric} by fold, per model",
                xname="fold_id",
            )

        payload = {
            "comparison/model_summary": summary_table,
            "comparison/model_fold_detail": detail_table,
            "comparison/bar_chart": self._wandb.plot.bar(
                summary_table,
                "model",
                "mean",
                title=f"Model comparison — mean {primary_metric} (outer folds)",
            ),
        }
        if line_chart is not None:
            payload["comparison/line_chart"] = line_chart
        with contextlib.suppress(Exception):
            self._child_run.log(payload)

    def _log_cv_result_tables(self) -> None:
        """Log one aggregated inner-CV results Table per model (across all outer folds)."""
        if self._child_run is None or not self._cv_results_by_model:
            return
        for model_name, fold_entries in self._cv_results_by_model.items():
            with contextlib.suppress(Exception):
                # Determine split column names from first entry.
                split_keys: list[str] = []
                for entry in fold_entries:
                    split_keys = entry.get("split_keys", [])
                    if split_keys:
                        break
                columns = [
                    "outer_fold",
                    "rank",
                    "mean_score",
                    "std_score",
                    "params",
                ] + split_keys
                table = self._wandb.Table(columns=columns)
                for entry in fold_entries:
                    outer_fold = entry["outer_fold"]
                    ranks = entry.get("rank_test_score", [])
                    means = entry.get("mean_test_score", [])
                    stds = entry.get("std_test_score", [])
                    params = entry.get("params_json", [])
                    split_scores = entry.get("split_scores", {})
                    n = entry.get("n_candidates", len(means))
                    for i in range(n):
                        split_vals = [
                            split_scores.get(k, [None] * n)[i] for k in split_keys
                        ]
                        table.add_data(
                            outer_fold,
                            ranks[i] if i < len(ranks) else None,
                            means[i] if i < len(means) else None,
                            stds[i] if i < len(stds) else None,
                            params[i] if i < len(params) else None,
                            *split_vals,
                        )
                self._child_run.log({f"{model_name}/cv_results": table})

    def log_config_task_summary(self, result: dict[str, Any]) -> None:
        """Log config task summary."""
        if self._child_run is None:
            return
        with contextlib.suppress(Exception):
            self._child_run.summary.update(
                {
                    "status": result.get("status"),
                    "task": result.get("task"),
                    "primary_metric": result.get("primary_metric"),
                    "primary_metric_mean": result.get("primary_metric_mean"),
                    "primary_metric_std": result.get("primary_metric_std"),
                    "primary_metric_bootstrap_ci": result.get(
                        "primary_metric_bootstrap_ci"
                    ),
                    "n_samples": result.get("n_samples"),
                    "n_fold_results": len(result.get("fold_results", [])),
                }
            )
        with contextlib.suppress(Exception):
            self._log_fold_table()
        with contextlib.suppress(Exception):
            self._log_model_comparison(str(result.get("primary_metric", "")))
        with contextlib.suppress(Exception):
            self._log_cv_result_tables()

    def log_artifact(self, path: str | Path, name: str, artifact_type: str) -> None:
        """Log artifact."""
        if self._parent_run is None:
            return
        p = Path(path)
        if not p.exists():
            return
        with contextlib.suppress(Exception):
            artifact = self._wandb.Artifact(name=name, type=artifact_type)
            artifact.add_file(str(p))
            self._parent_run.log_artifact(artifact)

    def finish_config_task(self) -> None:
        """Finish config task."""
        if self._child_run is None:
            return
        with contextlib.suppress(Exception):
            self._child_run.finish()
        self._child_run = None
        self._fold_rows = []
        self._cv_results_by_model = {}
        self._child_context = {}

    def finish_run(self, all_results: list[dict[str, Any]]) -> None:
        """Finish run."""
        if self._parent_run is None:
            return
        with contextlib.suppress(Exception):
            summary_table = self._wandb.Table(
                columns=[
                    "event_type",
                    "transcript_filter",
                    "variability",
                    "group_col",
                    "task",
                    "status",
                    "primary_metric",
                    "primary_metric_mean",
                ]
            )
            for row in all_results:
                cfg = row.get("config", {})
                summary_table.add_data(
                    str(cfg.get("event_type", "")),
                    str(cfg.get("transcript_filter", "")),
                    str(cfg.get("variability", "")),
                    str(cfg.get("group_col", "")),
                    str(row.get("task", "")),
                    str(row.get("status", "")),
                    str(row.get("primary_metric", "")),
                    float(row.get("primary_metric_mean", np.nan)),
                )
            self._parent_run.log({"summary/by_config_task": summary_table})

        with contextlib.suppress(Exception):
            self._parent_run.finish()
        self._parent_run = None


def make_tracker(run_cfg: Any) -> NullTracker | WandbTracker:
    """Create a tracker from runtime config with graceful fallback."""
    if not bool(getattr(run_cfg, "use_wandb", False)):
        return NullTracker()

    try:
        import wandb
    except Exception:
        warnings.warn(
            "use_wandb=True but wandb is unavailable; proceeding with NullTracker",
            RuntimeWarning,
            stacklevel=2,
        )
        return NullTracker()

    project = str(getattr(run_cfg, "wandb_project", "splicing-ml")).strip()
    if not project:
        raise ValueError("W&B project must be non-empty when --wandb is enabled")

    require_auth = bool(getattr(run_cfg, "wandb_require_auth", True))
    has_env_api_key = bool(os.environ.get("WANDB_API_KEY"))
    has_stored_auth = False
    if not has_env_api_key:
        # Accept credentials loaded via `wandb login` (e.g., ~/.netrc).
        with contextlib.suppress(Exception):
            has_stored_auth = bool(getattr(wandb.Api(), "api_key", None))

    has_auth = has_env_api_key or has_stored_auth
    if require_auth and not has_auth:
        raise RuntimeError(
            "W&B is enabled but no credentials were found. "
            "Run `wandb login`, export WANDB_API_KEY, or disable strict auth "
            "with --wandb-no-require-auth."
        )
    if not require_auth and not has_auth:
        warnings.warn(
            "use_wandb=True but W&B credentials are missing; proceeding with NullTracker "
            "because strict auth is disabled",
            RuntimeWarning,
            stacklevel=2,
        )
        return NullTracker()

    return WandbTracker(wandb_module=wandb, run_cfg=run_cfg)
