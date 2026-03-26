from __future__ import annotations

"""W&B tracking adapter with a NullTracker fallback.

This module keeps experiment tracking optional and non-fatal. Call sites can
use a single tracker object without guard conditionals.
"""

import contextlib
import warnings
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np


__all__ = ["NullTracker", "WandbTracker", "make_tracker"]


class NullTracker:
    """No-op tracker matching the W&B tracker interface."""

    def start_run(self, *args: Any, **kwargs: Any) -> None:
        return None

    def start_config_task(self, *args: Any, **kwargs: Any) -> None:
        return None

    def log_fold_result(self, *args: Any, **kwargs: Any) -> None:
        return None

    def log_config_task_summary(self, *args: Any, **kwargs: Any) -> None:
        return None

    def log_artifact(self, *args: Any, **kwargs: Any) -> None:
        return None

    def finish_config_task(self, *args: Any, **kwargs: Any) -> None:
        return None

    def finish_run(self, *args: Any, **kwargs: Any) -> None:
        return None


class WandbTracker:
    """W&B-backed tracker.

    Creates one parent orchestrator run and one child run per
    (subset configuration x task) workload.
    """

    def __init__(self, wandb_module: Any, project: str, entity: str | None = None):
        self._wandb = wandb_module
        self._project = project
        self._entity = entity
        self._parent_run = None
        self._child_run = None

    def start_run(self, run_cfg: Any) -> None:
        self._parent_run = self._wandb.init(
            project=self._project,
            entity=self._entity,
            job_type="orchestrator",
            config=asdict(run_cfg),
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
        group = getattr(self._parent_run, "name", None)
        run_name = f"{cfg.event_type}-{cfg.transcript_filter}-{cfg.variability}-{cfg.group_col}-{task}"
        self._child_run = self._wandb.init(
            project=self._project,
            entity=self._entity,
            group=group,
            job_type="subset-task",
            name=run_name,
            config={
                "subset_config": asdict(cfg),
                "task": task,
                "run_config": asdict(run_cfg),
                "n_samples": int(n_samples),
            },
            reinit=True,
        )
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
        x_test: Any,
        train_scores: dict[str, float] | None = None,
    ) -> None:
        if self._child_run is None:
            return

        # Prefix model metrics so multiple model families can share one run.
        payload = {
            f"{model_name}/{k}": float(v)
            for k, v in scores.items()
            if isinstance(v, (int, float)) and np.isfinite(float(v))
        }
        payload[f"{model_name}/fold_id"] = int(fold_id)
        # Log training scores (prefixed train/) to expose overfitting gap.
        if train_scores:
            for k, v in train_scores.items():
                if isinstance(v, (int, float)) and np.isfinite(float(v)):
                    payload[f"{model_name}/train_{k}"] = float(v)
        self._child_run.log(payload)

        with contextlib.suppress(Exception):
            self._child_run.log({f"{model_name}/tuning": tuning_info})

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
            model = estimator.named_steps.get("model", estimator)
            prep = estimator.named_steps.get("prep")
            feature_names = None
            if prep is not None:
                feature_names = prep.get_feature_names_out()
            self._wandb.sklearn.plot_feature_importances(model, feature_names)

    def log_config_task_summary(self, result: dict[str, Any]) -> None:
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

    def log_artifact(self, path: str | Path, name: str, artifact_type: str) -> None:
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
        if self._child_run is None:
            return
        with contextlib.suppress(Exception):
            self._child_run.finish()
        self._child_run = None

    def finish_run(self, all_results: list[dict[str, Any]]) -> None:
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

    return WandbTracker(
        wandb_module=wandb,
        project=str(getattr(run_cfg, "wandb_project", "splicing-ml")),
        entity=getattr(run_cfg, "wandb_entity", None),
    )
