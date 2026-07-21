"""Tests for HTML report payload construction in splicing_ml.reporting.payload."""

import numpy as np
import pytest

from splicing_ml.reporting.payload import (
    _build_metric_fold_data,
    _build_metric_heatmap,
    _build_transformed_feature_distributions,
    _direction_sort_key,
    _is_lower_is_better,
    build_task_plot_payload,
    important_params_table_rows,
)


def _fold_row(model_name, outer_fold, scores, y_true=None, y_pred=None, **extra):
    row = {
        "model_name": model_name,
        "outer_fold": outer_fold,
        "scores": scores,
        "y_true": y_true or [],
        "y_pred": y_pred or [],
    }
    row.update(extra)
    return row


class TestBuildMetricHeatmap:
    """Test the model x metric heatmap builder."""

    def test_shapes_and_values_match_model_row_metric_col_layout(self):
        """z/sd/text are indexed [model][metric], one row per model_order entry."""
        fold_rows = [
            _fold_row("linear", 0, {"rmse": 0.10, "ccc": 0.80}),
            _fold_row("linear", 1, {"rmse": 0.20, "ccc": 0.60}),
            _fold_row("xgb", 0, {"rmse": 0.05, "ccc": 0.90}),
            _fold_row("xgb", 1, {"rmse": 0.15, "ccc": 0.70}),
        ]
        hm = _build_metric_heatmap(fold_rows, ["linear", "xgb"], ["rmse", "ccc"])

        assert hm["models"] == ["linear", "xgb"]
        assert hm["metrics"] == ["rmse", "ccc"]
        assert len(hm["z"]) == 2, "one row per model"
        assert len(hm["z"][0]) == 2, "one column per metric"

        assert hm["z"][0][0] == pytest.approx(0.15)  # linear rmse mean
        assert hm["z"][1][0] == pytest.approx(0.10)  # xgb rmse mean
        assert hm["z"][0][1] == pytest.approx(0.70)  # linear ccc mean
        assert hm["sd"][0][0] == pytest.approx(np.std([0.10, 0.20]))
        assert "+/-" in hm["text"][0][0]
        assert hm["directions"] == [True, False], "rmse is lower-is-better, ccc is not"

    def test_missing_metric_yields_none_and_na_text(self):
        """A model with no valid scores for a metric gets None/NA, not a crash."""
        fold_rows = [_fold_row("linear", 0, {"rmse": 0.1})]
        hm = _build_metric_heatmap(fold_rows, ["linear"], ["rmse", "auroc"])

        assert hm["z"][0][0] == pytest.approx(0.1)
        assert hm["z"][0][1] is None
        assert hm["sd"][0][1] is None
        assert hm["text"][0][1] == "NA"

    def test_nan_scores_are_excluded_from_mean(self):
        """NaN score entries are dropped, not averaged in as zero."""
        fold_rows = [
            _fold_row("linear", 0, {"rmse": 0.1}),
            _fold_row("linear", 1, {"rmse": float("nan")}),
        ]
        hm = _build_metric_heatmap(fold_rows, ["linear"], ["rmse"])
        assert hm["z"][0][0] == pytest.approx(0.1)


class TestBuildMetricFoldData:
    """Test per-metric fold-level box-plot point/line construction."""

    def test_points_grouped_by_metric_and_model(self):
        fold_rows = [
            _fold_row("linear", 0, {"rmse": 0.1, "ccc": 0.8}),
            _fold_row("xgb", 0, {"rmse": 0.2, "ccc": 0.7}),
            _fold_row("linear", 1, {"rmse": 0.15, "ccc": 0.75}),
        ]
        out = _build_metric_fold_data(fold_rows, ["linear", "xgb"], ["rmse", "ccc"])

        assert set(out.keys()) == {"rmse", "ccc"}
        rmse_points = out["rmse"]["points"]
        assert {p["model_name"] for p in rmse_points} == {"linear", "xgb"}
        linear_rmse = [p["score"] for p in rmse_points if p["model_name"] == "linear"]
        assert sorted(linear_rmse) == pytest.approx([0.1, 0.15])

    def test_fold_lines_only_include_models_present_in_that_fold(self):
        """A fold missing a model should not synthesize a score for it."""
        fold_rows = [
            _fold_row("linear", 0, {"rmse": 0.1}),
            _fold_row("xgb", 0, {"rmse": 0.2}),
            _fold_row("linear", 1, {"rmse": 0.3}),  # xgb absent in fold 1
        ]
        out = _build_metric_fold_data(fold_rows, ["linear", "xgb"], ["rmse"])
        lines = {line["outer_fold"]: line for line in out["rmse"]["lines"]}

        assert lines[0]["x"] == ["linear", "xgb"]
        assert lines[1]["x"] == ["linear"], "xgb missing from fold 1 must not appear"

    def test_nan_score_excluded_from_points(self):
        fold_rows = [_fold_row("linear", 0, {"rmse": float("nan")})]
        out = _build_metric_fold_data(fold_rows, ["linear"], ["rmse"])
        assert out["rmse"]["points"] == []
        assert out["rmse"]["lines"] == []


class TestBuildTransformedFeatureDistributions:
    """Test pooling of post-preprocessing feature samples for distribution plots."""

    def test_pools_values_across_folds_and_models_per_feature(self):
        fold_rows = [
            _fold_row(
                "linear",
                0,
                {"rmse": 0.1},
                transformed_features={
                    "feature_names": ["num__a", "num__b"],
                    "values": [[1.0, 10.0], [2.0, 20.0]],
                },
            ),
            _fold_row(
                "xgb",
                0,
                {"rmse": 0.2},
                transformed_features={
                    "feature_names": ["num__a", "num__b"],
                    "values": [[3.0, 30.0]],
                },
            ),
        ]
        out = _build_transformed_feature_distributions(fold_rows)
        assert set(out.keys()) == {"num__a", "num__b"}
        assert sorted(out["num__a"]) == [1.0, 2.0, 3.0]
        assert sorted(out["num__b"]) == [10.0, 20.0, 30.0]

    def test_rows_without_transformed_features_are_skipped(self):
        fold_rows = [_fold_row("linear", 0, {"rmse": 0.1})]
        assert _build_transformed_feature_distributions(fold_rows) == {}

    def test_empty_fold_rows_returns_empty_dict(self):
        assert _build_transformed_feature_distributions([]) == {}

    def test_downsamples_large_pooled_samples(self):
        big_values = [[float(i)] for i in range(10000)]
        fold_rows = [
            _fold_row(
                "linear",
                0,
                {"rmse": 0.1},
                transformed_features={"feature_names": ["num__a"], "values": big_values},
            ),
        ]
        out = _build_transformed_feature_distributions(fold_rows)
        assert len(out["num__a"]) <= 3000


class TestBuildTaskPlotPayload:
    """Test the top-level per-task payload builder."""

    def test_no_fold_results_returns_safe_empty_payload(self):
        payload = build_task_plot_payload({"task": "regression", "fold_results": []})
        assert payload["model_order"] == []
        assert payload["metric_fold_data"] == {}
        assert payload["transformed_feature_distributions"] == {}
        assert "No fold results available" in payload["warnings"]

    def test_all_nan_primary_scores_returns_safe_empty_payload(self):
        task_result = {
            "task": "regression",
            "primary_metric": "rmse",
            "fold_results": [_fold_row("linear", 0, {"rmse": float("nan")})],
        }
        payload = build_task_plot_payload(task_result)
        assert payload["model_order"] == []
        assert payload["metric_fold_data"] == {}
        assert "All fold scores are NaN" in payload["warnings"]

    def test_classification_populates_roc_pr_and_confusion(self):
        fold_rows = [
            _fold_row(
                "linear",
                0,
                {"auroc": 0.9, "ba": 0.8},
                y_true=[0, 1, 0, 1],
                y_pred=[0.1, 0.9, 0.2, 0.8],
                threshold=0.5,
            ),
            _fold_row(
                "linear",
                1,
                {"auroc": 0.85, "ba": 0.75},
                y_true=[0, 1, 1, 0],
                y_pred=[0.2, 0.7, 0.6, 0.3],
                threshold=0.5,
            ),
        ]
        task_result = {
            "task": "classification",
            "primary_metric": "auroc",
            "fold_results": fold_rows,
            "response_distribution": {},
        }
        payload = build_task_plot_payload(task_result)

        assert payload["model_order"] == ["linear"]
        assert set(payload["metric_fold_data"].keys()) == {"auroc", "ba"}
        assert "linear" in payload["confusion_by_model"]
        assert "linear" in payload["roc_by_model"]
        assert "linear" in payload["pr_by_model"]
        assert payload["pr_prevalence"] == pytest.approx(0.5)
        # Regression-only heatmaps must stay empty for classification.
        assert payload["metric_heatmap_original"]["models"] == []
        assert payload["metric_heatmap"]["models"] == ["linear"]

    def test_regression_splits_original_and_logit_heatmaps_and_fold_data(self):
        fold_rows = [
            _fold_row(
                "linear",
                0,
                {"original_rmse": 0.1, "logit_rmse": 1.2},
                y_true=[0.2, 0.5],
                y_pred=[0.25, 0.45],
            ),
            _fold_row(
                "linear",
                1,
                {"original_rmse": 0.15, "logit_rmse": 1.3},
                y_true=[0.3, 0.6],
                y_pred=[0.28, 0.55],
            ),
        ]
        task_result = {
            "task": "regression",
            "primary_metric": "original_rmse",
            "fold_results": fold_rows,
            "response_distribution": {},
        }
        payload = build_task_plot_payload(task_result)

        # Mixed heatmap is suppressed in favor of scale-specific ones.
        assert payload["metric_heatmap"]["models"] == []
        assert payload["metric_heatmap_original"]["metrics"] == ["original_rmse"]
        assert payload["metric_heatmap_logit"]["metrics"] == ["logit_rmse"]
        assert set(payload["metric_fold_data"].keys()) == {
            "original_rmse",
            "logit_rmse",
        }

    def test_regression_backfills_logit_metrics_when_absent(self):
        """Older artifacts without logit_* scores get them computed from y_true/y_pred."""
        fold_rows = [
            _fold_row(
                "linear",
                0,
                {"rmse": 0.1, "mad": 0.05, "r2_rss": 0.9, "ccc": 0.95},
                y_true=[0.2, 0.5, 0.8],
                y_pred=[0.25, 0.45, 0.75],
            ),
        ]
        task_result = {
            "task": "regression",
            "primary_metric": "rmse",
            "fold_results": fold_rows,
            "response_distribution": {},
        }
        payload = build_task_plot_payload(task_result)

        # Lower-is-better metrics (mad, rmse) grouped before higher-is-better
        # ones (ccc, r2_rss), alphabetical within each group.
        assert payload["metric_heatmap_original"]["metrics"] == [
            "mad",
            "rmse",
            "ccc",
            "r2_rss",
        ]
        assert payload["metric_heatmap_original"]["directions"] == [
            True,
            True,
            False,
            False,
        ]
        assert any(
            m.startswith("logit_") for m in payload["metric_heatmap_logit"]["metrics"]
        )
        assert any(k.startswith("logit_") for k in payload["metric_fold_data"])

    def test_transformed_feature_distributions_populated_end_to_end(self):
        fold_rows = [
            _fold_row(
                "linear",
                0,
                {"original_rmse": 0.1},
                y_true=[0.2, 0.5],
                y_pred=[0.25, 0.45],
                transformed_features={
                    "feature_names": ["num__width"],
                    "values": [[1.0], [2.0]],
                },
            ),
        ]
        task_result = {
            "task": "regression",
            "primary_metric": "original_rmse",
            "fold_results": fold_rows,
            "response_distribution": {},
        }
        payload = build_task_plot_payload(task_result)
        assert payload["transformed_feature_distributions"] == {"num__width": [1.0, 2.0]}


class TestMetricDirection:
    """Test lower/higher-is-better classification used for heatmap coloring."""

    @pytest.mark.parametrize(
        "metric,expected",
        [
            ("rmse", True),
            ("mad", True),
            ("original_rmse", True),
            ("logit_mad", True),
            ("ccc", False),
            ("r2_rss", False),
            ("original_ccc", False),
            ("logit_r2_rss", False),
            ("auroc", False),
            ("balanced_accuracy", False),
            ("mcc", False),
        ],
    )
    def test_is_lower_is_better(self, metric, expected):
        assert _is_lower_is_better(metric) is expected

    def test_direction_sort_key_groups_lower_is_better_first(self):
        metrics = ["ccc", "mad", "r2_rss", "rmse"]
        assert sorted(metrics, key=_direction_sort_key) == [
            "mad",
            "rmse",
            "ccc",
            "r2_rss",
        ]

    def test_classification_metrics_all_sort_after_none_are_lower_is_better(self):
        metrics = ["mcc", "auroc", "balanced_accuracy", "f1", "auprc"]
        # No classification metric is lower-is-better, so direction grouping
        # is a no-op and alphabetical order is preserved.
        assert sorted(metrics, key=_direction_sort_key) == sorted(metrics)


class TestImportantParamsTableRows:
    def test_most_frequent_params_summarized_per_model(self):
        task_result = {
            "fold_results": [
                {
                    "model_name": "linear",
                    "tuning": {"best_params": {"C": 1.0}},
                },
                {
                    "model_name": "linear",
                    "tuning": {"best_params": {"C": 1.0}},
                },
                {
                    "model_name": "linear",
                    "tuning": {"best_params": {"C": 10.0}},
                },
            ]
        }
        rows = important_params_table_rows(task_result)
        assert len(rows) == 1
        row = rows[0]
        assert row["model_name"] == "linear"
        assert row["selected_in_folds"] == "3"
        assert row["frequency"] == "2"
        assert '"C": 1.0' in row["most_frequent_best_params"]

    def test_no_fold_results_returns_empty_list(self):
        assert important_params_table_rows({"fold_results": []}) == []
