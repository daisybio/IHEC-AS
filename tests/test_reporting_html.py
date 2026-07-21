"""Tests for HTML report file writing in splicing_ml.reporting.html."""

import json

import pytest

from splicing_ml.reporting.html import (
    _PLOTLY_JS_HELPERS,
    config_key,
    generate_html_reports,
    slugify_config_key,
    write_subset_html_report,
)


def _make_task_result(event_type="RI", transcript_filter="biotype_filtered",
                       variability="Low", group_col="seqnames", task="regression"):
    return {
        "config": {
            "event_type": event_type,
            "transcript_filter": transcript_filter,
            "variability": variability,
            "group_col": group_col,
        },
        "task": task,
        "status": "ok",
        "primary_metric": "original_rmse",
        "n_samples": 42,
        "fold_results": [
            {
                "model_name": "linear",
                "outer_fold": 0,
                "scores": {"original_rmse": 0.1, "logit_rmse": 1.0},
                "y_true": [0.2, 0.5, 0.8],
                "y_pred": [0.25, 0.45, 0.75],
            },
            {
                "model_name": "linear",
                "outer_fold": 1,
                "scores": {"original_rmse": 0.15, "logit_rmse": 1.1},
                "y_true": [0.3, 0.6],
                "y_pred": [0.28, 0.55],
            },
        ],
        "response_distribution": {"psi_sample": [0.1, 0.5, 0.9]},
    }


class TestConfigKeyHelpers:
    def test_config_key_extracts_tuple_in_order(self):
        result = _make_task_result()
        assert config_key(result) == ("RI", "biotype_filtered", "Low", "seqnames")

    def test_config_key_defaults_missing_fields_to_na(self):
        assert config_key({"config": {}}) == ("NA", "NA", "NA", "NA")

    def test_slugify_replaces_unsafe_characters(self):
        slug = slugify_config_key(("RI", "biotype filtered", "Low/High", "seqnames"))
        assert " " not in slug and "/" not in slug
        assert slug == "RI_biotype_filtered_Low_High_seqnames"


class TestWriteSubsetHtmlReport:
    def test_writes_html_and_sibling_json_payload(self, tmp_path):
        out_path = tmp_path / "subset_report_test.html"
        write_subset_html_report(
            ("RI", "biotype_filtered", "Low", "seqnames"),
            [_make_task_result()],
            out_path,
            run_datetime="2026-07-14T00:00:00Z",
        )

        assert out_path.exists()
        json_path = out_path.with_suffix(".json")
        assert json_path.exists()

        html_text = out_path.read_text(encoding="utf-8")
        assert json_path.name in html_text, "HTML must reference the sibling JSON payload"
        assert "2026-07-14T00:00:00Z" in html_text

        payload = json.loads(json_path.read_text(encoding="utf-8"))
        assert len(payload) == 1
        task_payload = payload[0]

        # Regression: mixed heatmap suppressed, scale-specific ones populated
        # with models on one axis and metrics on the other (order is a JS-side
        # transpose concern; here we just check both axes are present/consistent).
        hm = task_payload["metric_heatmap_original"]
        assert hm["models"] == ["linear"]
        assert hm["metrics"] == ["original_rmse"]
        assert len(hm["z"]) == len(hm["models"])
        assert len(hm["z"][0]) == len(hm["metrics"])

        # Per-metric box-plot data must be present for every heatmap metric.
        assert "original_rmse" in task_payload["metric_fold_data"]
        assert "logit_rmse" in task_payload["metric_fold_data"]

    def test_creates_parent_directories(self, tmp_path):
        out_path = tmp_path / "nested" / "dir" / "report.html"
        write_subset_html_report(
            ("RI", "biotype_filtered", "Low", "seqnames"),
            [_make_task_result()],
            out_path,
        )
        assert out_path.exists()

    def test_missing_run_datetime_displays_na(self, tmp_path):
        out_path = tmp_path / "report.html"
        write_subset_html_report(
            ("RI", "biotype_filtered", "Low", "seqnames"),
            [_make_task_result()],
            out_path,
            run_datetime=None,
        )
        assert "N/A" in out_path.read_text(encoding="utf-8")


class TestGenerateHtmlReports:
    def test_groups_results_by_config_key_into_separate_reports(self, tmp_path):
        results = [
            _make_task_result(event_type="RI"),
            _make_task_result(event_type="SE"),
            _make_task_result(event_type="RI", task="classification"),
        ]
        generate_html_reports(results, tmp_path)

        reports_dir = tmp_path / "reports"
        html_files = sorted(p.name for p in reports_dir.glob("*.html"))
        assert html_files == [
            "subset_report_RI_biotype_filtered_Low_seqnames.html",
            "subset_report_SE_biotype_filtered_Low_seqnames.html",
        ]

        # Both RI tasks (regression + classification) land in the same report.
        ri_json = json.loads(
            (reports_dir / "subset_report_RI_biotype_filtered_Low_seqnames.json")
            .read_text(encoding="utf-8")
        )
        assert {t["task"] for t in ri_json} == {"regression", "classification"}

    def test_no_results_creates_empty_reports_dir(self, tmp_path):
        generate_html_reports([], tmp_path)
        assert (tmp_path / "reports").is_dir()
        assert list((tmp_path / "reports").glob("*.html")) == []


class TestPlotlyJsHelpersRegressions:
    """Guard the embedded JS against two Plotly rendering bugs found by hand.

    These only grep the JS source string (no JS runtime here), but they pin
    the fix so a careless edit to the template doesn't silently reintroduce
    either bug.
    """

    def test_confusion_matrix_yaxis_is_reversed_for_diagonal_layout(self):
        """Without this, the y-axis renders True 0 below True 1, putting the
        correct-prediction cells (TN, TP) on the anti-diagonal instead of the
        conventional top-left-to-bottom-right diagonal."""
        cm_fn_start = _PLOTLY_JS_HELPERS.index("function plotConfusionMatrix")
        cm_fn = _PLOTLY_JS_HELPERS[cm_fn_start:]
        assert "autorange: 'reversed'" in cm_fn


def test_box_grid_two_pass_ordering_in_template(tmp_path):
    """The box-plot grid must build every container div in one pass, then
    call Plotly.newPlot in a second pass -- interleaving the two (creating a
    div and immediately plotting into it before its siblings exist) was the
    root cause of clipped titles / missing categories on the first couple of
    per-metric box plots."""
    out_path = tmp_path / "report.html"
    write_subset_html_report(
        ("RI", "biotype_filtered", "Low", "seqnames"),
        [_make_task_result()],
        out_path,
    )
    html_text = out_path.read_text(encoding="utf-8")

    div_creation_idx = html_text.index("boxDivByMetric[metricName] = boxDiv")
    first_plot_call_idx = html_text.index("const boxDiv = boxDivByMetric[metricName]")
    assert div_creation_idx < first_plot_call_idx, (
        "container-div creation loop must appear before the Plotly.newPlot "
        "loop in the rendered template"
    )


def test_confusion_matrices_are_gathered_in_one_section_after_the_task_loop(tmp_path):
    """Confusion matrices must render as one shared section at the bottom of
    the report (appended to `container`, after every per-task card), not
    interleaved inside each task's own card."""
    out_path = tmp_path / "report.html"
    write_subset_html_report(
        ("RI", "biotype_filtered", "Low", "seqnames"),
        [_make_task_result(task="classification")],
        out_path,
    )
    html_text = out_path.read_text(encoding="utf-8")

    # Accumulation happens per-task, inside the render loop.
    accumulate_idx = html_text.index("allConfusionMatrices.push")
    # The shared section is built after the loop closes and appended to the
    # top-level `container`, not to any per-task `card`.
    loop_close_idx = html_text.index(
        "// ── Confusion matrices, gathered in one section at the report bottom"
    )
    section_append_idx = html_text.index("container.appendChild(cmSection)")

    assert accumulate_idx < loop_close_idx < section_append_idx, (
        "confusion-matrix accumulation must happen inside the per-task loop, "
        "and the shared section must be built and appended after it"
    )
    # The old per-card header/wrapper must be gone -- otherwise matrices
    # would render both inline and in the gathered section.
    assert "Confusion matrices by model" not in html_text
    assert "Confusion matrices (all models)" in html_text


class TestTransformedFeatureDistributions:
    """Test the post-preprocessing feature-distribution section."""

    def test_section_present_when_payload_has_distributions(self, tmp_path):
        task_result = _make_task_result()
        task_result["fold_results"][0]["transformed_features"] = {
            "feature_names": ["num__width"],
            "values": [[1.0], [2.0], [3.0]],
        }
        out_path = tmp_path / "report.html"
        write_subset_html_report(
            ("RI", "biotype_filtered", "Low", "seqnames"), [task_result], out_path
        )
        html_text = out_path.read_text(encoding="utf-8")
        assert "Transformed feature distributions (after preprocessing)" in html_text
        assert "plotFeatureHistogram" in html_text

        payload = json.loads(out_path.with_suffix(".json").read_text(encoding="utf-8"))
        assert payload[0]["transformed_feature_distributions"] == {
            "num__width": [1.0, 2.0, 3.0]
        }

    def test_section_omitted_when_no_transformed_features(self, tmp_path):
        out_path = tmp_path / "report.html"
        write_subset_html_report(
            ("RI", "biotype_filtered", "Low", "seqnames"), [_make_task_result()], out_path
        )
        payload = json.loads(out_path.with_suffix(".json").read_text(encoding="utf-8"))
        assert payload[0]["transformed_feature_distributions"] == {}

    def test_plot_helper_function_exists(self):
        assert "function plotFeatureHistogram" in _PLOTLY_JS_HELPERS

    def test_two_pass_ordering_in_template(self, tmp_path):
        """Same CSS-grid/Plotly timing fix as the box-plot and confusion-matrix
        grids: all container divs created before any Plotly.newPlot call."""
        out_path = tmp_path / "report.html"
        write_subset_html_report(
            ("RI", "biotype_filtered", "Low", "seqnames"), [_make_task_result()], out_path
        )
        html_text = out_path.read_text(encoding="utf-8")

        div_creation_idx = html_text.index("const fdDivs = featureNames.map")
        first_plot_call_idx = html_text.index(
            "featureNames.forEach((name, i) => {"
        )
        assert div_creation_idx < first_plot_call_idx, (
            "container-div creation loop must appear before the "
            "Plotly.newPlot loop in the rendered template"
        )
