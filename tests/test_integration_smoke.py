"""Smoke integration tests for the splicing ML pipeline."""

import pytest
import tempfile
from pathlib import Path

from splicing_ml.pipeline import run_all_configurations
from splicing_ml.config import RunConfig


@pytest.mark.integration
class TestIntegrationSmoke:
    """Basic integration test to validate end-to-end execution."""

    def test_smoke_run_minimal(self):
        """Run a minimal pipeline on subset of data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            run_cfg = RunConfig(
                data_path="processed_data/aggregated_dt_filtered.validation300k.csv.gz",
                output_dir=tmpdir,
                only_event_type="SE",
                only_transcript_filter="transcripts",
                only_variability="High",
                outer_splits=2,
                inner_splits=2,
                include_models=("linear",),
                run_regression=True,
                run_classification=True,
                smoke_mode=True,
                smoke_max_rows=5000,
                verbose=False,
            )

            bundle = run_all_configurations(run_cfg)

            # Verify output structure.
            assert "results" in bundle, "Bundle missing results"
            assert "run_config" in bundle, "Bundle missing run_config"
            assert len(bundle["results"]) > 0, "No results generated"

            # Verify result structure for each config-task.
            for result in bundle["results"]:
                assert "status" in result, "Result missing status"
                if result["status"] == "ok":
                    assert "fold_results" in result, "Result missing fold_results"
                    assert len(result["fold_results"]) > 0, "No fold results"

                    # Verify fold structure.
                    for fold in result["fold_results"]:
                        assert "model_name" in fold, "Fold missing model_name"
                        assert "scores" in fold, "Fold missing scores"
                        assert "outer_fold" in fold, "Fold missing outer_fold"

    def test_smoke_run_compact_output(self):
        """Verify compact output mode reduces payload size."""
        with tempfile.TemporaryDirectory() as tmpdir:
            run_cfg = RunConfig(
                data_path="processed_data/aggregated_dt_filtered.validation300k.csv.gz",
                output_dir=tmpdir,
                only_event_type="SE",
                only_transcript_filter="transcripts",
                only_variability="High",
                outer_splits=2,
                inner_splits=2,
                include_models=("linear",),
                run_regression=True,
                run_classification=False,
                smoke_mode=True,
                smoke_max_rows=5000,
                output_level="compact",
                verbose=False,
            )

            bundle = run_all_configurations(run_cfg)

            # Verify compact output excludes predictions.
            for result in bundle["results"]:
                if result["status"] == "ok":
                    for fold in result["fold_results"]:
                        # Compact mode should NOT have y_true/y_pred.
                        assert (
                            "y_true" not in fold
                        ), "Compact mode should exclude y_true"
                        assert (
                            "y_pred" not in fold
                        ), "Compact mode should exclude y_pred"
                        # But should still have scores.
                        assert "scores" in fold, "Should have scores"

    def test_smoke_run_no_calibration(self):
        """Verify --no-calibration flag is respected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            run_cfg = RunConfig(
                data_path="processed_data/aggregated_dt_filtered.validation300k.csv.gz",
                output_dir=tmpdir,
                only_event_type="SE",
                only_transcript_filter="transcripts",
                only_variability="High",
                outer_splits=2,
                inner_splits=2,
                include_models=("linear",),
                run_regression=False,
                run_classification=True,
                smoke_mode=True,
                smoke_max_rows=5000,
                calibrate_classifiers=False,
                verbose=False,
            )

            bundle = run_all_configurations(run_cfg)

            # Should complete successfully even without calibration.
            assert len(bundle["results"]) > 0, "No results generated"
            for result in bundle["results"]:
                assert result["status"] in [
                    "ok",
                    "skipped",
                    "failed",
                ], f"Unexpected status: {result['status']}"

    def test_smoke_run_inner_splits_respected(self):
        """Verify inner_splits config is respected (not derived from outer_k)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            run_cfg = RunConfig(
                data_path="processed_data/aggregated_dt_filtered.validation300k.csv.gz",
                output_dir=tmpdir,
                only_event_type="SE",
                only_transcript_filter="transcripts",
                only_variability="High",
                outer_splits=3,
                inner_splits=4,  # Explicitly set to 4, not derived from outer_splits.
                include_models=("linear",),
                run_regression=True,
                run_classification=False,
                smoke_mode=True,
                smoke_max_rows=5000,
                verbose=True,  # Enable verbose to see inner split count in logs.
            )

            bundle = run_all_configurations(run_cfg)

            # Verify successful completion.
            assert len(bundle["results"]) > 0, "No results generated"
            for result in bundle["results"]:
                if result["status"] == "ok":
                    # Should have successfully tuned models, indicating inner CV worked.
                    for fold in result["fold_results"]:
                        assert "tuning" in fold, "Fold missing tuning info"
