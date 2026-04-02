"""Tests for metric utilities in splicing_ml.metrics."""

import numpy as np
import pytest

from splicing_ml.metrics import (
    bootstrap_ci,
    classification_metrics,
    regression_metrics,
    regression_metrics_by_psi_bin,
    r2_from_rss,
    concordance_correlation_coefficient,
)


class TestRegressionMetrics:
    """Test regression metric computation."""

    def test_rmse_basic(self):
        """Test RMSE calculation."""
        y_true = np.array([1.0, 2.0, 3.0])
        y_pred = np.array([1.1, 2.1, 2.9])
        metrics = regression_metrics(y_true, y_pred)

        expected_rmse = np.sqrt(((0.1**2 + 0.1**2 + 0.1**2) / 3))
        assert np.isclose(
            metrics["rmse"], expected_rmse
        ), f"RMSE mismatch: {metrics['rmse']} vs {expected_rmse}"

    def test_r2_rss(self):
        """Test R² computed from residual sum of squares."""
        y_true = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y_pred = np.array([1.1, 2.0, 2.9, 4.1, 4.9])
        metrics = regression_metrics(y_true, y_pred)

        # R² = 1 - RSS/TSS
        rss = np.sum((y_true - y_pred) ** 2)
        tss = np.sum((y_true - np.mean(y_true)) ** 2)
        expected_r2 = 1.0 - (rss / tss)

        assert np.isclose(
            metrics["r2_rss"], expected_r2
        ), f"R² mismatch: {metrics['r2_rss']} vs {expected_r2}"

    def test_mad(self):
        """Test mean absolute deviation."""
        y_true = np.array([1.0, 2.0, 3.0])
        y_pred = np.array([1.2, 2.0, 2.7])
        metrics = regression_metrics(y_true, y_pred)

        expected_mad = np.mean(np.abs([0.2, 0.0, 0.3]))
        assert np.isclose(
            metrics["mad"], expected_mad
        ), f"MAD mismatch: {metrics['mad']} vs {expected_mad}"

    def test_perfect_prediction(self):
        """Test metrics when prediction is perfect."""
        y = np.array([1.0, 2.0, 3.0, 4.0])
        metrics = regression_metrics(y, y)

        assert metrics["rmse"] == 0.0, "RMSE should be 0 for perfect prediction"
        assert metrics["mad"] == 0.0, "MAD should be 0 for perfect prediction"
        assert metrics["r2_rss"] == 1.0, "R² should be 1.0 for perfect prediction"

    def test_constant_prediction(self):
        """Test metrics when prediction is constant (constant baseline)."""
        y_true = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y_pred = np.array([3.0, 3.0, 3.0, 3.0, 3.0])  # mean
        metrics = regression_metrics(y_true, y_pred)

        assert metrics["r2_rss"] == 0.0, "R² should be 0 for constant mean prediction"


class TestRegressionMetricsByPsiBin:
    """Tests for per-PSI-bin regression metrics."""

    def test_bin_keys_match_range(self):
        """Test bin keys match range."""
        rng = np.random.default_rng(0)
        y = rng.uniform(0.2, 0.8, 300)
        result = regression_metrics_by_psi_bin(y, y, psi_low=0.2, psi_high=0.8)
        assert len(result) == 3
        assert "bin_0.20_0.40" in result
        assert "bin_0.40_0.60" in result
        assert "bin_0.60_0.80" in result

    def test_perfect_prediction_gives_zero_rmse(self):
        """Test perfect prediction gives zero rmse."""
        rng = np.random.default_rng(1)
        y = rng.uniform(0.2, 0.8, 300)
        result = regression_metrics_by_psi_bin(y, y, psi_low=0.2, psi_high=0.8)
        for bin_metrics in result.values():
            assert bin_metrics["rmse"] == pytest.approx(0.0, abs=1e-10)

    def test_sparse_bin_returns_nan(self):
        # Only one sample in [0.2, 0.4)
        """Test sparse bin returns nan."""
        y_true = np.array([0.25, 0.5, 0.5, 0.5, 0.7])
        y_pred = np.array([0.30, 0.5, 0.5, 0.5, 0.7])
        result = regression_metrics_by_psi_bin(y_true, y_pred, psi_low=0.2, psi_high=0.8)
        assert result["bin_0.20_0.40"]["n"] == 1
        assert np.isnan(result["bin_0.20_0.40"]["rmse"])

    def test_n_counts_sum_to_total(self):
        """Test n counts sum to total."""
        rng = np.random.default_rng(2)
        y = rng.uniform(0.2, 0.8, 200)
        result = regression_metrics_by_psi_bin(y, y, psi_low=0.2, psi_high=0.8)
        total_n = sum(v["n"] for v in result.values())
        assert total_n == len(y)

    def test_custom_psi_range(self):
        """Test custom psi range."""
        rng = np.random.default_rng(3)
        y = rng.uniform(0.1, 0.9, 300)
        result = regression_metrics_by_psi_bin(y, y, psi_low=0.1, psi_high=0.9, n_bins=4)
        assert len(result) == 4
        assert "bin_0.10_0.30" in result


class TestClassificationMetrics:
    """Test classification metric computation."""

    def test_balanced_accuracy_perfect(self):
        """Test balanced accuracy when perfect."""
        y_true = np.array([0, 0, 1, 1, 0, 1])
        y_prob = np.array([0.1, 0.2, 0.8, 0.9, 0.1, 0.9])
        threshold = 0.5

        metrics = classification_metrics(y_true, y_prob, threshold)
        assert (
            metrics["balanced_accuracy"] == 1.0
        ), "balanced_accuracy should be 1.0 for perfect predictions"

    def test_balanced_accuracy_random_baseline(self):
        """Test balanced accuracy for random predictions."""
        # Random binary predictions should give ~0.5 balanced accuracy.
        np.random.seed(42)
        y_true = np.random.randint(0, 2, size=100)
        y_prob = np.random.rand(100)
        threshold = 0.5

        metrics = classification_metrics(y_true, y_prob, threshold)
        # Should be close to 0.5 but allow some variance.
        assert (
            0.3 < metrics["balanced_accuracy"] < 0.7
        ), f"random baseline should be ~0.5, got {metrics['balanced_accuracy']}"

    def test_auroc_requires_both_classes(self):
        """Test AUROC returns NaN when only one class present."""
        y_true = np.array([1, 1, 1, 1])
        y_prob = np.array([0.9, 0.8, 0.7, 0.6])
        threshold = 0.5

        metrics = classification_metrics(y_true, y_prob, threshold)
        assert np.isnan(
            metrics["auroc"]
        ), "AUROC should be NaN when only one class present"

    def test_auprc_single_class(self):
        """Test AUPRC with single positive class."""
        y_true = np.array([1, 1, 1, 1])
        y_prob = np.array([0.9, 0.8, 0.7, 0.6])
        threshold = 0.5

        metrics = classification_metrics(y_true, y_prob, threshold)
        # When all labels are 1, AUPRC should be 1.0 (all predictions correct).
        assert (
            metrics["auprc"] == 1.0
        ), "AUPRC should be 1.0 when all labels are positive"

    def test_f1_zero_division(self):
        """Test F1 score with no positive predictions."""
        y_true = np.array([1, 1, 0, 0])
        y_prob = np.array([0.1, 0.2, 0.9, 0.8])  # predictions: [0, 0, 1, 1]
        threshold = 0.5

        metrics = classification_metrics(y_true, y_prob, threshold)
        # Precision undefined (no positive predictions), F1 should handle gracefully.
        assert isinstance(metrics["f1"], float), "F1 should be a float"
        assert 0.0 <= metrics["f1"] <= 1.0, "F1 should be in [0, 1]"

    def test_mcc_perfect_prediction(self):
        """Test MCC for perfect prediction."""
        y_true = np.array([0, 0, 1, 1])
        y_prob = np.array([0.1, 0.2, 0.8, 0.9])
        threshold = 0.5

        metrics = classification_metrics(y_true, y_prob, threshold)
        assert metrics["mcc"] == 1.0, "MCC should be 1.0 for perfect prediction"


class TestBootstrapCI:
    """Test bootstrap confidence interval computation."""

    def test_bootstrap_ci_mean_in_interval(self):
        """Verify sample mean is typically within bootstrap CI."""
        np.random.seed(42)
        values = np.random.normal(loc=5.0, scale=1.0, size=1000)
        low, high = bootstrap_ci(values, n_bootstrap=100, seed=42)

        sample_mean = np.mean(values)
        assert (
            low <= sample_mean <= high
        ), f"sample mean {sample_mean} outside CI [{low}, {high}]"

    def test_bootstrap_ci_width(self):
        """Verify bootstrap CI width is reasonable."""
        np.random.seed(42)
        values = np.random.normal(loc=0.0, scale=1.0, size=500)
        low, high = bootstrap_ci(values, n_bootstrap=100, seed=42)

        width = high - low
        # For N≈500, std≈1, 95% CI width should be ~0.1-0.2.
        assert 0.05 < width < 0.5, f"CI width unreasonable: {width}"

    def test_bootstrap_ci_empty_input(self):
        """Test CI for empty input."""
        values = np.array([])
        low, high = bootstrap_ci(values, n_bootstrap=100, seed=42)

        assert np.isnan(low) and np.isnan(high), "CI should be NaN for empty input"

    def test_bootstrap_ci_single_value(self):
        """Test CI for single value."""
        values = np.array([5.0])
        low, high = bootstrap_ci(values, n_bootstrap=100, seed=42)

        # Single value resampled always gives that value.
        assert low == 5.0 and high == 5.0, "CI should be [5.0, 5.0] for single value"
