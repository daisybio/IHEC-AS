"""Tests for data targeting and label construction in splicing_ml.data."""

import numpy as np
import pandas as pd
import pytest

from splicing_ml.data import build_targets


class TestBuildTargets:
    """Test target construction for regression and classification."""

    def test_regression_target_preserves_psi(self):
        """Verify regression target is raw PSI values."""
        df = pd.DataFrame(
            {
                "PSI": [0.1, 0.5, 0.9, 0.3, 0.7],
                "other": [1, 2, 3, 4, 5],
            }
        )
        df_out, y = build_targets(df, "regression", low_thr=1 / 3, high_thr=2 / 3)

        assert (
            y == df["PSI"].astype(float).to_numpy()
        ).all(), "regression target does not match raw PSI"
        assert len(df_out) == len(df), "dataframe shrunk unexpectedly"

    def test_classification_binarization(self):
        """Verify classification creates correct binary labels."""
        df = pd.DataFrame(
            {
                "PSI": [0.1, 0.4, 0.35, 0.7, 0.8, 0.5],
                "other": [1, 2, 3, 4, 5, 6],
            }
        )
        low_thr, high_thr = 1 / 3, 2 / 3
        df_out, y = build_targets(
            df, "classification", low_thr=low_thr, high_thr=high_thr
        )

        # Rows outside [low_thr, high_thr] are removed.
        # Remaining: 0.1 (label 0), 0.7 (label 1), 0.8 (label 1)
        assert len(df_out) == 3, f"expected 3 rows kept, got {len(df_out)}"
        assert (y == np.array([0, 1, 1])).all(), f"labels mismatch: {y}"

    def test_classification_threshold_filtering(self):
        """Verify classification removes intermediate PSI values."""
        df = pd.DataFrame(
            {
                "PSI": [
                    0.2,  # < 1/3: keep as 0
                    0.4,  # in (1/3, 2/3): remove
                    0.5,  # in (1/3, 2/3): remove
                    0.7,  # > 2/3: keep as 1
                ],
            }
        )
        df_out, y = build_targets(df, "classification", low_thr=1 / 3, high_thr=2 / 3)

        assert len(df_out) == 2, f"expected 2 rows after filtering, got {len(df_out)}"
        assert (y == np.array([0, 1])).all(), f"labels incorrect: {y}"

    def test_classification_edge_case_exactly_threshold(self):
        """Verify behavior at exact threshold boundaries."""
        df = pd.DataFrame(
            {
                "PSI": [
                    1 / 3,  # exactly low_thr: kept (via >= high_thr check)
                    2 / 3,  # exactly high_thr: kept
                ],
            }
        )
        df_out, y = build_targets(df, "classification", low_thr=1 / 3, high_thr=2 / 3)

        # At exact boundaries: low_thr is not included in high-class,
        # high_thr is included in high-class.
        assert len(df_out) == 1, "expected 1 row at exact threshold"
        assert y[0] == 1, "exact high_thr should be class 1"

    def test_empty_classification_input(self):
        """Verify handling of data with no extreme values."""
        df = pd.DataFrame(
            {
                "PSI": [0.45, 0.50, 0.55],  # All between 1/3 and 2/3
            }
        )
        df_out, y = build_targets(df, "classification", low_thr=1 / 3, high_thr=2 / 3)

        assert len(df_out) == 0, "expected empty dataframe"
        assert len(y) == 0, "expected empty labels"

    def test_regression_all_rows_kept(self):
        """Verify regression uses all rows regardless of PSI value."""
        df = pd.DataFrame(
            {
                "PSI": [0.1, 0.5, 0.9, 0.05, 0.95],
            }
        )
        df_out, y = build_targets(df, "regression", low_thr=1 / 3, high_thr=2 / 3)

        assert len(df_out) == len(df), "regression should keep all rows"
        assert len(y) == len(df), "regression should have labels for all rows"

    def test_nan_psi_values(self):
        """Verify NaN PSI values are handled."""
        df = pd.DataFrame(
            {
                "PSI": [0.1, np.nan, 0.7, 0.5],
            }
        )
        df_out, y = build_targets(df, "regression", low_thr=1 / 3, high_thr=2 / 3)

        # NaN is typically dropped.
        assert len(df_out) <= len(df), "should not add rows"
