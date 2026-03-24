"""Tests for grouped split logic in splicing_ml.modeling."""

import numpy as np
import pandas as pd
import pytest

from splicing_ml.modeling import (
    balanced_group_split_indices,
    bounded_group_splits,
)


class TestBalancedGroupSplits:
    """Test grouped split construction and balance."""

    def test_split_integrity_no_overlap(self):
        """Verify train/test splits have no overlap."""
        groups = pd.Series(["A", "A", "B", "B", "C", "C", "D", "D"] * 10)
        splits = balanced_group_split_indices(groups, n_splits=3, seed=42)

        for train_idx, test_idx in splits:
            overlap = set(train_idx) & set(test_idx)
            assert len(overlap) == 0, f"train/test overlap: {overlap}"

    def test_split_coverage(self):
        """Verify all observations are covered exactly once."""
        groups = pd.Series(["A", "A", "B", "B", "C", "C"] * 5)
        splits = balanced_group_split_indices(groups, n_splits=2, seed=42)

        all_train = set()
        all_test = set()
        for train_idx, test_idx in splits:
            all_train.update(train_idx)
            all_test.update(test_idx)

        total_coverage = all_train | all_test
        assert len(total_coverage) == len(
            groups
        ), f"coverage mismatch: {len(total_coverage)} vs {len(groups)}"

    def test_group_exclusivity(self):
        """Verify entire groups are assigned to either train or test, never both."""
        groups = pd.Series(["A"] * 5 + ["B"] * 5 + ["C"] * 5)
        splits = balanced_group_split_indices(groups, n_splits=2, seed=42)

        for train_idx, test_idx in splits:
            train_groups = set(groups.iloc[train_idx].unique())
            test_groups = set(groups.iloc[test_idx].unique())
            overlap = train_groups & test_groups
            assert (
                len(overlap) == 0
            ), f"group not exclusive: {overlap} appears in both train and test"

    def test_balanced_fold_sizes(self):
        """Verify folds have roughly balanced observation counts."""
        groups = pd.Series(["A"] * 100 + ["B"] * 80 + ["C"] * 60)
        splits = balanced_group_split_indices(groups, n_splits=3, seed=42)

        fold_sizes = [len(te_idx) for _, te_idx in splits]
        mean_size = np.mean(fold_sizes)
        max_imbalance = max(abs(s - mean_size) / mean_size for s in fold_sizes)
        # Allow up to 50% imbalance due to group-boundary constraints.
        assert (
            max_imbalance < 0.5
        ), f"fold size imbalance too high: {fold_sizes}, mean={mean_size}"

    def test_at_least_two_splits(self):
        """Verify at least 2 splits are generated."""
        groups = pd.Series(["A"] * 10 + ["B"] * 10)
        splits = balanced_group_split_indices(groups, n_splits=5, seed=42)
        assert len(splits) >= 2, "fewer than 2 splits generated"

    def test_single_group_raises(self):
        """Verify error when only one group is available."""
        groups = pd.Series(["A"] * 20)
        with pytest.raises(ValueError, match="at least 2 unique groups"):
            balanced_group_split_indices(groups, n_splits=2, seed=42)


class TestBoundedGroupSplits:
    """Test split count bounding logic."""

    def test_clamp_to_max(self):
        """Verify split count is clamped to max_splits."""
        result = bounded_group_splits(
            requested_splits=20,
            n_groups=15,
            max_splits=10,
        )
        assert result == 10, f"not clamped to max: {result}"

    def test_clamp_to_groups(self):
        """Verify split count does not exceed number of groups."""
        result = bounded_group_splits(
            requested_splits=10,
            n_groups=3,
            max_splits=10,
        )
        assert result == 3, f"exceeds groups: {result}"

    def test_minimum_two(self):
        """Verify minimum of 2 splits."""
        result = bounded_group_splits(
            requested_splits=1,
            n_groups=10,
            max_splits=10,
        )
        assert result == 2, f"below minimum: {result}"

    def test_respects_valid_request(self):
        """Verify valid requests pass through unmodified."""
        result = bounded_group_splits(
            requested_splits=5,
            n_groups=10,
            max_splits=10,
        )
        assert result == 5, f"modified valid request: {result}"
