"""Tests for splicing_ml.models.evaluator._compute_transformed_feature_sample."""

import numpy as np
import pandas as pd
import pytest
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from splicing_ml.models.evaluator import _compute_transformed_feature_sample


def _make_prep_pipeline(n=10):
    """A tiny fitted Pipeline with a "prep" step mirroring build_preprocessor's
    naming convention: a plain numeric transformer and a "cat" one-hot step."""
    df = pd.DataFrame(
        {
            "num1": np.arange(n, dtype=float),
            "cat1": (["a", "b"] * ((n // 2) + 1))[:n],
        }
    )
    prep = ColumnTransformer(
        transformers=[
            ("num", StandardScaler(), ["num1"]),
            ("cat", OneHotEncoder(sparse_output=False), ["cat1"]),
        ]
    )
    prep.set_output(transform="pandas")
    pipe = Pipeline(steps=[("prep", prep)])
    pipe.fit(df)
    return pipe, df


class TestComputeTransformedFeatureSample:
    def test_excludes_one_hot_columns_and_keeps_numeric(self):
        pipe, df = _make_prep_pipeline()
        result = _compute_transformed_feature_sample(
            pipe, df, verbose=False, max_samples=0
        )
        assert result is not None
        assert all(not f.startswith("cat__") for f in result["feature_names"])
        assert "num__num1" in result["feature_names"]
        assert result["values"].shape == (len(df), len(result["feature_names"]))

    def test_subsamples_when_max_samples_below_row_count(self):
        pipe, df = _make_prep_pipeline(n=20)
        result = _compute_transformed_feature_sample(
            pipe, df, verbose=False, max_samples=5
        )
        assert result["values"].shape[0] == 5

    def test_no_subsampling_when_max_samples_le_zero(self):
        pipe, df = _make_prep_pipeline(n=20)
        result = _compute_transformed_feature_sample(
            pipe, df, verbose=False, max_samples=0
        )
        assert result["values"].shape[0] == 20

    def test_no_subsampling_when_row_count_below_max_samples(self):
        pipe, df = _make_prep_pipeline(n=5)
        result = _compute_transformed_feature_sample(
            pipe, df, verbose=False, max_samples=5000
        )
        assert result["values"].shape[0] == 5

    def test_returns_none_without_prep_step(self):
        pipe = Pipeline(steps=[("scaler", StandardScaler())])
        pipe.fit(pd.DataFrame({"a": [1.0, 2.0]}))
        result = _compute_transformed_feature_sample(
            pipe, pd.DataFrame({"a": [1.0, 2.0]}), verbose=False
        )
        assert result is None

    def test_returns_none_for_non_pipeline_estimator(self):
        result = _compute_transformed_feature_sample(
            StandardScaler(), pd.DataFrame({"a": [1.0, 2.0]}), verbose=False
        )
        assert result is None

    def test_returns_none_if_only_categorical_columns_survive(self):
        df = pd.DataFrame({"cat1": ["a", "b", "a", "b"]})
        prep = ColumnTransformer(
            transformers=[("cat", OneHotEncoder(sparse_output=False), ["cat1"])]
        )
        prep.set_output(transform="pandas")
        pipe = Pipeline(steps=[("prep", prep)])
        pipe.fit(df)
        result = _compute_transformed_feature_sample(pipe, df, verbose=False)
        assert result is None
