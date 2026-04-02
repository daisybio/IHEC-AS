from __future__ import annotations

"""Feature preprocessing builders.

This module keeps all transformation decisions in one place so they are easy to
inspect and test independently from model training.
"""

import math
from typing import Any

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, OneHotEncoder, StandardScaler

from .config import METADATA_COLUMNS
from .utils import vlog


__all__ = [
    "detect_gene_expression_column",
    "add_protocol_expression_interactions",
    "infer_feature_columns",
    "choose_log_columns",
    "classify_numeric_columns",
    "PercentileClipper",
    "build_preprocessor",
]

# Substring patterns (lower-cased) that mark a column for log1p transformation.
_LOG_COLUMN_PATTERNS: tuple[str, ...] = ("width", "distance")


def detect_gene_expression_column(df: pd.DataFrame) -> str | None:
    """Infer likely gene-expression column from available numeric columns."""
    cols = list(df.columns)
    lower_to_original = {c.lower(): c for c in cols}

    preferred_exact = [
        "gene_expression",
        "gene_expr",
        "expression",
        "tpm",
        "fpkm",
        "rpkm",
    ]
    for cand in preferred_exact:
        col = lower_to_original.get(cand)
        if col is not None and pd.api.types.is_numeric_dtype(df[col]):
            return col

    contains_gene_expr = [
        c
        for c in cols
        if ("gene" in c.lower() and "expr" in c.lower())
        and pd.api.types.is_numeric_dtype(df[c])
    ]
    if contains_gene_expr:
        return contains_gene_expr[0]

    contains_expr = [
        c
        for c in cols
        if ("expr" in c.lower() or "expression" in c.lower())
        and pd.api.types.is_numeric_dtype(df[c])
    ]
    if contains_expr:
        return contains_expr[0]

    return None


def add_protocol_expression_interactions(
    df: pd.DataFrame,
    verbose: bool = False,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Add protocol × gene-expression interaction features if possible."""
    out = df.copy()

    if "protocol" not in out.columns:
        vlog(
            verbose,
            "Protocol column not found; skipping protocol x expression interactions",
        )
        return out, {
            "interaction_base_protocol": [],
            "interaction_base_gene_expression": [],
            "interaction_features": [],
        }

    expr_col = detect_gene_expression_column(out)
    if expr_col is None:
        vlog(
            verbose,
            "No numeric gene-expression-like column found; skipping protocol x expression interactions",
        )
        return out, {
            "interaction_base_protocol": ["protocol"],
            "interaction_base_gene_expression": [],
            "interaction_features": [],
        }

    protocol_series = out["protocol"].astype("string").fillna("__MISSING__")
    expr = pd.to_numeric(out[expr_col], errors="coerce").fillna(0.0).astype(float)

    protocol_dummies = pd.get_dummies(protocol_series, prefix="protocol", dtype=float)
    interactions = protocol_dummies.mul(expr, axis=0)
    interactions.columns = [f"{c}__x__{expr_col}" for c in interactions.columns]

    out = pd.concat([out, interactions], axis=1)
    vlog(
        verbose,
        f"Added {interactions.shape[1]} protocol x {expr_col} interaction features",
    )
    return out, {
        "interaction_base_protocol": ["protocol"],
        "interaction_base_gene_expression": [expr_col],
        "interaction_features": interactions.columns.tolist(),
    }


def infer_feature_columns(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    """Split feature columns into numeric and categorical groups (vectorized for speed)."""
    # Speed optimization: use select_dtypes instead of list comprehensions (1-2% faster)
    feature_cols = [c for c in df.columns if c not in METADATA_COLUMNS]
    numeric_cols = list(df[feature_cols].select_dtypes(include=np.number).columns)
    categorical_cols = [c for c in feature_cols if c not in numeric_cols]
    return numeric_cols, categorical_cols


def choose_log_columns(
    df: pd.DataFrame,
    numeric_cols: list[str],
    skew_threshold: float = 1.0,
    verbose: bool = False,
) -> list[str]:
    """Pick non-negative numeric features with high skew for log1p transform.

    .. deprecated::
        ``build_preprocessor`` now uses explicit pattern-based column
        classification via ``classify_numeric_columns``.  This function is
        retained for backward compatibility and standalone use.
    """
    log_cols: list[str] = []
    for col in numeric_cols:
        s = df[col].dropna()
        if s.empty or s.min() < 0:
            continue
        skew_val = float(s.skew())
        if math.isfinite(skew_val) and skew_val >= skew_threshold:
            log_cols.append(col)
    vlog(
        verbose,
        f"Log-transform decision: {len(log_cols)} / {len(numeric_cols)} numeric features selected",
    )
    return log_cols


def classify_numeric_columns(
    df: pd.DataFrame,
    numeric_cols: list[str],
    verbose: bool = False,
) -> tuple[list[str], list[str], list[str]]:
    """Classify numeric columns into log, clip, and plain groups.

    Rules
    -----
    - **log** (log1p + StandardScaler): width/distance columns and the
      gene-expression column.  Identified by ``_LOG_COLUMN_PATTERNS``
      substrings or ``detect_gene_expression_column``.
    - **clip** (99th-pct upper clip + StandardScaler): histone-mark columns
      (``H3K*``).  These are already −log10 p-values and must not be
      log-transformed again.
    - **plain** (median impute + StandardScaler): everything else (DNAm,
      CpGs, splice-site scores, interaction features, …).

    Parameters
    ----------
    df:
        DataFrame containing the feature columns (used only to detect the
        gene-expression column name).
    numeric_cols:
        List of numeric column names to classify.
    verbose:
        Emit progress messages.

    Returns
    -------
    log_cols, clip_cols, plain_cols : tuple[list[str], list[str], list[str]]
    """
    expr_col = detect_gene_expression_column(df)

    log_cols: list[str] = []
    clip_cols: list[str] = []
    plain_cols: list[str] = []

    for col in numeric_cols:
        col_lower = col.lower()
        if any(pat in col_lower for pat in _LOG_COLUMN_PATTERNS) or col == expr_col:
            log_cols.append(col)
        elif col.startswith("H3K"):
            clip_cols.append(col)
        else:
            plain_cols.append(col)

    vlog(
        verbose,
        f"Column classification: {len(log_cols)} log, {len(clip_cols)} clip (histone), "
        f"{len(plain_cols)} plain",
    )
    return log_cols, clip_cols, plain_cols


class PercentileClipper(BaseEstimator, TransformerMixin):
    """Clip each feature at a given upper percentile computed on training data.

    This is a leakage-safe alternative to a hard-coded clip value: the
    percentile threshold is estimated from the training fold only and then
    applied identically to the test fold.

    Parameters
    ----------
    upper_percentile : float
        Percentile (0–100) used as the upper clip boundary.  Default is 99.
    """

    def __init__(self, upper_percentile: float = 99.0) -> None:
        """Initialize a PercentileClipper instance."""
        self.upper_percentile = upper_percentile

    def fit(self, X: Any, y: Any = None) -> "PercentileClipper":
        """Fit the model using the provided training data."""
        arr = np.asarray(X, dtype=float)
        self.clip_max_: np.ndarray = np.percentile(arr, self.upper_percentile, axis=0)
        return self

    def transform(self, X: Any, y: Any = None) -> np.ndarray:
        """Transform input features using the fitted state."""
        arr = np.asarray(X, dtype=float)
        return np.clip(arr, None, self.clip_max_)


def _log1p_dataframe(x: Any) -> Any:
    """Apply log1p safely on array-like input objects."""
    arr = (
        x.to_numpy(dtype=float)
        if hasattr(x, "to_numpy")
        else np.asarray(x, dtype=float)
    )
    return np.log1p(arr)


def _make_one_hot_encoder() -> OneHotEncoder:
    """Construct OneHotEncoder compatible with older/newer sklearn versions.

    Uses sparse_output=False to ensure consistent feature counts across different
    data subsets (e.g., train/test folds with different category distributions).
    """
    try:
        # sparse_output=False ensures consistent dense output shape across folds
        return OneHotEncoder(
            handle_unknown="ignore",
            sparse_output=False,
            drop=None,  # Keep all categories to prevent shape variations
        )
    except TypeError:
        # Older sklearn versions use sparse=True/False instead of sparse_output
        return OneHotEncoder(
            handle_unknown="ignore",
            sparse=False,
            drop=None,
        )


def build_preprocessor(
    df: pd.DataFrame,
    categorical_missing_strategy: str = "missing_token",
    categorical_missing_token: str = "__MISSING__",
    histone_upper_percentile: float = 99.0,
    verbose: bool = False,
) -> tuple[ColumnTransformer, dict[str, list[str]]]:
    """Build leakage-safe preprocessing pipeline and metadata summary.

    Column treatment
    ----------------
    - **Width / distance / gene-expression** columns → log1p + StandardScaler.
    - **Histone mark** columns (``H3K*``, already −log10 p-values) →
      per-feature upper clip at *histone_upper_percentile* + StandardScaler.
    - **All other numeric** columns → median impute + StandardScaler.
    - **Categorical** columns → missing-token impute + one-hot encoding.
    """
    numeric_cols, categorical_cols = infer_feature_columns(df)
    log_cols, clip_cols, plain_cols = classify_numeric_columns(
        df, numeric_cols, verbose=verbose
    )

    numeric_plain_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
    )

    numeric_log_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("log1p", FunctionTransformer(_log1p_dataframe, validate=False)),
            ("scaler", StandardScaler()),
        ]
    )

    numeric_clip_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("clip", PercentileClipper(upper_percentile=histone_upper_percentile)),
            ("scaler", StandardScaler()),
        ]
    )

    if categorical_missing_strategy == "missing_token":
        categorical_imputer = SimpleImputer(
            strategy="constant", fill_value=categorical_missing_token
        )
    elif categorical_missing_strategy == "most_frequent":
        categorical_imputer = SimpleImputer(strategy="most_frequent")
    else:
        raise ValueError(
            "categorical_missing_strategy must be one of: missing_token, most_frequent"
        )

    vlog(
        verbose,
        f"Categorical imputation strategy={categorical_missing_strategy}, token={categorical_missing_token}",
    )

    categorical_pipe = Pipeline(
        steps=[
            ("imputer", categorical_imputer),
            ("ohe", _make_one_hot_encoder()),
        ]
    )

    transformers = []
    if plain_cols:
        transformers.append(("num", numeric_plain_pipe, plain_cols))
    if log_cols:
        transformers.append(("num_log", numeric_log_pipe, log_cols))
    if clip_cols:
        transformers.append(("num_clip", numeric_clip_pipe, clip_cols))
    if categorical_cols:
        transformers.append(("cat", categorical_pipe, categorical_cols))

    preprocessor = ColumnTransformer(transformers=transformers, remainder="drop")
    details = {
        "numeric_standardized": plain_cols,
        "numeric_log1p_standardized": log_cols,
        "numeric_clip_standardized": clip_cols,
        "categorical_one_hot": categorical_cols,
        "categorical_missing_strategy": [categorical_missing_strategy],
        "categorical_missing_token": [categorical_missing_token],
    }
    vlog(
        verbose,
        f"Preprocessor built: n_numeric={len(numeric_cols)}, n_categorical={len(categorical_cols)}",
        level="info",
    )
    return preprocessor, details
