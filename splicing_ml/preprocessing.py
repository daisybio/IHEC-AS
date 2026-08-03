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
    "FEATURE_GROUPS",
    "select_feature_group_columns",
]

# Substring patterns (lower-cased) that mark a column for log1p transformation.
_LOG_COLUMN_PATTERNS: tuple[str, ...] = ("width", "distance")

# Ablation-study feature groups (§ splicing_ml ablations, 2026-07-14): every
# feature column must fall into exactly one group so `select_feature_group_columns`
# can silently-drop nothing by accident. Regexes match the literal column
# names in aggregated_dt_filtered_{tf}.csv.gz (see 05-create-aggregated-dt.Rmd).
FEATURE_GROUPS: dict[str, str] = {
    # H3K*;{3down,3up,5down,5up} + H3K*_source (observed/imputed flag per mark)
    "histone": r"^H3K",
    # WGBS methylation + CpG coverage windows
    "dnam": r"^(DNAm|CpGs);",
    # Sequence-only splice-site strength: MaxEntScan (3ss/3ssdown/5ss/5ssup) +
    # pangolin (sequence+tissue CNN, no per-sample RNA-seq data used) + GC content
    # of the four 200bp splice-site-flanking windows (gc_5up/gc_5down/gc_3up/gc_3down,
    # added by 03-prepare-aggregation.Rmd's §3.7 block). GC belongs here for the same
    # reason pangolin does: derived from genome sequence alone, no per-sample data.
    "sequence": r"^(pangolin_|[35]ss|gc_[35])",
    # RNA-binding-protein eCLIP binding-site aggregates
    "rbp": r"^rbp_",
    # Core spliceosome component expression
    "spliceosome": r"^spliceosome_",
    # Host gene expression (getmm/vst)
    "gene_expression": r"^gene_expression_",
    # Genomic position / structural context (not a biological signal per se)
    "positional": r"^(distance_|width;)",
    # Confounds / context: protocol, ontology, project, their interaction
    # term, and qc_flag_count (currently all-NaN pipeline-wide, see memory
    # project_qc_flag_count_all_nan.md -- kept as its own group so ablating it
    # is a no-op today rather than silently missing from every group).
    "confound": r"^(protocol|ontology|project|qc_flag_count|protocol_x_)",
}


def select_feature_group_columns(
    columns: list[str],
    groups: tuple[str, ...] | None,
    verbose: bool = False,
) -> list[str]:
    """Restrict `columns` to the requested FEATURE_GROUPS, for ablation runs.

    `groups=None` (or containing "all") returns `columns` unchanged. Any
    column matching none of FEATURE_GROUPS's patterns is treated as an error
    (not silently kept or dropped) -- ablation results are meaningless if an
    untracked column leaks through unfiltered.
    """
    import re

    if groups is None or "all" in groups:
        return list(columns)

    unknown = set(groups) - set(FEATURE_GROUPS)
    if unknown:
        raise ValueError(
            f"Unknown feature group(s) {sorted(unknown)}; "
            f"valid groups: {sorted(FEATURE_GROUPS)} (or 'all')"
        )

    patterns = [re.compile(FEATURE_GROUPS[g]) for g in groups]
    kept = [c for c in columns if any(p.match(c) for p in patterns)]

    untagged = [
        c for c in columns if not any(re.match(p, c) for p in FEATURE_GROUPS.values())
    ]
    if untagged:
        raise ValueError(
            f"Column(s) {untagged} match no FEATURE_GROUPS pattern -- update "
            "FEATURE_GROUPS before running a feature-group ablation, or "
            "results silently miscount what's included."
        )

    vlog(
        verbose,
        f"Feature-group ablation: kept {len(kept)}/{len(columns)} columns for groups={groups}",
        level="info",
    )
    return kept


def detect_gene_expression_column(df: pd.DataFrame) -> str | None:
    """Infer likely gene-expression column from available numeric columns.

    §4.12b routing: `05-create-aggregated-dt.Rmd` emits BOTH
    `gene_expression_getmm` (cross-gene comparable, pooled route — what
    splicing_ml uses) and `gene_expression_vst` (event-specific route, used by
    06/07/09-1 instead). Both match the old fuzzy `contains_gene_expr` fallback
    below (("gene" in c) and ("expr" in c)) and it just returned whichever
    happened to sort first in `df.columns` — non-deterministic w.r.t. which of
    the two is picked, and could silently pick vst (wrong for the pooled route:
    vst can be negative, breaking the log1p classification in
    `classify_numeric_columns`). Checking `gene_expression_getmm` exactly,
    before the fuzzy fallback, makes this deterministic.
    """
    cols = list(df.columns)
    lower_to_original = {c.lower(): c for c in cols}

    preferred_exact = [
        "gene_expression_getmm",
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

    # drop_first: computed once globally over the whole dataset (before any CV
    # split), so unlike the per-fold OneHotEncoder in _make_one_hot_encoder(),
    # there's no fold-schema-drift risk here — only the collinearity of
    # building both protocol dummies (they sum to 1), which fed straight into
    # both interaction columns and made their coefficients non-identifiable.
    protocol_dummies = pd.get_dummies(
        protocol_series, prefix="protocol", dtype=float, drop_first=True
    )
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
        self.n_features_in_ = arr.shape[1]
        if hasattr(X, "columns"):
            self.feature_names_in_ = np.asarray(X.columns, dtype=object)
        return self

    def transform(self, X: Any, y: Any = None) -> np.ndarray:
        """Transform input features using the fitted state."""
        arr = np.asarray(X, dtype=float)
        return np.clip(arr, None, self.clip_max_)

    def get_feature_names_out(self, input_features: Any = None) -> np.ndarray:
        """Return feature names unchanged (clip is a 1:1 transform)."""
        if input_features is not None:
            return np.asarray(input_features, dtype=object)
        if hasattr(self, "feature_names_in_"):
            return np.asarray(self.feature_names_in_, dtype=object)
        return np.asarray(
            [f"x{i}" for i in range(self.n_features_in_)], dtype=object
        )


def _log1p_dataframe(x: Any) -> Any:
    """Apply log1p safely on array-like input objects."""
    arr = (
        x.to_numpy(dtype=float)
        if hasattr(x, "to_numpy")
        else np.asarray(x, dtype=float)
    )
    return np.log1p(arr)


def _make_one_hot_encoder(categories: list[np.ndarray] | str = "auto") -> OneHotEncoder:
    """Construct OneHotEncoder compatible with older/newer sklearn versions.

    Uses sparse_output=False to ensure consistent feature counts across different
    data subsets (e.g., train/test folds with different category distributions).

    drop="first" (§2.5b fix): categories_ is fixed at fit time and
    handle_unknown="ignore" already guarantees fold-shape consistency on its
    own (unseen categories at transform time become all-zero rows, not new
    columns) — drop has no bearing on that. Without it, a 2-level column like
    `protocol` produced BOTH dummy columns, which are perfectly anti-correlated
    (r=1.00 co-occurrence), inflating/splitting feature importance between them.

    categories (fix, 2026-07-16): passing the explicit global per-column
    category lists (computed by build_preprocessor from the FULL dataset,
    before outer-fold splitting) means every fold's encoder already knows
    every real category up front, so grouped CV (ontology/seqnames) folds
    where a category is absent from that fold's train split no longer hit
    "unknown category" at transform time — a real fix (every real category
    gets a proper one-hot column instead of being silently zeroed out), not
    a suppression of the sklearn UserWarning that used to fire there.
    """
    try:
        # sparse_output=False ensures consistent dense output shape across folds
        return OneHotEncoder(
            categories=categories,
            handle_unknown="ignore",
            sparse_output=False,
            drop="first",
        )
    except TypeError:
        # Older sklearn versions use sparse=True/False instead of sparse_output
        return OneHotEncoder(
            categories=categories,
            handle_unknown="ignore",
            sparse=False,
            drop="first",
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
            (
                "log1p",
                FunctionTransformer(
                    _log1p_dataframe, validate=False, feature_names_out="one-to-one"
                ),
            ),
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

    # Global per-column category vocabulary, computed from the FULL dataset
    # (df here is every fold's data combined -- build_preprocessor runs once
    # before outer-fold splitting) rather than left for each fold's
    # OneHotEncoder to infer from just its own training subset. Must include
    # the missing-token fill value whenever a column has any NaN in the full
    # dataset, since the imputer (upstream in this same pipe) can introduce
    # that value even where the raw column never had it.
    cat_categories: list[np.ndarray] = []
    for col in categorical_cols:
        # Drop real NaNs BEFORE stringifying -- astype(str) turns NaN into
        # the literal string "nan", which would silently defeat this check.
        has_na = bool(df[col].isna().any())
        values = np.asarray(df[col].dropna().unique(), dtype=str)
        if categorical_missing_strategy == "missing_token" and has_na:
            values = np.append(values, categorical_missing_token)
        cat_categories.append(np.sort(np.unique(values)))

    categorical_pipe = Pipeline(
        steps=[
            ("imputer", categorical_imputer),
            (
                "ohe",
                _make_one_hot_encoder(cat_categories if categorical_cols else "auto"),
            ),
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
    # Real, stable feature names all the way through .transform() output (not
    # plain numpy) -- fixes the root cause of LightGBM's sklearn wrapper
    # auto-generating placeholder "Column_0" names on fit and then warning
    # about a feature-name mismatch on every later predict/eval_set call
    # (confirmed 2026-07-14: the warning only fires when the underlying array
    # lacks real column names; it disappears once fit and every later call
    # consistently receive the same real-named DataFrame). Every transformer
    # in this pipeline already implements get_feature_names_out (PercentileClipper
    # and the log1p FunctionTransformer explicitly; SimpleImputer/StandardScaler/
    # OneHotEncoder natively), so this is safe.
    preprocessor.set_output(transform="pandas")
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
