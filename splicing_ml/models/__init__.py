"""splicing_ml.models — model building, tuning, and evaluation subpackage.

Re-exports the full public API so callers can do:
    from splicing_ml.models import fit_best_estimator, evaluate_outer_fold
"""

from .beta import BetaRegressor, _inverse_logit, _logit_transform
from .deep import MLPClassifier, MLPRegressor
from .evaluator import evaluate_outer_fold, tune_threshold_balanced_accuracy
from .grids import (
    build_param_candidates,
    choose_param_grid,
    choose_param_lhs_candidates,
)
from .search import fit_best_estimator, metric_key, scorer_name
from .cuml_utils import cuml_gpu_available
from .xgb_utils import (
    _set_xgb_cpu_predictor_for_inference,
    _unwrap_model_step,
    xgb_gpu_available,
)

__all__ = [
    "BetaRegressor",
    "_inverse_logit",
    "_logit_transform",
    "MLPClassifier",
    "MLPRegressor",
    "evaluate_outer_fold",
    "tune_threshold_balanced_accuracy",
    "build_param_candidates",
    "choose_param_grid",
    "choose_param_lhs_candidates",
    "fit_best_estimator",
    "metric_key",
    "scorer_name",
    "cuml_gpu_available",
    "_set_xgb_cpu_predictor_for_inference",
    "_unwrap_model_step",
    "xgb_gpu_available",
]
