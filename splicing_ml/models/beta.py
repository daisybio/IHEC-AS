from __future__ import annotations

"""BetaRegressor: statsmodels beta regression wrapped as a sklearn estimator.

This module also provides the logit/inverse-logit scale helpers used in
training and evaluation of PSI regression targets.
"""

import warnings
from typing import Any

import numpy as np
from sklearn.base import BaseEstimator, RegressorMixin

from ..config import RNG_SEED

__all__ = ["BetaRegressor", "_inverse_logit", "_logit_transform"]


def _inverse_logit(logit_p: np.ndarray) -> np.ndarray:
    """Transform logit scale back to [0, 1].

    Parameters
    ----------
    logit_p : np.ndarray
        Unbounded logit values.

    Returns
    -------
    np.ndarray
        Values clipped to [0, 1].
    """
    logit_p = np.asarray(logit_p, dtype=float)
    result = 1.0 / (1.0 + np.exp(-logit_p))
    return np.clip(result, 0.0, 1.0)


def _logit_transform(psi: np.ndarray, epsilon: float = 1e-7) -> np.ndarray:
    """Transform PSI from [0, 1] to unbounded logit scale."""
    psi = np.asarray(psi, dtype=float)
    psi_clipped = np.clip(psi, epsilon, 1.0 - epsilon)
    return np.log(psi_clipped / (1.0 - psi_clipped))


class BetaRegressor(BaseEstimator, RegressorMixin):
    """Statsmodels-backed beta regression wrapped as an sklearn regressor.

    This model expects a target on the PSI scale (open interval (0, 1)).
    Predictions are returned on the same PSI scale.

    Parameters
    ----------
    epsilon : float
        Clipping value to keep targets strictly within (0, 1).
    maxiter : int
        Maximum number of optimization iterations per solver attempt.
    max_sparse_columns : int
        Upper bound on sparse-matrix columns before column truncation.
    max_dense_elements : int
        Upper bound on total dense-matrix elements (rows × cols).
    max_train_rows : int
        Random subset size used when training data exceeds this limit.
    """

    def __init__(
        self,
        epsilon: float = 1e-7,
        maxiter: int = 200,
        max_sparse_columns: int = 2000,
        max_dense_elements: int = 30_000_000,
        max_train_rows: int = 5000,
    ):
        self.epsilon = float(epsilon)
        self.maxiter = int(maxiter)
        self.max_sparse_columns = int(max_sparse_columns)
        self.max_dense_elements = int(max_dense_elements)
        self.max_train_rows = int(max_train_rows)

    # ------------------------------------------------------------------
    # Private helpers — each handles one step of fit()
    # ------------------------------------------------------------------

    @staticmethod
    def _as_dense_2d(x: Any) -> np.ndarray:
        """Convert preprocessor output to a dense 2D float array."""
        if hasattr(x, "toarray"):
            x_arr = x.toarray()
        else:
            x_arr = np.asarray(x)
        x_arr = np.asarray(x_arr, dtype=float)
        if x_arr.ndim == 1:
            x_arr = x_arr.reshape(-1, 1)
        return x_arr

    def _filter_sparse_columns(self, x: Any) -> tuple[Any, np.ndarray | None]:
        """Truncate to the most informative columns when input is sparse.

        Returns the (possibly truncated) input and the kept column indices
        (None if no truncation was applied).
        """
        if not (hasattr(x, "getnnz") and hasattr(x, "shape")):
            return x, None
        n_sparse_cols = int(x.shape[1])
        max_sparse_cols = max(1, int(self.max_sparse_columns))
        if n_sparse_cols <= max_sparse_cols:
            return x, None
        # Keep the columns with most non-zeros (highest information density).
        nnz = np.asarray(x.getnnz(axis=0)).ravel()
        keep_idx = np.argsort(nnz)[-max_sparse_cols:]
        keep_idx = np.sort(np.asarray(keep_idx, dtype=int))
        return x[:, keep_idx], keep_idx

    def _remove_degenerate_columns(
        self, x_arr: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
        """Remove constant, non-finite, and linearly dependent columns.

        Returns (filtered_array, feature_mask, qr_keep_idx).
        feature_mask is a boolean mask over the original column set.
        qr_keep_idx is an index array for the QR-selected columns (None if
        QR was not applied or failed).
        """
        # Step 1: remove constant and non-finite columns.
        finite_col_mask = np.all(np.isfinite(x_arr), axis=0)
        var_col_mask = np.std(x_arr, axis=0) > 1e-12
        feature_mask = finite_col_mask & var_col_mask
        if not np.any(feature_mask):
            # Fallback: keep at least one column to avoid empty design matrix.
            feature_mask = np.zeros(x_arr.shape[1], dtype=bool)
            feature_mask[0] = True
        x_use = x_arr[:, feature_mask]

        # Step 2: QR factorisation to remove linear dependencies (common after
        # one-hot expansions) and stabilise the beta-regression Hessian solves.
        qr_keep_idx = None
        try:
            from scipy.linalg import qr as scipy_qr

            _, r_mat, piv = scipy_qr(x_use, mode="economic", pivoting=True)
            if r_mat.size:
                diag = np.abs(np.diag(r_mat))
                tol = max(
                    np.max(diag) * 1e-8,
                    np.max(diag) * max(x_use.shape) * np.finfo(float).eps,
                )
                rank = int(np.sum(diag > tol))
                qr_keep_idx = np.sort(np.asarray(piv[: max(1, rank)], dtype=int))
                x_use = x_use[:, qr_keep_idx]
        except Exception:
            pass

        return x_use, feature_mask, qr_keep_idx

    def _optimize_model(self, x_design: np.ndarray, y_clipped: np.ndarray) -> Any:
        """Fit BetaModel trying multiple optimisers and covariance types.

        Raises RuntimeError if all attempts fail.
        """
        try:
            from statsmodels.othermod.betareg import BetaModel
        except Exception as exc:
            raise RuntimeError(
                "beta model requested but statsmodels beta regression support is unavailable"
            ) from exc

        fit_kw_variants = [
            {"cov_type": "none", "skip_hessian": True},
            {"cov_type": "none"},
            {},
        ]
        fit_errors: list[str] = []

        # Strict beta-only policy: retry with multiple optimizers but never
        # fall back to a non-beta estimator.
        for method in ("bfgs", "lbfgs", "newton"):
            for fit_kw in fit_kw_variants:
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        result = BetaModel(endog=y_clipped, exog=x_design).fit(
                            method=method,
                            maxiter=self.maxiter,
                            disp=False,
                            **fit_kw,
                        )
                    self.fit_backend_ = (
                        f"beta:{method}:{'+'.join(sorted(fit_kw.keys()))}"
                        if fit_kw
                        else f"beta:{method}:default"
                    )
                    return result
                except TypeError:
                    # Older statsmodels may not accept newer fit kwargs.
                    continue
                except Exception as exc:
                    fit_errors.append(f"{method}/{fit_kw}: {exc}")
            else:
                continue
            break  # noqa: SIM105 — explicit loop break after inner match

        raise RuntimeError(
            "beta regression failed to converge "
            f"(n={y_clipped.shape[0]}, p={x_design.shape[1]}). "
            f"attempts={'; '.join(fit_errors)}"
        )

    # ------------------------------------------------------------------
    # sklearn interface
    # ------------------------------------------------------------------

    def fit(self, x: Any, y: Any) -> "BetaRegressor":
        """Fit beta regression on PSI-scale targets.

        Applies sparse-column filtering, degenerate-column removal, optional
        row subsampling, and iterative optimizer fallback.
        """
        from statsmodels.tools.tools import add_constant

        # Sparse column filtering (stores index for predict).
        x_in, self._sparse_keep_idx_ = self._filter_sparse_columns(x)

        # Convert to dense 2D float array.
        x_arr = self._as_dense_2d(x_in)

        # Validate dense matrix size before allocating.
        dense_elements = int(x_arr.shape[0]) * int(x_arr.shape[1])
        if dense_elements > max(1, int(self.max_dense_elements)):
            raise RuntimeError(
                "beta regression design matrix too large for dense statsmodels fit "
                f"(n={x_arr.shape[0]}, p={x_arr.shape[1]}, elements={dense_elements}, "
                f"limit={self.max_dense_elements})"
            )

        # Clip targets to open (0, 1) interval required by beta distribution.
        y_arr = np.asarray(y, dtype=float)
        y_arr = np.clip(y_arr, self.epsilon, 1.0 - self.epsilon)

        # Optional row subsampling to control fit time for large datasets.
        if x_arr.shape[0] > max(1, int(self.max_train_rows)):
            rng = np.random.default_rng(RNG_SEED)
            fit_idx = rng.choice(
                x_arr.shape[0],
                size=max(1, int(self.max_train_rows)),
                replace=False,
            )
            fit_idx = np.sort(np.asarray(fit_idx, dtype=int))
            x_arr = x_arr[fit_idx]
            y_arr = y_arr[fit_idx]

        # Degenerate column removal and QR rank reduction.
        x_use, self._feature_mask_, self._qr_keep_idx_ = (
            self._remove_degenerate_columns(x_arr)
        )

        x_design = add_constant(x_use, has_constant="skip")
        self._beta_result_ = self._optimize_model(x_design, y_arr)
        self.n_features_in_ = x_arr.shape[1]
        return self

    def predict(self, x: Any) -> np.ndarray:
        """Predict PSI values in (0, 1) for new observations."""
        from statsmodels.tools.tools import add_constant

        if not hasattr(self, "_beta_result_"):
            raise RuntimeError("BetaRegressor is not fitted")

        x_in = x
        if getattr(self, "_sparse_keep_idx_", None) is not None and hasattr(x_in, "shape"):
            x_in = x_in[:, self._sparse_keep_idx_]

        x_arr = self._as_dense_2d(x_in)
        if hasattr(self, "_feature_mask_"):
            x_use = x_arr[:, self._feature_mask_]
        else:
            x_use = x_arr
        if getattr(self, "_qr_keep_idx_", None) is not None:
            x_use = x_use[:, self._qr_keep_idx_]

        x_design = add_constant(x_use, has_constant="skip")
        pred = np.asarray(self._beta_result_.predict(x_design), dtype=float)
        return np.clip(pred, self.epsilon, 1.0 - self.epsilon)
