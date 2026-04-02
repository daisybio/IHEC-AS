from __future__ import annotations

"""PyTorch MLP estimators with a scikit-learn compatible interface.

Two public estimators are provided:

- ``MLPRegressor``    — feeds into the regression nested-CV path.
- ``MLPClassifier``   — feeds into the classification nested-CV path.

Both are drop-in replacements for any other sklearn estimator in the
existing ``sklearn.Pipeline(preprocessor → model)`` pattern.  They are
registered as model type ``"mlp"`` in ``config.ALL_MODEL_TYPES`` and their
hyperparameter candidates are built by ``grids.build_param_candidates``.

Design principles
-----------------
* **Device-agnostic** — uses ``torch.accelerator.is_available()`` (PyTorch ≥ 2.4)
  with a graceful fallback for older installs (CUDA → MPS → CPU).
* **Regularisation** — L2 weight-decay via Adam, Dropout, BatchNorm1d.
* **Early stopping** — patience-based, restores best-validation-loss weights.
* **Internal train/val split** — ``val_fraction`` of the training data is held
  out for early-stopping validation; this is separate from the outer/inner CV
  folds managed by the pipeline.
* **Strict sklearn parameter convention** — every ``__init__`` argument is
  stored as an identical attribute so that ``get_params()`` / ``set_params()``
  / ``clone()`` and ``GridSearchCV`` work correctly.
"""

import copy
from typing import Any, Sequence

import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin, RegressorMixin
from sklearn.utils.validation import check_is_fitted

from ..config import RNG_SEED
from ._contexts import _mlp_val_indices_context, _es_pp_val_context

__all__ = ["MLPRegressor", "MLPClassifier"]


# ---------------------------------------------------------------------------
# Early-stopping helper
# ---------------------------------------------------------------------------


class _EarlyStopping:
    """Patience-based early stopping with best-weight restoration.

    Parameters
    ----------
    patience : int
        Number of epochs without improvement before stopping.
    """

    def __init__(self, patience: int = 5) -> None:
        """Initialize a _EarlyStopping instance."""
        self.patience = patience
        self.best_loss: float = float("inf")
        self.counter: int = 0
        self.best_weights: dict | None = None

    def step(self, val_loss: float, model: Any) -> bool:
        """Update state and return True when training should stop."""
        if val_loss < self.best_loss:
            self.best_loss = val_loss
            self.counter = 0
            # Deep copy so subsequent weight updates don't corrupt the snapshot.
            self.best_weights = copy.deepcopy(model.state_dict())
        else:
            self.counter += 1
        return self.counter >= self.patience

    def restore_best(self, model: Any) -> None:
        """Load the best-seen weights back into *model* (in-place)."""
        if self.best_weights is not None:
            model.load_state_dict(self.best_weights)


# ---------------------------------------------------------------------------
# PyTorch network
# ---------------------------------------------------------------------------


def _build_mlp(
    n_features: int, hidden_sizes: Sequence[int], dropout_rate: float
) -> Any:
    """Build and return a ``torch.nn.Module`` MLP.

    Architecture: Input → (Linear → BatchNorm1d → ReLU → Dropout) × N → Linear.

    The output layer has no activation; the loss functions handle the final
    transformation (sigmoid inside BCEWithLogitsLoss; identity for MSELoss).

    Parameters
    ----------
    n_features : int
        Number of input features.
    hidden_sizes : sequence of int
        Width of each hidden layer.
    dropout_rate : float
        Dropout probability applied after each hidden ReLU.
    """
    import torch.nn as nn

    layers: list[nn.Module] = []
    in_dim = int(n_features)
    for h in hidden_sizes:
        layers.append(nn.Linear(in_dim, int(h)))
        layers.append(nn.BatchNorm1d(int(h)))
        layers.append(nn.ReLU())
        layers.append(nn.Dropout(p=float(dropout_rate)))
        in_dim = int(h)
    layers.append(nn.Linear(in_dim, 1))
    return nn.Sequential(*layers)


# ---------------------------------------------------------------------------
# sklearn-compatible base
# ---------------------------------------------------------------------------


class _BaseMLP(BaseEstimator):
    """Base sklearn estimator wrapping a PyTorch MLP.

    Subclasses must implement ``_criterion()`` (returns the loss function) and
    ``_n_outputs()`` (always 1 for current tasks).
    """

    def __init__(
        self,
        hidden_sizes: tuple[int, ...] = (128, 64),
        dropout_rate: float = 0.2,
        learning_rate: float = 1e-3,
        weight_decay: float = 1e-4,
        batch_size: int = 128,
        max_epochs: int = 150,
        early_stopping_patience: int = 5,
        val_fraction: float = 0.15,
        random_state: int = RNG_SEED,
    ) -> None:
        # All params stored exactly as received — required by sklearn convention.
        """Initialize a _BaseMLP instance."""
        self.hidden_sizes = hidden_sizes
        self.dropout_rate = dropout_rate
        self.learning_rate = learning_rate
        self.weight_decay = weight_decay
        self.batch_size = batch_size
        self.max_epochs = max_epochs
        self.early_stopping_patience = early_stopping_patience
        self.val_fraction = val_fraction
        self.random_state = random_state

    # ------------------------------------------------------------------
    # Device detection
    # ------------------------------------------------------------------

    def _get_device(self) -> Any:
        """Return the best available ``torch.device``.

        Tries ``torch.accelerator.is_available()`` first (PyTorch ≥ 2.4).
        Falls back to the legacy CUDA → MPS → CPU chain.
        """
        import torch

        try:
            if torch.accelerator.is_available():
                # current_accelerator() returns e.g. device("cuda") or device("mps")
                return torch.device(torch.accelerator.current_accelerator())
        except AttributeError:
            # torch.accelerator not present in PyTorch < 2.4
            pass

        if torch.cuda.is_available():
            return torch.device("cuda")
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")

    # ------------------------------------------------------------------
    # Input conversion
    # ------------------------------------------------------------------

    @staticmethod
    def _to_numpy(X: Any) -> np.ndarray:
        """Convert arbitrary array-like (including scipy sparse) to float32 ndarray."""
        if hasattr(X, "toarray"):  # scipy sparse matrix
            X = X.toarray()
        arr = np.asarray(X, dtype=np.float32)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        return arr

    # ------------------------------------------------------------------
    # Subclass hooks
    # ------------------------------------------------------------------

    def _criterion(self) -> Any:
        """Internal helper for criterion."""
        raise NotImplementedError

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------

    def fit(self, X: Any, y: Any, val_indices: np.ndarray | None = None) -> "_BaseMLP":
        """Fit the MLP on *X* / *y*.

        Internally carves out ``val_fraction`` of the data as a validation
        set for early stopping. If ``val_indices`` is provided (e.g., from
        cross-validation), uses those specific indices; otherwise uses a
        sequential split at the end of the data to respect any implicit
        grouping (e.g., seqnames or ontology clustering).

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Feature matrix.  Scipy sparse matrices are accepted.
        y : array-like of shape (n_samples,)
            Target values (float for regression; 0/1 float for classification).
        val_indices : array-like of int, optional
            Indices of samples to use for validation (early stopping).
            If provided, overrides ``val_fraction``. Typically passed from
            GridSearchCV/OptunaSearchCV inner CV to avoid data leakage.
        """
        import torch

        rng = np.random.default_rng(int(self.random_state))
        torch.manual_seed(int(self.random_state))

        X_arr = self._to_numpy(X)
        y_arr = np.asarray(y, dtype=np.float32).ravel()
        n_samples, n_features = X_arr.shape
        self.n_features_in_ = n_features

        # --- train / val split ---
        if val_indices is not None:
            # Explicit indices take priority (e.g. final refit after search).
            val_idx = np.asarray(val_indices, dtype=int)
            train_idx = np.setdiff1d(np.arange(n_samples), val_idx)
        elif (
            hasattr(_es_pp_val_context, "X_val_pp")
            and _es_pp_val_context.X_val_pp is not None
        ):
            # Preprocessed val data injected by _ESPipeline — use all X as train
            X_tr, y_tr = X_arr, y_arr
            X_va = self._to_numpy(_es_pp_val_context.X_val_pp)
            y_va = np.asarray(_es_pp_val_context.y_val, dtype=np.float32).ravel()
        elif (
            hasattr(_mlp_val_indices_context, "indices")
            and _mlp_val_indices_context.indices is not None
        ):
            # Legacy fallback (kept for backward compat; should not fire in pipeline)
            val_idx = np.asarray(_mlp_val_indices_context.indices, dtype=int)
            train_idx = np.setdiff1d(np.arange(n_samples), val_idx)
        else:
            # Standalone mode: sequential split respecting data order.
            val_size = max(2, int(n_samples * float(self.val_fraction)))
            train_idx = np.arange(n_samples - val_size)
            val_idx = np.arange(n_samples - val_size, n_samples)

        if val_indices is None and not (
            hasattr(_es_pp_val_context, "X_val_pp")
            and _es_pp_val_context.X_val_pp is not None
        ):
            X_tr, y_tr = X_arr[train_idx], y_arr[train_idx]
            X_va, y_va = X_arr[val_idx], y_arr[val_idx]

        # --- model, optimiser, scheduler, criterion ----------------------
        # hidden_sizes may arrive as a comma-separated string when sampled by
        # Optuna (CategoricalDistribution requires scalar choices, not tuples).
        _hidden = self.hidden_sizes
        if isinstance(_hidden, str):
            _hidden = tuple(int(x) for x in _hidden.split(",") if x.strip())
        device = self._get_device()
        net = _build_mlp(
            n_features=n_features,
            hidden_sizes=_hidden,
            dropout_rate=self.dropout_rate,
        ).to(device)

        optimizer = torch.optim.Adam(
            net.parameters(),
            lr=float(self.learning_rate),
            weight_decay=float(self.weight_decay),
        )
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode="min", factor=0.5, patience=3, min_lr=1e-6
        )
        criterion = self._criterion()
        stopper = _EarlyStopping(patience=int(self.early_stopping_patience))

        # Pre-allocate tensors on device once (avoids per-epoch transfers).
        X_tr_t = torch.from_numpy(X_tr).to(device)
        y_tr_t = torch.from_numpy(y_tr).to(device)
        X_va_t = torch.from_numpy(X_va).to(device)
        y_va_t = torch.from_numpy(y_va).to(device)

        n_tr = X_tr_t.shape[0]
        bs = max(2, int(self.batch_size))  # ensure batch >= 2 for BatchNorm

        # Mixed precision: ~1.5-2× speedup on CUDA with negligible accuracy cost.
        use_amp = device.type == "cuda"
        scaler = torch.amp.GradScaler(enabled=use_amp)

        # Per-epoch generator so shuffle is reproducible per epoch.
        gen = torch.Generator(device="cpu")

        # --- training loop -----------------------------------------------
        for epoch in range(int(self.max_epochs)):
            net.train()
            gen.manual_seed(int(self.random_state) + epoch)
            perm = torch.randperm(n_tr, generator=gen)

            for start in range(0, n_tr, bs):
                batch_idx = perm[start : start + bs]
                if batch_idx.shape[0] < 2:
                    # BatchNorm1d in train mode requires >= 2 samples per batch.
                    continue
                xb = X_tr_t[batch_idx]
                yb = y_tr_t[batch_idx]
                optimizer.zero_grad()
                with torch.amp.autocast(device_type=device.type, enabled=use_amp):
                    out = net(xb).squeeze(-1)
                    loss = criterion(out, yb)
                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()

            # Validation pass — eval mode uses running BN statistics.
            net.eval()
            with torch.no_grad():
                with torch.amp.autocast(device_type=device.type, enabled=use_amp):
                    val_out = net(X_va_t).squeeze(-1)
                val_loss = criterion(val_out.float(), y_va_t).item()

            scheduler.step(val_loss)

            if stopper.step(val_loss, net):
                break

        # Always restore the best-validation weights.
        stopper.restore_best(net)
        self.model_ = net
        self.device_ = device
        self.n_epochs_trained_ = epoch + 1  # Store for diagnostics
        return self

    # ------------------------------------------------------------------
    # Inference helper
    # ------------------------------------------------------------------

    def _forward(self, X: Any) -> np.ndarray:
        """Run inference and return raw logits / regression values as numpy."""
        import torch

        check_is_fitted(self, "model_")
        self.model_.eval()
        with torch.no_grad():
            X_arr = self._to_numpy(X)
            X_t = torch.from_numpy(X_arr).to(self.device_)
            out = self.model_(X_t).squeeze(-1).cpu().numpy()
        return out


# ---------------------------------------------------------------------------
# Public estimators
# ---------------------------------------------------------------------------


class MLPRegressor(RegressorMixin, _BaseMLP):
    """MLP regressor with PyTorch backend and sklearn interface.

    Predicts on the same (possibly logit-transformed) scale as the targets
    passed to ``fit``; the pipeline's ``evaluate_outer_fold`` handles any
    inverse transformation.

    Attributes (set during fit)
    ---------------------------
    n_epochs_trained_ : int
        Number of epochs completed before early stopping (or max_epochs if
        no early stopping triggered). Use this to assess whether the epoch
        budget is sufficient: if consistently equals max_epochs, consider
        increasing it.

    Parameters
    ----------
    hidden_sizes : tuple of int, default (128, 64)
        Width of each hidden layer.
    dropout_rate : float, default 0.2
        Dropout probability after each hidden activation.
    learning_rate : float, default 1e-3
        Initial Adam learning rate.
    weight_decay : float, default 1e-4
        L2 penalty (Adam ``weight_decay``).
    batch_size : int, default 128
        Mini-batch size.  Clamped to >= 2 during training.
    max_epochs : int, default 75
        Maximum training epochs.
    early_stopping_patience : int, default 5
        Stop after this many epochs without validation-loss improvement.
    val_fraction : float, default 0.15
        Fraction of training data reserved for early-stopping validation.
    random_state : int, default 42
        Seed for weight initialisation and data shuffling.
    """

    def _criterion(self) -> Any:
        """Internal helper for criterion."""
        import torch.nn as nn

        return nn.MSELoss()

    def predict(self, X: Any) -> np.ndarray:
        """Return continuous predictions for *X*."""
        return self._forward(X).ravel()


class MLPClassifier(ClassifierMixin, _BaseMLP):
    """MLP binary classifier with PyTorch backend and sklearn interface.

    Outputs probabilities compatible with ``CalibratedClassifierCV`` and
    the pipeline's threshold-tuning logic.

    Attributes (set during fit)
    ---------------------------
    n_epochs_trained_ : int
        Number of epochs completed before early stopping (or max_epochs if
        no early stopping triggered). Use this to assess whether the epoch
        budget is sufficient: if consistently equals max_epochs, consider
        increasing it.

    Parameters
    ----------
    hidden_sizes : tuple of int, default (128, 64)
        Width of each hidden layer.
    dropout_rate : float, default 0.2
        Dropout probability after each hidden activation.
    learning_rate : float, default 1e-3
        Initial Adam learning rate.
    weight_decay : float, default 1e-4
        L2 penalty (Adam ``weight_decay``).
    batch_size : int, default 128
        Mini-batch size.  Clamped to >= 2 during training.
    max_epochs : int, default 75
        Maximum training epochs.
    early_stopping_patience : int, default 5
        Stop after this many epochs without validation-loss improvement.
    val_fraction : float, default 0.15
        Fraction of training data reserved for early-stopping validation.
    random_state : int, default 42
        Seed for weight initialisation and data shuffling.
    """

    def fit(self, X: Any, y: Any) -> "MLPClassifier":
        """Fit the classifier.  Sets ``classes_`` before delegating to base."""
        y_arr = np.asarray(y)
        self.classes_ = np.unique(y_arr)
        return super().fit(X, y_arr)

    def _criterion(self) -> Any:
        """Internal helper for criterion."""
        import torch.nn as nn

        # Targets must be float32 in [0, 1]; the base fit() already converts y.
        return nn.BCEWithLogitsLoss()

    def predict_proba(self, X: Any) -> np.ndarray:
        """Return class probabilities of shape (n_samples, 2).

        Column 0 is P(class=0), column 1 is P(class=1).
        """
        import torch

        logits = self._forward(X)  # shape (n,)
        probs_pos = torch.sigmoid(torch.from_numpy(logits)).numpy()
        probs_neg = 1.0 - probs_pos
        return np.column_stack([probs_neg, probs_pos])

    def predict(self, X: Any) -> np.ndarray:
        """Return binary class labels using a 0.5 threshold."""
        return (self.predict_proba(X)[:, 1] >= 0.5).astype(int)
