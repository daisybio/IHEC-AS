from __future__ import annotations

"""Small shared utilities (logging and JSON serialization helpers)."""

from datetime import datetime
from typing import Any

import numpy as np
import pandas as pd


_LOG_LEVEL_RANK = {"none": 0, "info": 1, "debug": 2}
_CURRENT_LOG_LEVEL = "none"


def set_log_level(level: str) -> None:
    """Set global log level for vlog filtering."""
    normalized = str(level).strip().lower()
    if normalized not in _LOG_LEVEL_RANK:
        raise ValueError(f"Unsupported log level: {level}")
    global _CURRENT_LOG_LEVEL
    _CURRENT_LOG_LEVEL = normalized


def vlog(verbose: bool, message: str, level: str = "debug") -> None:
    """Print log line when verbose mode and level filter allow it."""
    if not verbose:
        return
    normalized = str(level).strip().lower()
    msg_level = _LOG_LEVEL_RANK.get(normalized, _LOG_LEVEL_RANK["debug"])
    current_level = _LOG_LEVEL_RANK.get(_CURRENT_LOG_LEVEL, _LOG_LEVEL_RANK["none"])
    if current_level >= msg_level:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] [{normalized}] {message}")


def progress_iter(
    iterable,
    *,
    total: int | None = None,
    desc: str | None = None,
    enabled: bool = True,
):
    """Return an iterable wrapped by tqdm progress bar when available.

    Falls back to the original iterable if tqdm is unavailable or disabled.
    """
    if not enabled:
        return iterable
    try:
        from tqdm.auto import tqdm

        return tqdm(iterable, total=total, desc=desc, leave=True)
    except Exception:
        return iterable


def safe_json(value: Any) -> Any:
    """Recursively convert objects to JSON-safe primitives."""
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (pd.Series, pd.Index)):
        return value.tolist()
    if isinstance(value, dict):
        return {k: safe_json(v) for k, v in value.items()}
    if isinstance(value, list):
        return [safe_json(v) for v in value]
    if isinstance(value, tuple):
        return [safe_json(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def sanitize_best_params(params: dict[str, Any]) -> dict[str, Any]:
    """Ensure GridSearch best-params payload can be serialized."""
    return {k: safe_json(v) for k, v in params.items()}
