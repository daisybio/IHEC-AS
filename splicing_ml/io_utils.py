from __future__ import annotations

"""Persistence helpers for compact artifact writing."""

import gzip
import json
import pickle
from pathlib import Path
from typing import Any

from .utils import safe_json


__all__ = ["save_pickle_gz", "save_json_gz", "load_pickle_gz", "load_json_gz"]


def save_pickle_gz(obj: Any, path: Path) -> None:
    """Save object as gzip-compressed pickle using highest protocol."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb", compresslevel=9) as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)


def save_json_gz(obj: Any, path: Path) -> None:
    """Save object as gzip-compressed UTF-8 JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=9) as f:
        json.dump(safe_json(obj), f, indent=2)


def load_pickle_gz(path: Path) -> Any:
    """Load a gzip-compressed pickle artifact."""
    with gzip.open(path, "rb") as f:
        return pickle.load(f)


def load_json_gz(path: Path) -> Any:
    """Load a gzip-compressed JSON artifact."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as f:
            return json.load(f)
    return json.loads(path.read_text(encoding="utf-8"))
