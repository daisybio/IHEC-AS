from __future__ import annotations

"""Data loading, subset generation, and target construction."""

import time
from typing import Any

import pandas as pd
import numpy as np

from .config import REQUIRED_COLUMNS, SubsetConfig
from .utils import vlog


__all__ = [
    "data_reader_capabilities",
    "data_reader_startup_warnings",
    "load_dataset",
    "benchmark_data_readers",
    "generate_subset_configs",
    "subset_dataframe",
    "build_targets",
]


def _logit_transform(psi: np.ndarray, epsilon: float = 1e-7) -> np.ndarray:
    """Transform PSI from [0,1] to unbounded logit scale.

    Parameters
    ----------
    psi : np.ndarray
        Array of values in [0, 1].
    epsilon : float
        Clipping value to avoid log(0) or log(∞).

    Returns
    -------
    np.ndarray
        Logit-transformed values, unbounded.
    """
    psi = np.asarray(psi, dtype=float)
    psi_clipped = np.clip(psi, epsilon, 1.0 - epsilon)
    return np.log(psi_clipped / (1.0 - psi_clipped))


def _is_module_available(module_name: str) -> bool:
    """Best-effort optional dependency presence check."""
    try:
        __import__(module_name)
        return True
    except Exception:
        return False


def data_reader_capabilities() -> dict[str, bool]:
    """Report availability of optional fast-reader dependencies."""
    return {
        "polars": _is_module_available("polars"),
        "pyarrow": _is_module_available("pyarrow"),
    }


def data_reader_startup_warnings(reader_backend: str) -> list[str]:
    """Return human-readable startup warnings for configured backend."""
    caps = data_reader_capabilities()
    warnings: list[str] = []

    if reader_backend in {"auto", "polars"} and not caps["polars"]:
        if reader_backend == "polars":
            warnings.append(
                "data reader backend is set to polars but polars is not installed; loading will fail unless backend is changed"
            )
        else:
            warnings.append(
                "polars is not installed; auto backend will skip the fastest path and fall back to pandas"
            )

    if reader_backend in {"auto", "pandas"} and not caps["pyarrow"]:
        warnings.append(
            "pyarrow is not installed; pandas will fall back to the C parser which may be slower on large CSV files"
        )

    return warnings


def _load_dataset_pandas(path: str, verbose: bool = False) -> pd.DataFrame:
    """Read CSV using pandas with a pyarrow-first strategy."""
    t0 = time.perf_counter()
    try:
        # PyArrow parser is often faster than the default C parser on large CSVs.
        df = pd.read_csv(path, compression="infer", engine="pyarrow")
        dt = time.perf_counter() - t0
        vlog(verbose, f"Loaded with pandas(pyarrow) in {dt:.2f}s")
        return df
    except Exception as exc:
        vlog(verbose, f"pandas(pyarrow) read failed, falling back to pandas(c): {exc}")
        df = pd.read_csv(path, compression="infer", low_memory=False)
        dt = time.perf_counter() - t0
        vlog(verbose, f"Loaded with pandas(c) in {dt:.2f}s")
        return df


def _load_dataset_polars(path: str, verbose: bool = False) -> pd.DataFrame:
    """Read CSV with polars and convert to pandas for sklearn compatibility."""
    t0 = time.perf_counter()
    try:
        import polars as pl
    except Exception as exc:
        raise RuntimeError(
            "polars backend requested but polars is not installed"
        ) from exc

    df_pl = pl.read_csv(path, infer_schema_length=10000, try_parse_dates=False)
    df = df_pl.to_pandas()
    dt = time.perf_counter() - t0
    vlog(verbose, f"Loaded with polars->pandas in {dt:.2f}s")
    return df


def load_dataset(
    path: str,
    verbose: bool = False,
    reader_backend: str = "auto",
) -> pd.DataFrame:
    """Load table from CSV/CSV.GZ and validate required columns."""
    vlog(
        verbose,
        f"Loading dataset from {path} (backend={reader_backend})",
        level="info",
    )

    if reader_backend not in {"auto", "pandas", "polars"}:
        raise ValueError("reader_backend must be one of: auto, pandas, polars")

    if reader_backend == "pandas":
        df = _load_dataset_pandas(path, verbose=verbose)
    elif reader_backend == "polars":
        df = _load_dataset_polars(path, verbose=verbose)
    else:
        try:
            df = _load_dataset_polars(path, verbose=verbose)
        except Exception as exc:
            vlog(verbose, f"Auto backend fallback to pandas due to: {exc}")
            df = _load_dataset_pandas(path, verbose=verbose)

    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    vlog(verbose, f"Dataset loaded with shape={df.shape}", level="info")
    return df


def benchmark_data_readers(
    path: str,
    runs: int = 1,
    verbose: bool = False,
) -> dict[str, Any]:
    """Benchmark auto/pandas/polars reader backends on one dataset."""
    runs = max(1, int(runs))
    backends = ["auto", "pandas", "polars"]
    results: list[dict[str, Any]] = []

    for backend in backends:
        times: list[float] = []
        shape: tuple[int, int] | None = None
        error_msg = ""
        for _ in range(runs):
            t0 = time.perf_counter()
            try:
                df = load_dataset(path, verbose=verbose, reader_backend=backend)
                elapsed = time.perf_counter() - t0
                times.append(elapsed)
                shape = (int(df.shape[0]), int(df.shape[1]))
            except Exception as exc:
                error_msg = str(exc)
                break

        if times:
            mean_sec = float(np.mean(times))
            min_sec = float(np.min(times))
            max_sec = float(np.max(times))
        else:
            mean_sec = float("inf")
            min_sec = float("inf")
            max_sec = float("inf")

        results.append(
            {
                "backend": backend,
                "ok": bool(times),
                "runs_completed": len(times),
                "runs_requested": runs,
                "mean_seconds": mean_sec,
                "min_seconds": min_sec,
                "max_seconds": max_sec,
                "shape": list(shape) if shape is not None else [],
                "error": error_msg,
            }
        )

    ok_results = [r for r in results if r["ok"]]
    winner = (
        min(ok_results, key=lambda r: r["mean_seconds"])["backend"]
        if ok_results
        else ""
    )

    return {
        "path": path,
        "runs": runs,
        "capabilities": data_reader_capabilities(),
        "results": results,
        "winner_backend": winner,
    }


def generate_subset_configs(
    df: pd.DataFrame, verbose: bool = False
) -> list[SubsetConfig]:
    """Create all requested subset combinations available in the data."""
    event_types = ["SE", "RI"]
    transcript_filters = ["transcripts", "tsl_filtered", "biotype_filtered"]
    variability_values = ["High", "Low", "both"]
    group_cols = ["seqnames", "ontology"]

    available_event_types = set(df["Event Type"].dropna().unique().tolist())
    available_transcript_filters = set(
        df["transcript_filter"].dropna().unique().tolist()
    )
    available_variability = set(df["Variability"].dropna().unique().tolist())

    configs: list[SubsetConfig] = []
    for event_type in event_types:
        if event_type not in available_event_types:
            continue
        for transcript_filter in transcript_filters:
            if transcript_filter not in available_transcript_filters:
                continue
            for variability in variability_values:
                if variability != "both" and variability not in available_variability:
                    continue
                for group_col in group_cols:
                    configs.append(
                        SubsetConfig(
                            event_type=event_type,
                            transcript_filter=transcript_filter,
                            variability=variability,
                            group_col=group_col,
                        )
                    )
    vlog(verbose, f"Generated {len(configs)} subset configurations", level="info")
    return configs


def subset_dataframe(
    df: pd.DataFrame, cfg: SubsetConfig, verbose: bool = False
) -> pd.DataFrame:
    """Filter a full table to one subset configuration."""
    sel = (df["Event Type"] == cfg.event_type) & (
        df["transcript_filter"] == cfg.transcript_filter
    )
    if cfg.variability != "both":
        sel = sel & (df["Variability"] == cfg.variability)
    out = df.loc[sel].copy().dropna(subset=["PSI", "seqnames", "ontology"])
    vlog(
        verbose,
        "Subset "
        f"event_type={cfg.event_type}, transcript_filter={cfg.transcript_filter}, "
        f"variability={cfg.variability}, group_col={cfg.group_col} -> rows={out.shape[0]}",
    )
    return out


def build_targets(
    df: pd.DataFrame,
    task: str,
    low_thr: float,
    high_thr: float,
    verbose: bool = False,
    use_logit: bool = True,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Build regression targets or binarized classification labels from PSI.

    For regression, optionally applies logit transform to model PSI on unbounded scale.
    """
    if task == "regression":
        y = df["PSI"].astype(float).to_numpy()
        if use_logit:
            y = _logit_transform(y)
            vlog(verbose, f"Built regression target (logit-scale) with n={y.size}")
        else:
            vlog(verbose, f"Built regression target with n={y.size}")
        return df, y

    y_raw = df["PSI"].astype(float)
    keep = (y_raw <= low_thr) | (y_raw >= high_thr)
    df2 = df.loc[keep].copy()
    y = (df2["PSI"].astype(float) >= high_thr).astype(int).to_numpy()
    vlog(
        verbose,
        f"Built classification target with thresholds=({low_thr:.4f}, {high_thr:.4f}), "
        f"kept_n={y.size}",
    )
    return df2, y
