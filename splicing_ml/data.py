from __future__ import annotations

"""Data loading, subset generation, and target construction."""

import importlib.util
import os
import time
from typing import Any

import pandas as pd
import numpy as np

from .config import REQUIRED_COLUMNS, SubsetConfig
from .utils import vlog


__all__ = [
    "data_reader_capabilities",
    "data_reader_startup_warnings",
    "recommend_polars_max_threads",
    "estimate_filter_fraction_for_path",
    "recommend_auto_data_reader_backend",
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
    return importlib.util.find_spec(module_name) is not None


def _available_memory_bytes() -> int | None:
    """Return best-effort available memory (bytes), accounting for cgroup limits."""
    mem_available: int | None = None

    # Linux primary signal for reclaimable memory without swapping.
    try:
        with open("/proc/meminfo", "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    parts = line.split()
                    mem_available = int(parts[1]) * 1024
                    break
    except Exception:
        mem_available = None

    # If the process is cgroup-constrained, available memory cannot exceed that.
    try:
        with open("/proc/self/cgroup", "r", encoding="utf-8") as f:
            cgroup_path = f.readline().strip().split(":", 2)[-1]
        base = os.path.join("/sys/fs/cgroup", cgroup_path.lstrip("/"))
        lim_path = os.path.join(base, "memory.max")
        cur_path = os.path.join(base, "memory.current")
        if os.path.exists(lim_path):
            with open(lim_path, "r", encoding="utf-8") as f:
                lim_raw = f.read().strip()
            if lim_raw != "max":
                cgroup_limit = int(lim_raw)
                cgroup_current = 0
                if os.path.exists(cur_path):
                    with open(cur_path, "r", encoding="utf-8") as f:
                        cgroup_current = int(f.read().strip())
                cgroup_available = max(0, cgroup_limit - cgroup_current)
                mem_available = (
                    cgroup_available
                    if mem_available is None
                    else min(mem_available, cgroup_available)
                )
    except Exception:
        pass

    # SLURM allocations may enforce memory limits even when cgroup memory.max is "max".
    try:
        slurm_limit_bytes: int | None = None
        mem_per_node = os.environ.get("SLURM_MEM_PER_NODE")
        if mem_per_node and mem_per_node.isdigit():
            slurm_limit_bytes = int(mem_per_node) * 1024 * 1024
        else:
            mem_per_cpu = os.environ.get("SLURM_MEM_PER_CPU")
            cpus_per_task = os.environ.get("SLURM_CPUS_PER_TASK", "1")
            if mem_per_cpu and mem_per_cpu.isdigit() and cpus_per_task.isdigit():
                slurm_limit_bytes = int(mem_per_cpu) * int(cpus_per_task) * 1024 * 1024

        if slurm_limit_bytes is not None:
            mem_available = (
                slurm_limit_bytes
                if mem_available is None
                else min(mem_available, slurm_limit_bytes)
            )
    except Exception:
        pass

    return mem_available


def recommend_polars_max_threads(
    max_cores: int,
    available_memory_bytes: int | None = None,
    reserve_gb: float = 4.0,
    gb_per_thread: float = 4.0,
    filter_fraction: float | None = None,
) -> int:
    """Recommend a conservative Polars thread cap from available memory.

    Heuristic:
    - keep `reserve_gb` for Python runtime/model state
    - budget about `gb_per_thread` per parser worker
    - always return in [1, max_cores]
    """
    cap = max(1, int(max_cores))
    if available_memory_bytes is None:
        available_memory_bytes = _available_memory_bytes()

    if available_memory_bytes is None:
        # Unknown memory budget: default to a conservative small pool.
        return min(cap, 2)

    available_gb = available_memory_bytes / (1024.0**3)
    budget_gb = max(0.0, available_gb - float(reserve_gb))
    if budget_gb <= 0:
        return 1

    effective_gb_per_thread = float(gb_per_thread)
    if filter_fraction is not None:
        # If a strong subset filter is requested, per-thread memory pressure usually
        # decreases; allow a larger thread budget while keeping a conservative floor.
        frac = max(0.1, min(1.0, float(filter_fraction)))
        effective_gb_per_thread = max(1.0, effective_gb_per_thread * frac)

    memory_limited_threads = max(1, int(budget_gb // effective_gb_per_thread))
    return max(1, min(cap, memory_limited_threads))


def _estimate_filter_fraction_csv(
    path: str,
    filter_event_type: str | None,
    filter_transcript_filter: str | None,
    filter_variability: str | None,
    sample_rows: int = 500000,
) -> float | None:
    """Estimate row retention fraction for CLI subset filters from a CSV sample."""
    has_filters = any(
        x is not None for x in [filter_event_type, filter_transcript_filter]
    ) or (filter_variability not in {None, "both"})
    if not has_filters:
        return 1.0

    usecols = ["Event Type", "transcript_filter", "Variability"]
    seen = 0
    matched = 0
    try:
        for chunk in pd.read_csv(
            path,
            compression="infer",
            usecols=usecols,
            chunksize=100000,
            low_memory=False,
        ):
            mask = pd.Series(True, index=chunk.index)
            if filter_event_type is not None:
                mask &= chunk["Event Type"] == filter_event_type
            if filter_transcript_filter is not None:
                mask &= chunk["transcript_filter"] == filter_transcript_filter
            if filter_variability is not None and filter_variability != "both":
                mask &= chunk["Variability"] == filter_variability

            seen += int(chunk.shape[0])
            matched += int(mask.sum())
            if seen >= sample_rows:
                break

        if seen == 0:
            return None
        return max(0.0, min(1.0, matched / seen))
    except Exception:
        return None


def estimate_filter_fraction_for_path(
    path: str,
    filter_event_type: str | None,
    filter_transcript_filter: str | None,
    filter_variability: str | None,
    sample_rows: int = 500000,
) -> float | None:
    """Public wrapper for estimating row retention under CLI subset filters."""
    return _estimate_filter_fraction_csv(
        path,
        filter_event_type=filter_event_type,
        filter_transcript_filter=filter_transcript_filter,
        filter_variability=filter_variability,
        sample_rows=sample_rows,
    )


def recommend_auto_data_reader_backend(
    path: str,
    max_cores: int,
    filter_event_type: str | None = None,
    filter_transcript_filter: str | None = None,
    filter_variability: str | None = None,
) -> tuple[str, str]:
    """Recommend auto backend with a memory safety gate for polars->pandas.

    Returns
    -------
    tuple
        (backend, reason)
    """
    if not _is_module_available("polars"):
        return ("pandas", "polars not installed")

    available = _available_memory_bytes()
    if available is None:
        return ("polars", "memory budget unknown")

    try:
        file_size = os.path.getsize(path)
    except Exception:
        file_size = None

    if file_size is None:
        return ("polars", "input size unknown")

    # Conservative peak estimate for compressed CSV -> polars DataFrame -> pandas DataFrame.
    # Keeping a fixed reserve helps avoid OOM-kill from transient allocator spikes.
    base_expand = 10.0
    if path.endswith(".gz"):
        base_expand = 12.0
    fraction = _estimate_filter_fraction_csv(
        path,
        filter_event_type=filter_event_type,
        filter_transcript_filter=filter_transcript_filter,
        filter_variability=filter_variability,
    )

    effective_size = file_size
    fraction_msg = ""
    if fraction is not None:
        # Keep a floor so estimates do not become unrealistically optimistic.
        effective_size = int(file_size * max(0.01, float(fraction)))
        fraction_msg = f", sampled_filter_fraction={fraction:.3f}"

    estimated_polars_plus_pandas = int(effective_size * base_expand * 2.0)
    reserve_bytes = 2 * 1024**3

    if available < (estimated_polars_plus_pandas + reserve_bytes):
        avail_gb = available / (1024.0**3)
        need_gb = (estimated_polars_plus_pandas + reserve_bytes) / (1024.0**3)
        return (
            "pandas",
            (
                "memory guard: "
                f"available={avail_gb:.1f}GB < estimated_polars_peak={need_gb:.1f}GB"
                f"{fraction_msg}"
            ),
        )

    polars_threads = recommend_polars_max_threads(
        max_cores=max_cores,
        available_memory_bytes=available,
        filter_fraction=fraction,
    )
    return (
        "polars",
        (
            "memory guard passed: "
            f"available={(available / (1024.0 ** 3)):.1f}GB, polars_threads={polars_threads}"
            f"{fraction_msg}"
        ),
    )


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


def _apply_row_filters(
    df: pd.DataFrame,
    filter_event_type: str | None,
    filter_transcript_filter: str | None,
    filter_variability: str | None,
) -> pd.DataFrame:
    """Apply optional row-level subset filters used by CLI single-config runs."""
    mask = pd.Series(True, index=df.index)
    if filter_event_type is not None:
        mask &= df["Event Type"] == filter_event_type
    if filter_transcript_filter is not None:
        mask &= df["transcript_filter"] == filter_transcript_filter
    if filter_variability is not None and filter_variability != "both":
        mask &= df["Variability"] == filter_variability
    return df.loc[mask]


def _load_dataset_pandas(
    path: str,
    verbose: bool = False,
    filter_event_type: str | None = None,
    filter_transcript_filter: str | None = None,
    filter_variability: str | None = None,
) -> pd.DataFrame:
    """Read CSV using pandas with a pyarrow-first strategy."""
    t0 = time.perf_counter()
    has_filters = any(
        x is not None for x in [filter_event_type, filter_transcript_filter]
    ) or (filter_variability not in {None, "both"})

    if has_filters:
        # Chunked parsing keeps peak memory bounded for very large datasets.
        chunks: list[pd.DataFrame] = []
        first_columns: list[str] | None = None
        for chunk in pd.read_csv(
            path,
            compression="infer",
            low_memory=False,
            chunksize=100000,
        ):
            if first_columns is None:
                first_columns = list(chunk.columns)
            filtered = _apply_row_filters(
                chunk,
                filter_event_type=filter_event_type,
                filter_transcript_filter=filter_transcript_filter,
                filter_variability=filter_variability,
            )
            if not filtered.empty:
                chunks.append(filtered)

        if chunks:
            df = pd.concat(chunks, ignore_index=True)
        else:
            df = pd.DataFrame(columns=(first_columns or []))
        dt = time.perf_counter() - t0
        vlog(
            verbose,
            (
                "Loaded with pandas(chunked-c) in "
                f"{dt:.2f}s after row filters "
                f"event_type={filter_event_type}, "
                f"transcript_filter={filter_transcript_filter}, "
                f"variability={filter_variability}"
            ),
        )
        return df

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


def _load_dataset_polars(
    path: str,
    verbose: bool = False,
    filter_event_type: str | None = None,
    filter_transcript_filter: str | None = None,
    filter_variability: str | None = None,
) -> pd.DataFrame:
    """Read CSV with polars and convert to pandas for sklearn compatibility."""
    t0 = time.perf_counter()
    try:
        import polars as pl
    except Exception as exc:
        raise RuntimeError(
            "polars backend requested but polars is not installed"
        ) from exc

    has_filters = any(
        x is not None for x in [filter_event_type, filter_transcript_filter]
    ) or (filter_variability not in {None, "both"})

    if has_filters:
        # Lazy CSV scan keeps memory lower by pushing row predicates into parsing.
        lf = pl.scan_csv(path, infer_schema_length=10000, try_parse_dates=False)
        if filter_event_type is not None:
            lf = lf.filter(pl.col("Event Type") == filter_event_type)
        if filter_transcript_filter is not None:
            lf = lf.filter(pl.col("transcript_filter") == filter_transcript_filter)
        if filter_variability is not None and filter_variability != "both":
            lf = lf.filter(pl.col("Variability") == filter_variability)
        df_pl = lf.collect()
    else:
        df_pl = pl.read_csv(path, infer_schema_length=10000, try_parse_dates=False)

    df = df_pl.to_pandas()
    dt = time.perf_counter() - t0
    if has_filters:
        vlog(
            verbose,
            (
                "Loaded with polars(lazy-filtered)->pandas in "
                f"{dt:.2f}s after row filters "
                f"event_type={filter_event_type}, "
                f"transcript_filter={filter_transcript_filter}, "
                f"variability={filter_variability}"
            ),
        )
    else:
        vlog(verbose, f"Loaded with polars->pandas in {dt:.2f}s")
    return df


def load_dataset(
    path: str,
    verbose: bool = False,
    reader_backend: str = "auto",
    filter_event_type: str | None = None,
    filter_transcript_filter: str | None = None,
    filter_variability: str | None = None,
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
        df = _load_dataset_pandas(
            path,
            verbose=verbose,
            filter_event_type=filter_event_type,
            filter_transcript_filter=filter_transcript_filter,
            filter_variability=filter_variability,
        )
    elif reader_backend == "polars":
        df = _load_dataset_polars(
            path,
            verbose=verbose,
            filter_event_type=filter_event_type,
            filter_transcript_filter=filter_transcript_filter,
            filter_variability=filter_variability,
        )
    else:
        try:
            df = _load_dataset_polars(
                path,
                verbose=verbose,
                filter_event_type=filter_event_type,
                filter_transcript_filter=filter_transcript_filter,
                filter_variability=filter_variability,
            )
        except Exception as exc:
            vlog(verbose, f"Auto backend fallback to pandas due to: {exc}")
            df = _load_dataset_pandas(
                path,
                verbose=verbose,
                filter_event_type=filter_event_type,
                filter_transcript_filter=filter_transcript_filter,
                filter_variability=filter_variability,
            )

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
        y_raw = df["PSI"].astype(float)
        keep = (y_raw > low_thr) & (y_raw < high_thr)
        df = df.loc[keep].copy()
        y = df["PSI"].astype(float).to_numpy()
        if use_logit:
            y = _logit_transform(y)
            vlog(
                verbose,
                f"Built regression target (logit-scale) with thresholds=({low_thr:.4f}, {high_thr:.4f}), "
                f"kept_n={y.size}",
            )
        else:
            vlog(
                verbose,
                f"Built regression target with thresholds=({low_thr:.4f}, {high_thr:.4f}), "
                f"kept_n={y.size}",
            )
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
