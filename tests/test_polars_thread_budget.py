"""Tests for memory-aware Polars thread recommendation."""

from splicing_ml.data import recommend_polars_max_threads


def _gib(n: int) -> int:
    """Internal helper for gib."""
    return n * (1024**3)


def test_recommend_polars_threads_memory_limited() -> None:
    """Test recommend polars threads memory limited."""
    threads = recommend_polars_max_threads(
        max_cores=16,
        available_memory_bytes=_gib(20),
        reserve_gb=4.0,
        gb_per_thread=4.0,
    )
    assert threads == 4


def test_recommend_polars_threads_respects_cpu_cap() -> None:
    """Test recommend polars threads respects cpu cap."""
    threads = recommend_polars_max_threads(
        max_cores=6,
        available_memory_bytes=_gib(200),
        reserve_gb=4.0,
        gb_per_thread=4.0,
    )
    assert threads == 6


def test_recommend_polars_threads_minimum_one() -> None:
    """Test recommend polars threads minimum one."""
    threads = recommend_polars_max_threads(
        max_cores=12,
        available_memory_bytes=_gib(3),
        reserve_gb=4.0,
        gb_per_thread=4.0,
    )
    assert threads == 1


def test_recommend_polars_threads_unknown_memory_defaults_to_two(monkeypatch) -> None:
    """Test recommend polars threads unknown memory defaults to two."""
    monkeypatch.setattr("splicing_ml.data._available_memory_bytes", lambda: None)
    assert recommend_polars_max_threads(max_cores=8, available_memory_bytes=None) == 2
    assert recommend_polars_max_threads(max_cores=1, available_memory_bytes=None) == 1


def test_recommend_polars_threads_filter_fraction_can_increase_budget() -> None:
    """Test recommend polars threads filter fraction can increase budget."""
    base = recommend_polars_max_threads(
        max_cores=8,
        available_memory_bytes=_gib(16),
        reserve_gb=4.0,
        gb_per_thread=4.0,
        filter_fraction=None,
    )
    filtered = recommend_polars_max_threads(
        max_cores=8,
        available_memory_bytes=_gib(16),
        reserve_gb=4.0,
        gb_per_thread=4.0,
        filter_fraction=0.1,
    )

    assert base == 3
    assert filtered == 8
