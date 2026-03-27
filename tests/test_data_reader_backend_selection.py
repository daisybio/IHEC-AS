"""Tests for memory-aware auto data-reader backend selection."""

from splicing_ml.data import recommend_auto_data_reader_backend


def test_auto_backend_uses_pandas_when_memory_tight(monkeypatch, tmp_path) -> None:
    data_path = tmp_path / "big.csv.gz"
    data_path.write_bytes(b"x")

    monkeypatch.setattr("splicing_ml.data._is_module_available", lambda name: True)
    monkeypatch.setattr(
        "splicing_ml.data._available_memory_bytes", lambda: 8 * (1024**3)
    )
    monkeypatch.setattr("os.path.getsize", lambda _: int(1.4 * (1024**3)))

    backend, reason = recommend_auto_data_reader_backend(str(data_path), max_cores=8)
    assert backend == "pandas"
    assert "memory guard" in reason


def test_auto_backend_uses_polars_when_memory_sufficient(monkeypatch, tmp_path) -> None:
    data_path = tmp_path / "small.csv.gz"
    data_path.write_bytes(b"x")

    monkeypatch.setattr("splicing_ml.data._is_module_available", lambda name: True)
    monkeypatch.setattr(
        "splicing_ml.data._available_memory_bytes", lambda: 128 * (1024**3)
    )
    monkeypatch.setattr("os.path.getsize", lambda _: int(0.2 * (1024**3)))

    backend, reason = recommend_auto_data_reader_backend(str(data_path), max_cores=8)
    assert backend == "polars"
    assert "memory guard passed" in reason


def test_auto_backend_without_polars_falls_back_to_pandas(
    monkeypatch, tmp_path
) -> None:
    data_path = tmp_path / "data.csv.gz"
    data_path.write_bytes(b"x")

    monkeypatch.setattr("splicing_ml.data._is_module_available", lambda name: False)

    backend, reason = recommend_auto_data_reader_backend(str(data_path), max_cores=8)
    assert backend == "pandas"
    assert reason == "polars not installed"


def test_auto_backend_filter_aware_can_select_polars(monkeypatch, tmp_path) -> None:
    data_path = tmp_path / "big.csv.gz"
    data_path.write_bytes(b"x")

    monkeypatch.setattr("splicing_ml.data._is_module_available", lambda name: True)
    monkeypatch.setattr(
        "splicing_ml.data._available_memory_bytes", lambda: 16 * (1024**3)
    )
    monkeypatch.setattr("os.path.getsize", lambda _: int(1.4 * (1024**3)))
    monkeypatch.setattr(
        "splicing_ml.data._estimate_filter_fraction_csv", lambda *args, **kwargs: 0.02
    )

    backend, reason = recommend_auto_data_reader_backend(
        str(data_path),
        max_cores=8,
        filter_event_type="RI",
        filter_transcript_filter="transcripts",
        filter_variability="Low",
    )
    assert backend == "polars"
    assert "sampled_filter_fraction=0.020" in reason
