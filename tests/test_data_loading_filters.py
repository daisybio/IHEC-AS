"""Tests for memory-friendly filtered pandas loading."""

from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from splicing_ml.data import _load_dataset_polars, generate_subset_configs, load_dataset


def _write_test_csv(path: Path) -> None:
    """Internal helper for write test csv."""
    rows = [
        {
            "PSI": 0.1,
            "Event Type": "RI",
            "transcript_filter": "transcripts",
            "Variability": "Low",
            "seqnames": "chr1",
            "ontology": "A",
            "feature1": 1.0,
        },
        {
            "PSI": 0.9,
            "Event Type": "RI",
            "transcript_filter": "transcripts",
            "Variability": "High",
            "seqnames": "chr1",
            "ontology": "A",
            "feature1": 2.0,
        },
        {
            "PSI": 0.2,
            "Event Type": "SE",
            "transcript_filter": "transcripts",
            "Variability": "Low",
            "seqnames": "chr2",
            "ontology": "B",
            "feature1": 3.0,
        },
    ]
    pd.DataFrame(rows).to_csv(path, index=False)


def _normalize_for_compare(df: pd.DataFrame) -> pd.DataFrame:
    """Sort deterministically to compare row-equivalent dataframes."""
    out = df.copy()
    sort_cols = list(out.columns)
    if sort_cols:
        out = out.sort_values(by=sort_cols, kind="mergesort", na_position="last")
    return out.reset_index(drop=True)


def test_load_dataset_pandas_applies_filters(tmp_path) -> None:
    """Test load dataset pandas applies filters."""
    csv_path = tmp_path / "toy.csv"
    _write_test_csv(csv_path)

    df = load_dataset(
        str(csv_path),
        reader_backend="pandas",
        filter_event_type="RI",
        filter_transcript_filter="transcripts",
        filter_variability="Low",
    )

    assert df.shape[0] == 1
    assert set(df["Event Type"].unique()) == {"RI"}
    assert set(df["Variability"].unique()) == {"Low"}


def test_load_dataset_pandas_variability_both_keeps_both_levels(tmp_path) -> None:
    """Test load dataset pandas variability both keeps both levels."""
    csv_path = tmp_path / "toy.csv"
    _write_test_csv(csv_path)

    df = load_dataset(
        str(csv_path),
        reader_backend="pandas",
        filter_event_type="RI",
        filter_transcript_filter="transcripts",
        filter_variability="both",
    )

    assert df.shape[0] == 2
    assert set(df["Variability"].unique()) == {"Low", "High"}


def test_load_dataset_polars_receives_filters(monkeypatch, tmp_path) -> None:
    """Test load dataset polars receives filters."""
    csv_path = tmp_path / "toy.csv"
    _write_test_csv(csv_path)

    captured: dict[str, str | None] = {
        "event": None,
        "transcript": None,
        "variability": None,
    }

    def fake_polars_loader(
        path: str,
        verbose: bool = False,
        filter_event_type: str | None = None,
        filter_transcript_filter: str | None = None,
        filter_variability: str | None = None,
    ) -> pd.DataFrame:
        """Fake polars loader."""
        captured["event"] = filter_event_type
        captured["transcript"] = filter_transcript_filter
        captured["variability"] = filter_variability
        return pd.DataFrame(
            {
                "PSI": [0.1],
                "Event Type": ["RI"],
                "transcript_filter": ["transcripts"],
                "Variability": ["Low"],
                "seqnames": ["chr1"],
                "ontology": ["A"],
            }
        )

    monkeypatch.setattr("splicing_ml.data._load_dataset_polars", fake_polars_loader)

    load_dataset(
        str(csv_path),
        reader_backend="polars",
        filter_event_type="RI",
        filter_transcript_filter="transcripts",
        filter_variability="Low",
    )

    assert captured == {
        "event": "RI",
        "transcript": "transcripts",
        "variability": "Low",
    }


def test_full_dataset_polars_and_pandas_filtering_match_for_all_subsets() -> None:
    """Test full dataset polars and pandas filtering match for all subsets."""
    pytest.importorskip("polars")

    project_root = Path(__file__).resolve().parents[1]
    data_path = project_root / "processed_data" / "aggregated_dt_filtered.csv.gz"
    if not data_path.exists():
        pytest.skip(f"Full dataset not found: {data_path}")

    # Load once with pandas and derive all subset combinations from actual data.
    df_full = load_dataset(str(data_path), reader_backend="pandas")
    configs = generate_subset_configs(df_full)
    combos = sorted(
        {(cfg.event_type, cfg.transcript_filter, cfg.variability) for cfg in configs}
    )
    assert combos, "No subset combinations were generated from full dataset"

    for event_type, transcript_filter, variability in combos:
        expected_mask = (df_full["Event Type"] == event_type) & (
            df_full["transcript_filter"] == transcript_filter
        )
        if variability != "both":
            expected_mask &= df_full["Variability"] == variability
        expected = df_full.loc[expected_mask].copy()

        actual = _load_dataset_polars(
            str(data_path),
            filter_event_type=event_type,
            filter_transcript_filter=transcript_filter,
            filter_variability=variability,
        )

        common_cols = sorted(set(expected.columns) & set(actual.columns))
        assert common_cols, (
            "No common columns to compare for combo "
            f"event_type={event_type}, transcript_filter={transcript_filter}, "
            f"variability={variability}"
        )
        expected_cmp = _normalize_for_compare(expected[common_cols])
        actual_cmp = _normalize_for_compare(actual[common_cols])

        assert_frame_equal(
            expected_cmp,
            actual_cmp,
            check_dtype=False,
            obj=(
                "pandas(full->in-memory-filter) vs polars(filtered-load) mismatch for "
                f"event_type={event_type}, transcript_filter={transcript_filter}, "
                f"variability={variability}, expected_rows={len(expected_cmp)}, "
                f"actual_rows={len(actual_cmp)}"
            ),
        )
