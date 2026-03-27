"""Tests for memory-friendly filtered pandas loading."""

from pathlib import Path

import pandas as pd

from splicing_ml.data import load_dataset


def _write_test_csv(path: Path) -> None:
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


def test_load_dataset_pandas_applies_filters(tmp_path) -> None:
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
