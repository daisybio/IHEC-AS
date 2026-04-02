"""Unit tests for hyperparameter candidate builders in grids.py."""

from __future__ import annotations

import pytest

from splicing_ml.models import grids as g


def _assert_nonempty_list_dict(out: list[dict], keys: list[str]) -> None:
    """Internal helper for assert nonempty list dict."""
    assert isinstance(out, list)
    assert len(out) > 0
    for cand in out:
        assert isinstance(cand, dict)
        for key in keys:
            assert key in cand
            assert isinstance(cand[key], list)
            assert len(cand[key]) > 0


def test_build_param_candidates_dispatches_correctly(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test build param candidates dispatches correctly."""
    calls: list[str] = []

    def _fake_grid(**kwargs):
        """Internal helper for fake grid."""
        calls.append("grid")
        return [{"model": ["grid"]}]

    def _fake_lhs(**kwargs):
        """Internal helper for fake lhs."""
        calls.append("lhs")
        return [{"model": ["lhs"]}]

    monkeypatch.setattr(g, "choose_param_grid", _fake_grid)
    monkeypatch.setattr(g, "choose_param_lhs_candidates", _fake_lhs)

    out1 = g.build_param_candidates(
        model_name="linear",
        task="regression",
        n_samples=1000,
        n_features=100,
        budget=8,
        strategy="grid",
    )
    out2 = g.build_param_candidates(
        model_name="rf",
        task="regression",
        n_samples=1000,
        n_features=100,
        budget=8,
        strategy="hybrid",
    )
    out3 = g.build_param_candidates(
        model_name="linear",
        task="regression",
        n_samples=1000,
        n_features=100,
        budget=8,
        strategy="random",
    )

    assert calls == ["grid", "lhs", "lhs"]
    assert out1[0]["model"] == ["grid"]
    assert out2[0]["model"] == ["lhs"]
    assert out3[0]["model"] == ["lhs"]


@pytest.mark.parametrize("task", ["regression", "classification"])
def test_grid_linear(task: str) -> None:
    """Test grid linear."""
    out = g.choose_param_grid(
        model_name="linear",
        task=task,
        n_samples=2000,
        n_features=50,
        grid_size=9,
        cuml_use_gpu=False,
    )
    _assert_nonempty_list_dict(out, ["model"])


@pytest.mark.parametrize("task", ["regression", "classification"])
def test_grid_svm(task: str) -> None:
    """Test grid svm."""
    out = g.choose_param_grid(
        model_name="svm",
        task=task,
        n_samples=5000,
        n_features=120,
        grid_size=9,
        cuml_use_gpu=False,
    )
    _assert_nonempty_list_dict(out, ["model", "model__C", "model__gamma"])
    c_vals = out[0]["model__C"]
    gamma_vals = out[0]["model__gamma"]
    assert min(c_vals) > 0
    assert min(gamma_vals) > 0


@pytest.mark.parametrize("task", ["regression", "classification"])
def test_grid_mlp(task: str) -> None:
    """Test grid mlp."""
    out = g.choose_param_grid(
        model_name="mlp",
        task=task,
        n_samples=5000,
        n_features=120,
        grid_size=16,
    )
    _assert_nonempty_list_dict(
        out,
        [
            "model",
            "model__hidden_sizes",
            "model__dropout_rate",
            "model__learning_rate",
            "model__weight_decay",
            "model__batch_size",
        ],
    )


@pytest.mark.parametrize("task", ["regression", "classification"])
def test_grid_rf_and_xgb(task: str) -> None:
    """Test grid rf and xgb."""
    pytest.importorskip("xgboost")

    rf = g.choose_param_grid(
        model_name="rf",
        task=task,
        n_samples=5000,
        n_features=120,
        grid_size=16,
        xgb_use_gpu=False,
    )
    _assert_nonempty_list_dict(
        rf,
        [
            "model",
            "model__n_estimators",
            "model__max_depth",
            "model__colsample_bynode",
            "model__min_child_weight",
        ],
    )

    xgb = g.choose_param_grid(
        model_name="xgb",
        task=task,
        n_samples=5000,
        n_features=120,
        grid_size=7,
        xgb_use_gpu=False,
    )
    assert len(xgb) == 7
    _assert_nonempty_list_dict(
        xgb,
        [
            "model",
            "model__n_estimators",
            "model__max_depth",
            "model__learning_rate",
            "model__subsample",
            "model__colsample_bytree",
        ],
    )


@pytest.mark.parametrize("model_name", ["rf", "xgb", "svm", "mlp"])
@pytest.mark.parametrize("task", ["regression", "classification"])
def test_lhs_candidates_shape(model_name: str, task: str) -> None:
    """Test lhs candidates shape."""
    if model_name in {"rf", "xgb"}:
        pytest.importorskip("xgboost")

    out = g.choose_param_lhs_candidates(
        model_name=model_name,
        task=task,
        n_samples=6000,
        n_features=150,
        budget=6,
        xgb_use_gpu=False,
        cuml_use_gpu=False,
        scale_mode="auto",
    )
    assert len(out) == 6
    for cand in out:
        assert isinstance(cand, dict)
        for key, value in cand.items():
            assert isinstance(value, list)
            assert len(value) == 1


@pytest.mark.parametrize("task", ["regression", "classification"])
def test_lhs_linear_is_singleton(task: str) -> None:
    """Test lhs linear is singleton."""
    out = g.choose_param_lhs_candidates(
        model_name="linear",
        task=task,
        n_samples=6000,
        n_features=150,
        budget=6,
        cuml_use_gpu=False,
    )
    _assert_nonempty_list_dict(out, ["model"])
    assert len(out) == 1


@pytest.mark.parametrize("task", ["regression", "classification"])
def test_optuna_distributions_svm_mlp(task: str) -> None:
    """Test optuna distributions svm mlp."""
    optuna = pytest.importorskip("optuna")
    dist = optuna.distributions

    base_svm, d_svm = g.build_optuna_distributions(
        model_name="svm",
        task=task,
        n_samples=5000,
        n_features=120,
        cuml_use_gpu=False,
    )
    assert hasattr(base_svm, "fit")
    assert set(d_svm.keys()) == {"model__C", "model__gamma"}
    assert isinstance(d_svm["model__C"], dist.FloatDistribution)
    assert isinstance(d_svm["model__gamma"], dist.FloatDistribution)

    base_mlp, d_mlp = g.build_optuna_distributions(
        model_name="mlp",
        task=task,
        n_samples=5000,
        n_features=120,
    )
    assert hasattr(base_mlp, "fit")
    assert set(d_mlp.keys()) == {
        "model__hidden_sizes",
        "model__dropout_rate",
        "model__learning_rate",
        "model__weight_decay",
    }


@pytest.mark.parametrize("task", ["regression", "classification"])
def test_optuna_distributions_rf_xgb(task: str) -> None:
    """Test optuna distributions rf xgb."""
    optuna = pytest.importorskip("optuna")
    dist = optuna.distributions
    pytest.importorskip("xgboost")

    base_rf, d_rf = g.build_optuna_distributions(
        model_name="rf",
        task=task,
        n_samples=5000,
        n_features=120,
        xgb_use_gpu=False,
    )
    assert hasattr(base_rf, "fit")
    assert set(d_rf.keys()) == {
        "model__n_estimators",
        "model__max_depth",
        "model__colsample_bynode",
        "model__min_child_weight",
    }
    assert isinstance(d_rf["model__n_estimators"], dist.IntDistribution)

    base_xgb, d_xgb = g.build_optuna_distributions(
        model_name="xgb",
        task=task,
        n_samples=5000,
        n_features=120,
        xgb_use_gpu=False,
    )
    assert hasattr(base_xgb, "fit")
    assert set(d_xgb.keys()) == {
        "model__n_estimators",
        "model__max_depth",
        "model__learning_rate",
        "model__subsample",
        "model__colsample_bytree",
        "model__min_child_weight",
    }
    assert isinstance(d_xgb["model__learning_rate"], dist.FloatDistribution)


def test_unsupported_model_raises() -> None:
    """Test unsupported model raises."""
    with pytest.raises(ValueError):
        g.choose_param_grid(
            model_name="does-not-exist",
            task="regression",
            n_samples=100,
            n_features=10,
            grid_size=5,
        )

    with pytest.raises(ValueError):
        g.choose_param_lhs_candidates(
            model_name="does-not-exist",
            task="regression",
            n_samples=100,
            n_features=10,
            budget=5,
        )

    with pytest.raises(ValueError):
        g.build_optuna_distributions(
            model_name="does-not-exist",
            task="regression",
            n_samples=100,
            n_features=10,
        )
