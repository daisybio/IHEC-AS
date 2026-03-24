"""Splicing ML package.

This package exposes the high-level runner for nested cross-validated splicing
prediction workflows.
"""

# Lazy import to avoid sys.modules side effects when running as __main__.
__all__ = ["run_all_configurations"]


def __getattr__(name: str):
    """Lazy-load run_all_configurations on first access."""
    if name == "run_all_configurations":
        from .pipeline import run_all_configurations

        return run_all_configurations
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
