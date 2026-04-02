"""Centralised thread-local storage for cross-module parameter passing.

These contexts allow passing data across sklearn's Pipeline abstraction
without modifying sklearn internals. Each is a threading.local() storage
object that is thread-safe (each worker thread has its own copy).
"""

import threading

__all__ = [
    "_mlp_val_indices_context",
    "_es_raw_val_context",
    "_es_pp_val_context",
]


# Legacy: MLP val-indices context (kept for backward compat; will be unused after fix).
# Set by _ContextAwareCV before each inner fold fit call.
_mlp_val_indices_context: threading.local = threading.local()

# Raw (non-preprocessed) val data set by _ContextAwareCV before each inner fold.
# Preprocessor has not yet been fitted at this point.
_es_raw_val_context: threading.local = threading.local()

# Preprocessed val data set by _ESPipeline after fitting the preprocessor.
# MLP reads from this to get the val set already transformed to the model's input space.
_es_pp_val_context: threading.local = threading.local()
