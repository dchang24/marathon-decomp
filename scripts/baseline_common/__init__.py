"""Shared helpers for the rank-1 baseline pipeline.

Imported by the digit-prefixed script dirs ``scripts/02_model_selection/baseline``
and ``scripts/03_model_fit/baseline`` (which cannot import each other directly —
``02_...`` / ``03_...`` are not valid Python module names). Each runnable script
adds ``scripts/`` to ``sys.path`` and does ``from baseline_common import ...``.

Modules
-------
- ``slices``  : the named production ``SliceSpec``s + spec/slug/CLI helpers.
- ``inits``   : reproducible (u, v) init / perturbation schemes.
- ``fitting`` : fit-with-iter-capture, convergence tables, registry save.
"""
