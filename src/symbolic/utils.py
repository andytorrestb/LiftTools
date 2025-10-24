"""
Symbolic utilities for LiftTools.

Provides a central SymPy import alias and helpers to interoperate between
symbolic expressions and fast numeric evaluation.

Exports
- sym: SymPy module alias
- symbols, Matrix, sin, cos, pi, lambdify, simplify: convenience re-exports
- ensure_symbol(value, name): Symbol-or-value normalizer
- eval_numeric(expr, subs): Evaluate a scalar SymPy expression to float
- to_callable(expr, args): Cached NumPy-backed callable for fast evaluation
"""
from __future__ import annotations

from collections import OrderedDict
from typing import Any, Dict, Hashable, Iterable, Mapping, Sequence, Tuple, Union, Callable

import sympy as sym
from sympy import symbols, Matrix, sin, cos, pi, lambdify, simplify
from sympy.core.basic import Basic
from numbers import Number

# Optional: used only if caller passes non-scalar to eval_numeric in the future
# import numpy as np

__all__ = [
    "sym",
    "symbols",
    "Matrix",
    "sin",
    "cos",
    "pi",
    "lambdify",
    "simplify",
    "ensure_symbol",
    "eval_numeric",
    "to_callable",
]


def ensure_symbol(value: Any, name: str) -> Union[Basic, Number, Any]:
    """Return a SymPy Symbol if value is None; otherwise return the input as-is.

    Notes
    - If you pass an existing SymPy expression, it is returned unchanged.
    - If you pass a numeric value (int/float), it is returned unchanged.
    - If you pass a string and value is not None, the string is returned as-is.
      Prefer passing None to create a Symbol instead.
    """
    if value is None:
        # Create a real-valued symbol for common numeric use
        return sym.Symbol(name, real=True)
    # Keep SymPy expressions untouched
    if isinstance(value, Basic):
        return value
    # Return numerics (or any other object) as provided
    return value


def _normalize_subs(subs: Mapping[Union[str, Basic], Any] | None) -> Dict[Basic, Any]:
    """Normalize substitution mapping keys to SymPy Symbols/expressions."""
    if not subs:
        return {}
    norm: Dict[Basic, Any] = {}
    for k, v in subs.items():
        if isinstance(k, Basic):
            key = k
        else:
            # Treat string keys as Symbols; other types -> try sympify
            key = sym.Symbol(str(k)) if isinstance(k, str) else sym.sympify(k)
        norm[key] = v
    return norm


def eval_numeric(expr: Any, subs: Mapping[Union[str, Basic], Any] | None = None) -> float:
    """Evaluate a SymPy scalar expression to a Python float.

    Parameters
    - expr: SymPy expression or numeric value
    - subs: dict mapping Symbols (or their names) to numeric values

    Returns
    - float value of the expression after substitution

    Raises
    - TypeError if the expression is not scalar-coercible to float
    """
    if isinstance(expr, (int, float)):
        return float(expr)
    if not isinstance(expr, Basic):
        # Try to coerce non-SymPy inputs into a SymPy object
        expr = sym.sympify(expr)
    subs_norm = _normalize_subs(subs)
    evaluated = expr.subs(subs_norm)
    # Use evalf to get a numeric value, then float()
    try:
        return float(evaluated.evalf())
    except TypeError as e:
        raise TypeError(
            "eval_numeric expects a scalar expression coercible to float"
        ) from e


def _expr_key(expr: Any) -> Hashable:
    # srepr tends to be a stable canonical form for expressions
    if isinstance(expr, Basic):
        return sym.srepr(expr)
    # For non-SymPy inputs, sympify first
    return sym.srepr(sym.sympify(expr))


def _arg_key(arg: Any) -> Hashable:
    if isinstance(arg, Basic):
        return sym.srepr(arg)
    if isinstance(arg, str):
        return f"Symbol('{arg}')"
    # Fallback to stringified sympify representation
    return sym.srepr(sym.sympify(arg))


def _to_symbol(a: Any) -> Basic:
    if isinstance(a, Basic):
        return a
    if isinstance(a, str):
        return sym.Symbol(a)
    # sympify as last resort (e.g., sympy function of Symbol)
    return sym.sympify(a)
