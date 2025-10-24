from .Airfoil import Airfoil
import math
import numpy as np
from typing import List, Tuple

# Require SymPy via our utilities; FlatPlate always drives from symbolic expressions
from symbolic.utils import sym, eval_numeric

class FlatPlate(Airfoil):
    def __init__(self, chord: float, pivot: float = 0.25):
        """Simple FlatPlate geometry container."""

        self.chord = {
            'value': float(chord),
            'symbol': sym.Symbol('c', real=True),
        }

        self.pivot = {
            'value': float(pivot),
            'symbol': sym.Symbol('p', real=True),
        }

        self.alpha = {
            'value': 0.0,
            'symbol': sym.Symbol('a', real=True),
        }

        self.x: List[float] = []
        self.y: List[float] = []

    def symbolic_endpoints(self, c=None, a=None, p=None) -> Tuple[Tuple[object, object], Tuple[object, object]]:
        """Return SymPy expressions for endpoints rotated about a pivot.

        Inputs
        - c: chord length (SymPy symbol or numeric)
        - a: angle of attack in radians (SymPy symbol or numeric)
        - p: pivot location along x; defaults to c/4

        Output
        - ((x0, y0), (x1, y1)) expressions (always SymPy types)
        """
        c_sym = sym.sympify(c) if c is not None else self.chord['symbol']
        a_sym = sym.sympify(a) if a is not None else self.alpha['symbol']
        p_sym = sym.sympify(p) if p is not None else self.pivot['symbol']

        # Calculate locations of rotated endpoints about pivot p
        # Convention: positive a_sym = CLOCKWISE rotation
        # => use cos(a) unchanged, flip sign on sin terms compared to CCW
        x0 = p_sym + (0 - p_sym) * sym.cos(a_sym)
        y0 = -(0 - p_sym) * sym.sin(a_sym)
        x1 = p_sym + (c_sym - p_sym) * sym.cos(a_sym)
        y1 = -(c_sym - p_sym) * sym.sin(a_sym)
        return (x0, y0), (x1, y1)

    def compute_endpoint_values(self):
        """Compute endpoints (rotated about quarter-chord) and assign self.x, self.y."""
        assert self.chord['value'] > 0.0, "Chord must be positive."
        assert -math.pi/2 <= self.alpha['value'] <= math.pi/2, (
            "Alpha must be between -90 and 90 degrees in radians."
        )

        (x0_expr, y0_expr), (x1_expr, y1_expr) = self.symbolic_endpoints(
            c=None, a=None, p=None
        )
        subs = {
            self.chord['symbol']: self.chord['value'],
            self.alpha['symbol']: self.alpha['value'],
            self.pivot['symbol']: self.pivot['value'],
        }

        x0 = eval_numeric(x0_expr, subs)
        y0 = eval_numeric(y0_expr, subs)
        x1 = eval_numeric(x1_expr, subs)
        y1 = eval_numeric(y1_expr, subs)

        self.x = [x0, x1]
        self.y = [y0, y1]
        return

    def sample(self, n: int) -> Tuple[List[float], List[float]]:
        """Return n points along the straight plate between endpoints.

        Requires endpoints to be available via orient_to_alpha().
        """
        assert self.x and self.y, (
            "Call orient_to_alpha() first to generate endpoints."
        )
        assert n >= 2, "n must be >= 2"
        # Delegate to parametric sampling for consistency
        xs, ys = self.sample_points(n)
        return xs, ys

    def param_line_expr(self, c=None, a=None, p=None):
        """Return a parametric line between endpoints as expressions x(s), y(s).

        s in [0, 1]. If c/a/p are None, uses the instance's symbols.
        Returns: (s_symbol, x_of_s_expr, y_of_s_expr)
        """
        s = sym.Symbol('s', real=True)
        (x0, y0), (x1, y1) = self.symbolic_endpoints(c=c, a=a, p=p)
        x_s = x0 + s * (x1 - x0)
        y_s = y0 + s * (y1 - y0)
        return s, x_s, y_s

    def sample_points(self, n: int) -> Tuple[List[float], List[float]]:
        """Sample the parametric line using lambdified callables for speed.

        Requires endpoint expressions and substitutes the current numeric
        values for chord, alpha, and pivot.
        """
        assert n >= 2, "n must be >= 2"
        s_sym, x_s_expr, y_s_expr = self.param_line_expr(c=None, a=None, p=None)
        # Build fast callables with a stable arg order
        args = [s_sym, self.chord['symbol'], self.alpha['symbol'], self.pivot['symbol']]
        fx = to_callable(x_s_expr, args)
        fy = to_callable(y_s_expr, args)

        s_vals = np.linspace(0.0, 1.0, n)
        c_val = self.chord['value']
        a_val = self.alpha['value']
        p_val = self.pivot['value']

        xs_arr = fx(s_vals, c_val, a_val, p_val)
        ys_arr = fy(s_vals, c_val, a_val, p_val)

        return list(np.asarray(xs_arr, dtype=float)), list(np.asarray(ys_arr, dtype=float))
