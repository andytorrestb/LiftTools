from .Airfoil import Airfoil
import math
from typing import List, Tuple

# Require SymPy via our utilities; FlatPlate always drives from symbolic expressions
from symbolic.utils import sym, eval_numeric


class FlatPlate(Airfoil):
    def __init__(self, chord: float, alpha_rad: float = 0.0):
        """Simple FlatPlate geometry container."""
        self.chord = float(chord)
        self.alpha_rad = float(alpha_rad)
        self.x: List[float] = []
        self.y: List[float] = []

    def set_alpha(self, alpha_rad: float):
        """Set angle of attack in radians."""
        self.alpha_rad = float(alpha_rad)

    @staticmethod
    def symbolic_endpoints(c=None, a=None, p=None) -> Tuple[Tuple[object, object], Tuple[object, object]]:
        """Return SymPy expressions for endpoints rotated about a pivot.

        Inputs
        - c: chord length (SymPy symbol or numeric)
        - a: angle of attack in radians (SymPy symbol or numeric)
        - p: pivot location along x; defaults to c/4

        Output
        - ((x0, y0), (x1, y1)) expressions (always SymPy types)
        """
        c_sym = sym.sympify(c) if c is not None else sym.Symbol('c', real=True)
        a_sym = sym.sympify(a) if a is not None else sym.Symbol('a', real=True)
        p_sym = sym.sympify(p) if p is not None else c_sym/4

        x0 = p_sym + (0 - p_sym)*sym.cos(a_sym)
        y0 = (0 - p_sym)*sym.sin(a_sym)
        x1 = p_sym + (c_sym - p_sym)*sym.cos(a_sym)
        y1 = (c_sym - p_sym)*sym.sin(a_sym)
        return (x0, y0), (x1, y1)

    def orient_to_alpha(self):
        """Compute endpoints (rotated about quarter-chord) and assign self.x, self.y."""
        assert self.chord > 0.0, "Chord must be positive."
        assert -math.pi/2 <= self.alpha_rad <= math.pi/2, (
            "Alpha must be between -90 and 90 degrees in radians."
        )

        (x0_expr, y0_expr), (x1_expr, y1_expr) = self.symbolic_endpoints(
            c=self.chord, a=self.alpha_rad, p=None
        )

        x0 = eval_numeric(x0_expr)
        y0 = eval_numeric(y0_expr)
        x1 = eval_numeric(x1_expr)
        y1 = eval_numeric(y1_expr)

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
        x0, x1 = self.x
        y0, y1 = self.y
        xs: List[float] = []
        ys: List[float] = []
        for i in range(n):
            s = i/(n-1)
            xs.append(x0 + s*(x1 - x0))
            ys.append(y0 + s*(y1 - y0))
        return xs, ys
