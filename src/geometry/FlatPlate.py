from .Airfoil import Airfoil
import math
import numpy as np
from typing import List, Tuple

# Require SymPy via our utilities; FlatPlate always drives from symbolic expressions
from symbolic.utils import sym, eval_numeric, to_callable

def rotate_point(point: Tuple[float, float], pivot: Tuple[float, float], angle_rad: float) -> Tuple[float, float]:
    """Rotate a point around a pivot by a given angle in radians.

    Parameters
    - point: (x, y) coordinates of the point to rotate
    - pivot: (x, y) coordinates of the pivot point
    - angle_rad: rotation angle in radians (positive = clockwise)

    Returns
    - (x_rotated, y_rotated): coordinates of the rotated point
    """
    x1 = point[0] - pivot[0]
    y1 = point[1] - pivot[1]

    angle_rad = -angle_rad  # Convert to counter-clockwise for standard rotation
    x_rotated = x1 * math.cos(angle_rad) + y1 * math.sin(angle_rad)
    y_rotated = x1 * math.sin(angle_rad) - y1 * math.cos(angle_rad)

    x_final = x_rotated + pivot[0]
    y_final = y_rotated + pivot[1]
    return (x_final, y_final)

class FlatPlate(Airfoil):
    # Optional flap parameters (may be set externally by models)
    flap_deflection_rad: float = 0.0  # positive = clockwise per class convention
    flap_length_le: float = 1.0       # fraction of chord measured from LE to hinge

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
        # Convention here: positive a_sym = CLOCKWISE rotation
        # This matches convention in aerodynamics.
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
    
    def symbolic_line_expr(self, c=None, a=None, p=None):
        """Return the line in slope-intercept form y = m*x + b for the rotated plate.

        Uses symbolic endpoints to derive:
        - m = (y1 - y0) / (x1 - x0) = tan(alpha)
        - b = y0 - m*x0 = -p * tan(alpha)

        Note: Not defined for vertical orientation (alpha = ±pi/2).

        If c/a/p are None, uses the instance's symbols.
        Returns: (x_symbol, y_of_x_expr, m_expr, b_expr)
        """
        # Use instance symbols if not provided
        c_sym = sym.sympify(c) if c is not None else self.chord['symbol']
        a_sym = sym.sympify(a) if a is not None else self.alpha['symbol']
        p_sym = sym.sympify(p) if p is not None else self.pivot['symbol']

        # Derive slope and intercept from the symbolic endpoints
        (x0, y0), (x1, y1) = self.symbolic_endpoints(c=c_sym, a=a_sym, p=p_sym)
        m = sym.simplify((y1 - y0) / (x1 - x0))
        b = sym.simplify(y0 - m * x0)

        # y(x) = m*x + b
        x = sym.Symbol('x', real=True)
        return sym.simplify(m * x + b)

    # --- Flapped geometry helpers -----------------------------------------
    def symbolic_hinge_point(self, c=None, a=None, p=None, l_le=None):
        """Return SymPy expressions (x_h, y_h) for the flap hinge point.

        Hinge is located at x = l_le * c on the unrotated chord and then
        rotated about pivot p by angle a (clockwise-positive convention).
        """
        c_sym = sym.sympify(c) if c is not None else self.chord['symbol']
        a_sym = sym.sympify(a) if a is not None else self.alpha['symbol']
        p_sym = sym.sympify(p) if p is not None else self.pivot['symbol']
        l_sym = sym.sympify(l_le) if l_le is not None else sym.Symbol('l_le', real=True)

        x_h = p_sym + (l_sym * c_sym - p_sym) * sym.cos(a_sym)
        y_h = -(l_sym * c_sym - p_sym) * sym.sin(a_sym)
        return x_h, y_h

    def symbolic_flapped_piecewise(self, c=None, a=None, p=None, l_le=None, delta=None):
        """Return y(x) for a flapped flat plate as a SymPy Piecewise.

        Main plate: from LE to hinge follows rotation by a.
        Flap: from hinge to new TE follows rotation by (a - delta) about hinge.

        Returns a tuple:
        - x symbol
        - y_piecewise(x)
        - (x0, y0): LE point after rotation by a
        - (x_te, y_te): TE point after flap deflection
        - (x_h, y_h): Hinge point
        """
        # Symbols
        c_sym = sym.sympify(c) if c is not None else self.chord['symbol']
        a_sym = sym.sympify(a) if a is not None else self.alpha['symbol']
        p_sym = sym.sympify(p) if p is not None else self.pivot['symbol']
        l_sym = sym.sympify(l_le) if l_le is not None else sym.Symbol('l_le', real=True)
        d_sym = sym.sympify(delta) if delta is not None else sym.Symbol('delta', real=True)

        x = sym.Symbol('x', real=True)

        # Endpoints of the base (unflapped) plate rotated by a
        (x0, y0), (x1, y1) = self.symbolic_endpoints(c=c_sym, a=a_sym, p=p_sym)

        # Hinge point and flap length
        x_h, y_h = self.symbolic_hinge_point(c=c_sym, a=a_sym, p=p_sym, l_le=l_sym)
        flap_c = (1 - l_sym) * c_sym

        # Lines: main and flap, as y = m x + b
        # Using the class convention (clockwise-positive): slope = -tan(angle)
        m_main = -sym.tan(a_sym)
        b_main = sym.simplify(y0 - m_main * x0)  # must pass through LE as well
        y_main = sym.simplify(m_main * x + b_main)

        m_flap = -sym.tan(a_sym - d_sym)
        b_flap = sym.simplify(y_h - m_flap * x_h)
        y_flap = sym.simplify(m_flap * x + b_flap)

        # New TE point located flap_c from hinge at orientation theta_flap = -(a - delta)
        theta_flap = -(a_sym - d_sym)
        x_te = sym.simplify(x_h + flap_c * sym.cos(theta_flap))
        y_te = sym.simplify(y_h + flap_c * sym.sin(theta_flap))

        # Piecewise by x relative to hinge x; assumes not vertical
        y_piece = sym.Piecewise((y_main, x <= x_h), (y_flap, True))
        return x, y_piece, (x0, y0), (x_te, y_te), (x_h, y_h)

    def sample_flapped_by_x(self, n: int) -> Tuple[List[float], List[float]]:
        """Sample the flapped geometry y(x) along x between LE and new TE.

        Uses current numeric values for chord, alpha, pivot, flap_length_le,
        and flap_deflection_rad. Returns x,y lists of length n.
        """
        assert n >= 2, "n must be >= 2"

        # Build symbolic piecewise once
        x_sym, y_piece, (x0, y0), (x_te, y_te), (x_h, y_h) = self.symbolic_flapped_piecewise(
            c=None, a=None, p=None, l_le=self.flap_length_le, delta=self.flap_deflection_rad
        )

        # Lambdify y(x) with explicit arg order
        args = [
            x_sym,
            self.chord['symbol'],
            self.alpha['symbol'],
            self.pivot['symbol'],
            sym.Symbol('l_le', real=True),
            sym.Symbol('delta', real=True),
        ]
        f_y = to_callable(y_piece, args)

        # Numeric endpoints for x-range
        subs = {
            self.chord['symbol']: self.chord['value'],
            self.alpha['symbol']: self.alpha['value'],
            self.pivot['symbol']: self.pivot['value'],
            sym.Symbol('l_le', real=True): float(getattr(self, 'flap_length_le', 1.0)),
            sym.Symbol('delta', real=True): float(getattr(self, 'flap_deflection_rad', 0.0)),
        }

        x0_v = eval_numeric(x0, subs)
        xte_v = eval_numeric(x_te, subs)
        x_start = min(x0_v, xte_v)
        x_end = max(x0_v, xte_v)

        xs = np.linspace(x_start, x_end, n)
        ys = f_y(xs, subs[self.chord['symbol']], subs[self.alpha['symbol']], subs[self.pivot['symbol']], subs[sym.Symbol('l_le', real=True)], subs[sym.Symbol('delta', real=True)])

        return list(np.asarray(xs, dtype=float)), list(np.asarray(ys, dtype=float))
    

    def plot_flapped_plate(self, subplot = False) -> None:
        """Plot the flapped flat plate using Matplotlib for visual verification.

        Uses current numeric values for chord, alpha, pivot, flap_length_le,
        and flap_deflection_rad.

        Plot is produced using hard-coded solution.
        """

        assert hasattr(self, 'flap_length_le'), "Flap length not set."
        assert hasattr(self, 'flap_deflection_rad'), "Flap deflection not set."
        assert hasattr(self, 'alpha'), "Alpha not set."
        assert self.chord['value'] > 0.0, "Chord must be positive."

        # Set up x coordinates for the wing and flap.
        k = self.flap_length_le
        c = self.chord['value']
        # n = 100
        # x_wing = np.linspace(0, k*c, k*n)
        # x_flap = np.linspace(k*c, c, (1-k)*n)

        # Slope of plates (main and flap)
        m_wing = math.tan(-self.alpha['value'])
        m_flap = math.tan(-(self.alpha['value'] + self.flap_deflection_rad))

        # For wing: find y intercept and x-bounds.
        # Wing equations (straight line: y=mx+b)
        b_wing = -m_wing * self.pivot['value']
        pivot = (0.25*c, 0.0)
        x_LE_wing, y_LE_wing = 0.0, 0.0
        LE_wing = (x_LE_wing, y_LE_wing)
        LE_wing_rotated = rotate_point(LE_wing, pivot, self.alpha['value'])

        x_TE_wing, y_TE_wing = k*c, 0.0
        TE_wing = (x_TE_wing, y_TE_wing)
        TE_wing_rotated = rotate_point(TE_wing, pivot, self.alpha['value'])

        # For flapp: find y intercept and x-bounds.
        b_flap = TE_wing_rotated[1] - m_flap * TE_wing_rotated[0]

        LE_flap_rotated = TE_wing_rotated
        x_TE_flap, y_TE_flap = c, 0.0
        TE_flap = (x_TE_flap, y_TE_flap)
        total_deflection = self.alpha['value'] + self.flap_deflection_rad
        TE_flap_rotated = rotate_point(TE_flap, (k*c, 0.0), total_deflection)

        # Create x-coordinate using caculated endpints.
        n = 100
        x_wing = np.linspace(LE_wing_rotated[0], TE_wing_rotated[0], int(k*n))
        x_flap = np.linspace(TE_wing_rotated[0], TE_flap_rotated[0], int((1-k)*n))

        # Calculate y-coordinates using line equations.
        y_wing = m_wing * x_wing + b_wing
        y_flap = m_flap * x_flap + b_flap

        if subplot:
            return (x_wing, y_wing), (x_flap, y_flap)
        else:
            # Plot the flapped flat plate.
            import matplotlib.pyplot as plt 
            plt.figure(figsize=(8, 4))
            plt.plot(x_wing, y_wing, label="Main Plate")
            plt.plot(x_flap, y_flap, label="Flap")
            plt.title("Flapped Flat Plate Geometry")
            plt.xlabel("x")
            plt.ylabel("y")
            plt.legend()
            plt.grid(True)
            plt.savefig("flapped_flat_plate.png")
            plt.clf()
        return

        






        # y-intercepts

