import math

import sympy as sp

from models.AirFoilModel import AirfoilModel
from models.WeisingersApprox import WeisingersApprox

class ThinAirfoilTheory(AirfoilModel):
    def set_naca_camber(self, m: float, p: float) -> None:
        """Set the camber line for a NACA 4-digit airfoil.

        Parameters
        - m: maximum camber as fraction of chord (e.g., 0.02 for 2%)
        - p: location of maximum camber as fraction of chord (e.g., 0.4 for 40%)
        """
        x = sp.Symbol('x', real=True, positive=True)
        c = self.geometry.chord['symbol']

        # Define piecewise camber line
        z_expr = sp.Piecewise(
            (m * (x / p**2) * (2*p - x/c), x < p*c),
            (m * ((c - x) / (1 - p)**2) * (1 + x/c - 2*p), x >= p*c)
        )
        self.z_expr = z_expr
        # Store parameters for downstream use (e.g., split integration)
        self.m = float(m)
        self.p = float(p)
        
    def alpha_0(self) -> sp.Expr:
        """Compute zero-lift angle using thin airfoil theory.

        Handles NACA-4 camber Piecewise robustly by splitting the integral
        at theta* = acos(1 - 2p), avoiding SymPy's inequality solver on cos.
        """
        c_sym = self.geometry.chord['symbol']

        x = sp.Symbol('x', real=True, positive=True)
        theta = sp.Symbol('theta', real=True, positive=True)
        x_sub = (c_sym/2.0) * (1 - sp.cos(theta))

        # If we have stored NACA parameters, integrate each branch explicitly
        if hasattr(self, 'm') and hasattr(self, 'p') and 0.0 <= self.p <= 1.0:
            m = sp.nsimplify(self.m)
            p = sp.nsimplify(self.p)
            c = c_sym

            # Branch definitions (same as set_naca_camber)
            z1 = m * (x / p**2) * (2*p - x/c)
            z2 = m * ((c - x) / (1 - p)**2) * (1 + x/c - 2*p)

            dz1dx = sp.diff(z1, x)
            dz2dx = sp.diff(z2, x)

            dz1dx_theta = sp.simplify(dz1dx.subs(x, x_sub))
            dz2dx_theta = sp.simplify(dz2dx.subs(x, x_sub))

            integrand1 = dz1dx_theta * (sp.cos(theta) - 1)
            integrand2 = dz2dx_theta * (sp.cos(theta) - 1)

            # Split angle where x = p*c => (c/2)(1 - cos th) = p*c => cos th = 1 - 2p
            theta_star = sp.acos(1 - 2*p)

            I1 = sp.integrate(integrand1, (theta, 0, theta_star))
            I2 = sp.integrate(integrand2, (theta, theta_star, sp.pi))
            alpha_L0 = (1 / sp.pi) * (I1 + I2)
            out = sp.simplify(sp.together(alpha_L0))
            # If chord symbol remains, substitute numeric chord value for a numeric result
            if c_sym in out.free_symbols:
                out = sp.simplify(out.subs(c_sym, self.geometry.chord['value']))
            return out

        # Fallback: try direct integration for non-piecewise or external z_expr
        z_expr = getattr(self, 'z_expr', None)
        if z_expr is None:
            return sp.S.Zero

        dzdx = sp.diff(z_expr, x)
        dzdx_theta = sp.simplify(dzdx.subs(x, x_sub))
        integrand = dzdx_theta * (sp.cos(theta) - 1)
        alpha_L0 = (1 / sp.pi) * sp.integrate(integrand, (theta, 0, sp.pi))
        out = sp.simplify(sp.together(alpha_L0))
        if c_sym in out.free_symbols:
            out = sp.simplify(out.subs(c_sym, self.geometry.chord['value']))
        return out
    
    def alpha_n(self, n: int) -> sp.Expr:
        """Compute the angle for the nth moment using thin airfoil theory.

        Handles NACA-4 camber Piecewise robustly by splitting the integral
        at theta* = acos(1 - 2p), avoiding SymPy's inequality solver on cos.
        """
        c_sym = self.geometry.chord['symbol']

        x = sp.Symbol('x', real=True, positive=True)
        theta = sp.Symbol('theta', real=True, positive=True)
        n = sp.Symbol('n', integer=True, positive=True)
        x_sub = (c_sym/2.0) * (1 - sp.cos(theta))

        # If we have stored NACA parameters, integrate each branch explicitly
        if hasattr(self, 'm') and hasattr(self, 'p') and 0.0 <= self.p <= 1.0:
            m = sp.nsimplify(self.m)
            p = sp.nsimplify(self.p)
            c = c_sym

            # Branch definitions (same as set_naca_camber)
            z1 = m * (x / p**2) * (2*p - x/c)
            z2 = m * ((c - x) / (1 - p)**2) * (1 + x/c - 2*p)

            dz1dx = sp.diff(z1, x)
            dz2dx = sp.diff(z2, x)

            dz1dx_theta = sp.simplify(dz1dx.subs(x, x_sub))
            dz2dx_theta = sp.simplify(dz2dx.subs(x, x_sub))

            integrand1 = dz1dx_theta * sp.cos(n*theta)
            integrand2 = dz2dx_theta * sp.cos(n*theta)

            # Split angle where x = p*c => (c/2)(1 - cos th) = p*c => cos th = 1 - 2p
            theta_star = sp.acos(1 - 2*p)

            I1 = sp.integrate(integrand1, (theta, 0, theta_star))
            I2 = sp.integrate(integrand2, (theta, theta_star, sp.pi))
            alpha_n_val = (1 / sp.pi) * (I1 + I2)
            out = sp.simplify(sp.together(alpha_n_val))
            # If chord symbol remains, substitute numeric chord value for a numeric result
            if c_sym in out.free_symbols:
                out = sp.simplify(out.subs(c_sym, self.geometry.chord['value']))
            return out

    def solve(self) -> dict:
        # Thin airfoil theory for a flat plate: cl = 2*pi*(alpha - alpha_0)
        # Use the model's alpha_rad directly to avoid requiring a geometry instance.
        alpha_0_expr = self.alpha_0()  # For flat plate, camber is zero.
        # Optional: keep z_expr if set, but don't require it
        try:
            alpha_0_val = float(alpha_0_expr)
        except Exception:
            # Fallback to numerical evaluation
            alpha_0_val = float(sp.N(alpha_0_expr))
        print(f"Computed alpha_0 (radians): {alpha_0_val}")
        cl = 2.0 * math.pi * (float(self.alpha_rad) - alpha_0_val)
        cm = (math.pi / 4.0) * (self.alpha_n(2) - self.alpha_n(1))
        results = {
            'cl': cl,
            'cm': cm,
            'alpha_rad': float(self.alpha_rad),
        }
        return results