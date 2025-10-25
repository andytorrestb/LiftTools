import math

import sympy as sp

from models.AirFoilModel import AirfoilModel

class ThinAirfoilTheory(AirfoilModel):

    def set_z_expr(self, foil_type="plate") -> None:
        if foil_type == "plate":
            # For a flat plate, use the rotated endpoint y-expression as z(x)=constant along the plate
            # Use the geometry's default symbols by calling without arguments
            x, theta = sp.symbols('x theta', real=True, positive=True)
            self.x_sym = x
            self.theta_sym = theta
            self.z_expr = self.geometry.symbolic_endpoints()[0][1]  # y0 expression
            return
        
    def alpha_0(self) -> sp.Expr:
        c = self.geometry.chord['symbol']

        x = self.x_sym
        theta = self.theta_sym
        x_sub = (c/2.0) * (1 - sp.cos(theta)
                           )
        z_expr = self.z_expr
        dzdx = sp.diff(z_expr, sp.Symbol('x', real=True))
        dzdx_theta = sp.simplify(dzdx.subs(x, x_sub))
        integrand = dzdx_theta * (sp.cos(theta) - 1)
        alpha_L0 = (1 / sp.pi) * sp.integrate(integrand, (theta, 0, sp.pi))
        return sp.simplify(sp.together(alpha_L0))
    
    def solve(self) -> dict:
        # Thin airfoil theory for a flat plate: cl = 2*pi*(alpha - alpha_0)
        # Use the model's alpha_rad directly to avoid requiring a geometry instance.
        alpha_0 = self.alpha_0()  # For flat plate, camber is zero.
        # Optional: keep z_expr if set, but don't require it
        print(f"Computed alpha_0 (radians): {alpha_0}")
        cl = 2.0 * math.pi * (self.alpha_rad - alpha_0)
        results = {
            'cl': cl,
            'alpha_rad': float(self.alpha_rad),
        }
        return results
    pass