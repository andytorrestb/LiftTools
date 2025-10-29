import math
import numpy as np
import matplotlib.pyplot as plt

from models.AirFoilModel import AirfoilModel

def rotate_point(point: tuple[float, float], pivot: tuple[float, float], theta: tuple[float, float]) -> tuple[float, float]:
    """Rotate a point (x, z) by angle alpha_rad around the origin."""
    x_rot = x * math.cos(alpha_rad) - z * math.sin(alpha_rad)
    z_rot = x * math.sin(alpha_rad) + z * math.cos(alpha_rad)
    return x_rot, z_rot

class WeisingersApprox(AirfoilModel):
    def set_panels(self, n_panels: int) -> None:
        assert self.geometry.z is not None, "Camber line must be set before discretizing."
        # establish local name references for important data.
        z = self.geometry.z
        c = self.geometry.chord['value']
        k = self.geometry.k

        # Determine number of panels on main wing and flap
        n_wing_panels = int(k*(n_panels+1))
        n_flap_panels = n_panels - n_wing_panels

        # Discretize and interpolate panels on main wing and flap
        wing_panel_x = np.linspace(0, k*c, n_wing_panels)
        wing_panel_z = np.interp(wing_panel_x, self.geometry.x_lower, z)

        flap_panel_x = np.linspace(k*c, c, n_flap_panels+1)
        flap_panel_z = np.interp(flap_panel_x, self.geometry.x_lower, z)

        # Set to object properties for later use.
        self.wing_panel_x = wing_panel_x
        self.wing_panel_z = wing_panel_z    
        self.flap_panel_x = flap_panel_x
        self.flap_panel_z = flap_panel_z

    def plot_panels(self) -> None:
        plt.figure(figsize=(10, 4))
        # plt.plot(self.geometry.x_lower, self.geometry.z, 'b-', label='Camber Line')
        plt.plot(self.wing_panel_x, self.wing_panel_z, 'ro--', label='Wing Panels')
        plt.plot(self.flap_panel_x, self.flap_panel_z, 'go--', label='Flap Panels')
        plt.xlabel('x (chordwise)')
        plt.ylabel('z (camber)')
        plt.ylim(-0.1, 0.1)
        plt.title('Airfoil Camber Line with Discretized Panels')
        plt.legend()
        plt.grid(True)
        plt.savefig(f"weisingers_approx_panels_naca{self.geometry.naca_code}.png")         

    def solve(self) -> dict:
        # Placeholder for Weisinger's Approximation solution method
        print("Solving using Weisinger's Approximation...")
        pass

    def set_flap(self, deflection_rad: float, length_le: float) -> None:
        # Placeholder for setting flap deflection specific to Weisinger's Approximation
        self.geometry.flap_deflection_rad = float(deflection_rad)
        self.geometry.flap_length_le = float(length_le)
        pass


    def orient_panels() -> None:
        """Update panel coordinates for the current alpha if supported.

        Subclasses that provide `compute_endpoint_values()` will have this
        method call it to refresh x and y based on the current parameters.
        """
        assert self.geometry.alpha is not None, "Geometry must have alpha set before orienting panels."
        assert self.geomeetry.delta is not None, "Geometry must have flap deflection set before orienting panels."



    def solve_plate(self) -> None:
        # Placeholder for setting control points specific to Weisinger's Approximation

        c = self.geometry.chord['value']
        p = self.geometry.pivot['value']
        alpha = self.geometry.alpha['value']

        delta = self.geometry.flap_deflection_rad
        l_le = self.geometry.flap_length_le

        # Calculate main and flap chord lengths
        main_c = l_le * c
        flap_c = (1 - l_le) * c

        # Calulculate Vortex source locations.
        G1_x = 0.25 * main_c
        G2_x = main_c + 0.25 * flap_c

        # Calculate location of flow tangency conditions.
        V1_x = 0.75 * main_c
        V2_x = main_c + 0.75 * flap_c

        # Calculate distances from vortices to control points.
        r11 = G1_x - V1_x
        r12 = G1_x - V2_x
        r21 = G2_x - V1_x
        r22 = G2_x - V2_x

        # Calculate flow tangency conditions at control points.
        dzdx1 = math.tan(-1*alpha)  # Flat plate has zero slope on main section
        dzdx2 = math.tan(-1*(alpha + delta))  # Flap slope

        # Define influence coefficient matrix
        A = (1/(2.0*math.pi)) * np.array([[-1/r11, 1/r21],
                      [-1/r12, -1/r22]])
        
        # Define right-hand side vector using flow tangency conditions.
        b = np.array([
            (dzdx1 * math.cosh(alpha) - math.sin(alpha)) * self.U_inf,
            (dzdx2 * math.cosh(alpha + delta) - math.sin(alpha + delta)) * self.U_inf
        ])  

        G = np.linalg.solve(A, b)

        L = self.rho_inf * self.U_inf * (G[0] + G[1])
        cl = L / (0.5 * self.rho_inf * self.U_inf**2 * c)
        cd = 0.0  # Placeholder for drag coefficient calculation
        cm_cg = 0.0  # Placeholder for moment coefficient calculation
        
        results = {
            'L': L,
            'cl': cl,
            'cd': cd,
            'cm_cg': cm_cg
        }

        return results
    