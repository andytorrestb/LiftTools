import apps_header

import math

import numpy as np
import matplotlib.pyplot as plt

from geometry import NACA
from models import ThinAirfoilTheory as TAT, WeisingersApprox

if __name__ == "__main__":
    # Step 0) Define study constants
    U_inf = 90.0  # m/s
    rho_inf = 1.2  # kg/m^3
    alpha = [math.radians(a_deg) for a_deg in np.linspace(-5, 15, 5)]  # angles of attack in radians
    alpha_deg = [math.degrees(a) for a in alpha]
    
    # Study 1) NACA 0012 airfoil (main + flap)
    naca0012 = NACA.NACA(naca_code="2412")
    naca0012.read_dat_file("../data/NACA2412.dat")
    naca0012.set_camber()
    naca0012.plot_airfoil()

    Cl_TAT = []
    Cm_cg_TAT = []
    CL_WA = []
    Cm_cg_WA = []
    for a in alpha:
        tat_naca0012 = TAT.ThinAirfoilTheory(geometry=naca0012)
        # tat_naca0012.set_naca_camber(m=naca0012.m, p=naca0012.p)
        tat_naca0012.set_naca_camber()
        tat_naca0012.set_flow_conditions(U_inf=U_inf, rho_inf=rho_inf)
        tat_naca0012.alpha_rad = a
        results = tat_naca0012.solve()
        Cl_TAT.append(results["cl"])
        Cm_cg_TAT.append(results['cm'])
        print(
            f"Angle of Attack (degrees): {math.degrees(a):.2f}, "
            f"Lift Coefficient (Cl) from Thin Airfoil Theory: {results['cl']:.4f}"
        )

        # Import above brings in the module; instantiate the class from it
        wa = WeisingersApprox.WeisingersApprox(geometry=naca0012)
        wa.set_flow_conditions(U_inf=U_inf, rho_inf=rho_inf)
        wa.geometry.set_alpha(a)
        wa.geometry.set_flap_properties(
            k=0.80,
            delta=math.radians(15)  # example flap deflection
        )
        wa.set_panels(n_panels=10)
        wa.orient_panels()
        wa.set_points()
        wa.compute_distances()
        wa.compute_panel_lengths()
        wa.compute_panel_normals()
        wa_results = wa.solve()
        print(
            f"Angle of Attack (degrees): {math.degrees(a):.2f}, "
            f"Lift Coefficient (Cl) from Weisinger's Approximation: {wa_results['cl']:.4f}, "
            f"Moment Coefficient about CG (Cm_cg): {wa_results['cm_cg']:.4f}"
        )
        CL_WA.append(wa_results["cl"])
        Cm_cg_WA.append(wa_results["cm_cg"])

    # Merge the two plots into a horizontal subplot layout
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True)

    # Left subplot: Lift Coefficient
    axes[0].plot(alpha_deg, Cl_TAT, label="Thin Airfoil Theory", marker="o")
    axes[0].plot(alpha_deg, CL_WA, label="Weisinger's Approximation", marker="x")
    axes[0].set_xlabel("Angle of Attack (degrees)")
    axes[0].set_ylabel("Lift Coefficient (Cl)")
    axes[0].set_title("Lift Coefficient vs Angle of Attack")
    axes[0].grid()
    axes[0].legend()

    # Right subplot: Moment Coefficient about CG
    axes[1].plot(alpha_deg, Cm_cg_TAT, label="Thin Airfoil Theory", marker="o")
    axes[1].plot(alpha_deg, Cm_cg_WA, label="Weisinger's Approximation", marker="x")
    axes[1].set_xlabel("Angle of Attack (degrees)")
    axes[1].set_ylabel("Moment Coefficient about CG (Cm_cg)")
    axes[1].set_title("Moment Coefficient about CG vs Angle of Attack")
    axes[1].grid()
    axes[1].legend()

    fig.tight_layout()
    fig.savefig("compare_tat_vs_wa_subplots.png")