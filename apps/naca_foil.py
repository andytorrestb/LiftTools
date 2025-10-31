import apps_header

import math

import numpy as np
import matplotlib.pyplot as plt

from geometry import NACA
from models import WeisingersApprox as wa

if __name__ == "__main__":
    # Step 0) Define study constants
    U_inf = 90.0  # m/s
    rho_inf = 1.2  # kg/m^3
    k = 0.80 # chordwise location of flap hinge
    c = 1.0  # chord length
    alpha = [math.radians(a_deg) for a_deg in np.linspace(-5, 15, 5)]  # angles of attack in radians
    delta = [math.radians(d_deg) for d_deg in [15, 30, 45, 60]]  # flap deflection angles in radians
    
    # Study 1) NACA 0012 airfoil (main + flap)
    # load from dat --> collapse to camber profile --> rotate for alpha and flap deflection
    # Run Weisinger's Approximation model.
    # Read in airfoil data. Plot surfaces and camber line.
    naca0012 = NACA.NACA(naca_code="0012")
    naca0012.read_dat_file("../data/NACA0012.dat")
    naca0012.set_camber()
    naca0012.plot_airfoil()

    Cl = []
    Cl_TAT = []
    Cm_cg = []
    Cm_cg_TAT = []
    for a in alpha:
        # Initialize Weisinger's Approximation model.
        wa_naca0012 = wa.WeisingersApprox(geometry=naca0012)

        # Set flap and attack angles. Rotate camberline accordingly.
        wa_naca0012.geometry.set_alpha(a)  # e.g., 5 degrees
        wa_naca0012.geometry.set_flap_properties(
            k=k,
            delta=delta[0]  # e.g., 10 degrees
        )

        # Discretize and orient panels.
        wa_naca0012.set_panels(n_panels=10)
        wa_naca0012.plot_panels()
        # Rotate main wing by alpha about pivot and flap by (alpha + delta) about hinge
        wa_naca0012.orient_panels()
        wa_naca0012.plot_panels(file_suffix="_oriented")

        # Set flow conditions and solve.
        wa_naca0012.set_flow_conditions(U_inf=U_inf, rho_inf=rho_inf)
        wa_naca0012.set_points()
        wa_naca0012.plot_points()
        wa_naca0012.compute_distances()
        wa_naca0012.compute_panel_lengths()
        wa_naca0012.compute_panel_normals()
        results = wa_naca0012.solve()
        Cl.append(results["cl"])
        Cm_cg.append(results["cm_cg"])

        Cl_TAT.append(2 * math.pi * a)
        Cm_cg_TAT.append(0.0)


    # Plot results.
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot([math.degrees(a) for a in alpha], Cl, marker='o', label="Weisinger's Approximation")
    plt.plot([math.degrees(a) for a in alpha], Cl_TAT, marker='x', label="Thin Airfoil Theory")
    plt.title("Lift Coefficient vs Angle of Attack")
    plt.xlabel("Angle of Attack (degrees)")
    plt.ylabel("Lift Coefficient (Cl)")
    plt.grid()
    plt.tight_layout()
    plt.legend()
    plt.subplot(1, 2, 2)
    plt.plot([math.degrees(a) for a in alpha], Cm_cg, marker='o', label="Weisinger's Approximation")
    plt.plot([math.degrees(a) for a in alpha], Cm_cg_TAT, marker='x', label="Thin Airfoil Theory")
    plt.title("Moment Coefficient about CG vs Angle of Attack")
    plt.xlabel("Angle of Attack (degrees)")
    plt.ylabel("Moment Coefficient (Cm_cg)")
    plt.grid()
    plt.tight_layout()
    plt.legend()
    plt.savefig(f"weisingers_approx_results_naca{naca0012.naca_code}.png")
    plt.clf()

    # Study 2) NACA 0012 airfoil (high resolution case N=100)
    wa_naca0012_highres = wa.WeisingersApprox(geometry=naca0012)
    wa_naca0012_highres.geometry.set_alpha(alpha[2])  # e.g., 5 degrees
    wa_naca0012_highres.geometry.set_flap_properties(
        k=k,
        delta=delta[0]  # e.g., 10 degrees
    )
    wa_naca0012_highres.set_flow_conditions(U_inf=U_inf, rho_inf=rho_inf)
    wa_naca0012_highres.set_panels(n_panels=100)
    wa_naca0012_highres.orient_panels()
    wa_naca0012_highres.set_points()
    wa_naca0012_highres.compute_distances()
    wa_naca0012_highres.compute_panel_lengths()
    wa_naca0012_highres.compute_panel_normals()
    results_highres = wa_naca0012_highres.solve()

    print(f"High resolution case (N=100) for NACA {naca0012.naca_code}:")
    print(f"Lift Coefficient (Cl): {results_highres['cl']:.4f}")
    print(f"Moment Coefficient about CG (Cm_cg): {results_highres['cm_cg']:.4f}")

    G = results_highres['G']
    # print(f"Circulation (G): {G:.4f} m^2/s")
    n_panels = np.linspace(0, 1, wa_naca0012_highres.n_panels)
    # Plot high resolution results.
    plt.figure(figsize=(10, 5))
    plt.plot(n_panels, G, marker='o', label="Circulation (G)")
    plt.xlabel("Panel Index")
    plt.ylabel("Circulation (m^2/s)")
    plt.title(f"Circulation Distribution for NACA {naca0012.naca_code} (N=100)")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"weisingers_approx_circulation_naca{naca0012.naca_code}.png")
    plt.clf()

    # Study 3) NACA 2412 airfoil (low resolution)
    naca2412 = NACA.NACA(naca_code="2412")
    naca2412.read_dat_file("../data/NACA2412.dat")
    naca2412.set_camber()
    naca2412.plot_airfoil()

    wa_naca2412 = wa.WeisingersApprox(geometry=naca2412)
    wa_naca2412.geometry.set_flap_properties(
        k=k,
        delta=delta[2]  # e.g., 10 degrees
    )
    wa_naca2412.geometry.set_alpha(alpha[2])  # e.g., 5 degrees 
    wa_naca2412.set_flow_conditions(U_inf=U_inf, rho_inf=rho_inf)
    wa_naca2412.set_panels(n_panels=10)
    wa_naca2412.plot_panels()







    pass