import apps_header

import math

import numpy as np

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

    # Initialize Weisinger's Approximation model.
    wa_naca0012 = wa.WeisingersApprox(geometry=naca0012)

    # Set flap and attack angles. Rotate camberline accordingly.
    wa_naca0012.geometry.set_alpha(alpha[2])  # e.g., 5 degrees
    wa_naca0012.geometry.set_flap_properties(
        k=k,
        delta=delta[2]  # e.g., 10 degrees
    )

    # Discretize and orient panels.
    wa_naca0012.set_panels(n_panels=10)
    wa_naca0012.plot_panels()
    # Rotate main wing by alpha about pivot and flap by (alpha + delta) about hinge
    wa_naca0012.orient_panels()
    # Do NOT re-discretize after orientation; it would overwrite the rotated coordinates
    wa_naca0012.plot_panels(file_suffix="_oriented")

    # Set flow conditions and solve.
    wa_naca0012.set_flow_conditions(U_inf=U_inf, rho_inf=rho_inf)
    wa_naca0012.set_points()
    wa_naca0012.plot_points()
    wa_naca0012.compute_distances()
    results = wa_naca0012.solve()


    print(wa_naca0012.geometry.k)
    print(wa_naca0012.geometry.delta)

    # Study 2) NACA 0012 airfoil (high resolution case N=100)
    
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