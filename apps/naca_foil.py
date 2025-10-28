import apps_header

import math

import numpy as np

from geometry import NACA

if __name__ == "__main__":
    # Step 0) Define study constants
    U_inf = 90.0  # m/s
    rho_inf = 1.2  # kg/m^3
    k = 0.75 # chordwise location of flap hinge
    c = 1.0  # chord length
    alpha = [math.radians(a_deg) for a_deg in np.linspace(-5, 15, 5)]  # angles of attack in radians
    delta = [math.radians(d_deg) for d_deg in [0, 5, 10, 15]]  # flap deflection angles in radians
    
    # Study 1) NACA 0012 airfoil (main + flap)
    # load from dat --> collapse to camber profile --> rotate for alpha and flap deflection
    # Run Weisinger's Approximation model.
    naca0012 = NACA.NACA(naca_code="0012")
    naca0012.read_dat_file("../data/NACA0012.dat")
    naca0012.plot_airfoil()

    # Study 2) NACA 0012 airfoil (high resolution case N=100)
    
    # Study 3) NACA 2412 airfoil (low resolution)
    naca2412 = NACA.NACA(naca_code="2412")
    naca2412.read_dat_file("../data/NACA2412.dat")
    naca2412.plot_airfoil()






    pass