import apps_header

import math

import numpy as np

if __name__ == "__main__":
    # Step 0) Define study constants
    naca0012 = "0012"
    naca2412 = "2412"
    U_inf = 90.0  # m/s
    rho_inf = 1.2  # kg/m^3
    k = 0.75 # chordwise location of flap hinge
    c = 1.0  # chord length
    alpha = [math.radians(a_deg) for a_deg in np.linspace(-5, 15, 5)]  # angles of attack in radians
    delta = [math.radians(d_deg) for d_deg in [0, 5, 10, 15]]  # flap deflection angles in radians
    
    # Study 1) NACA 0012 airfoil (main + flap)
    
    # Study 2) NACA 2412 airfoil (high resolution case N=100)
    
    # Study 3) NACA 2412 airfoil (low resolution)






    pass