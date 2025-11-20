import apps_header

import math

import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
import aerosandbox.tools.pretty_plots as p

from geometry import NACA
from models import WeisingersApprox as wa

if __name__ == "__main__":
    # Step 0) Define study constants
    U_inf = 90.0  # m/s
    rho_inf = 1.2  # kg/m^3
    k = 0.80 # chordwise location of flap hinge
    c = 1.0  # chord length
    alpha = [math.radians(a_deg) for a_deg in np.linspace(-5, 15, 5)]  # angles of attack in radians
    delta = [math.radians(d_deg) for d_deg in [15, 30, 45, 60]]  # flap deflection angles
    
    # Study 1) NACA 0012 airfoil (main + flap)
    # load from dat --> collapse to camber profile --> rotate for alpha and flap deflection
    # Run Weisinger's Approximation model.
    # Read in airfoil data. Plot surfaces and camber line.
    naca = NACA.NACA(naca_code="0012")
    naca.read_dat_file("../data/NACA0012.dat")
    naca.set_camber()
    naca.plot_airfoil()

    # Initialize Weisinger's Approximation model.
    wa_naca = wa.WeisingersApprox(geometry=naca)

    # Initialize neuralfoil analysis.
    af = asb.Airfoil(name="NACA0012")

    # Initialize result storage
    Cl_nf = []
    Cl_WA = []
    Cm_nf = []
    Cm_WA = []
    for a in alpha:

        # Set flap and attack angles. Rotate camberline accordingly.
        wa_naca.geometry.set_alpha(a)  # e.g., 5 degrees
        wa_naca.geometry.set_flap_properties(
            k=k,
            delta=delta[0]  # e.g., 10 degrees
        )

        # Discretize and orient panels.
        wa_naca.set_panels(n_panels=100)
        wa_naca.orient_panels()
        wa_naca.set_flow_conditions(U_inf=U_inf, rho_inf=rho_inf)
        wa_naca.set_points()
        wa_naca.compute_distances()
        wa_naca.compute_panel_lengths()
        wa_naca.compute_panel_normals()
        wa_naca.solve()
        results_wa = wa_naca.solve()
        Cl_WA.append(results_wa["cl"])
        Cm_WA.append(results_wa["cm_cg"])

        # Compute Reynolds number and Mach number for neuralfoil
        Re = (rho_inf * U_inf * c) / (1.81e-5)  # Dynamic viscosity of air at ~25C
        mach = U_inf / 343.0  # Speed of sound at ~20C

        print(f'Reynolds number: {Re:.2e}, Mach number: {mach:.4f}')

        # Neuralfoil analysis
        # NeuralFoil expects alpha in DEGREES; convert from radians used elsewhere
        results_af = af.get_aero_from_neuralfoil(
            alpha=math.degrees(a),
            Re=Re,
            mach=mach
        )

        # Ensure scalar values from potential array outputs
        cl_nf_val = float(np.array(results_af["CL"]).reshape(-1)[0])
        cm_nf_val = float(np.array(results_af["CM"]).reshape(-1)[0])

        Cl_nf.append(cl_nf_val)
        Cm_nf.append(cm_nf_val)
        print(
            f"Angle of Attack (degrees): {math.degrees(a):.2f}, "
            f"Lift Coefficient (Cl) from Neuralfoil: {cl_nf_val:.4f}, "
            f"Lift Coefficient (Cl) from Weisinger's Approximation: {results_wa['cl']:.4f}"
        )

    # Plot results.
    # Prepare angle array in degrees for plotting
    alpha_deg = [math.degrees(a) for a in alpha]

    # Create a single figure with two subplots and save once
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    # Left subplot: Lift Coefficient
    axs[0].plot(alpha_deg, Cl_nf, label="Neuralfoil Cl")
    axs[0].plot(alpha_deg, Cl_WA, label="Weisinger's Cl")
    axs[0].set_xlabel("Angle of Attack (degrees)")
    axs[0].set_ylabel("Lift Coefficient (Cl)")
    axs[0].set_title("Lift Coefficient Comparison")
    axs[0].grid(True)
    axs[0].legend()

    # Right subplot: Moment Coefficient
    axs[1].plot(alpha_deg, Cm_nf, label="Neuralfoil Cm")
    axs[1].plot(alpha_deg, Cm_WA, label="Weisinger's Cm")
    axs[1].set_xlabel("Angle of Attack (degrees)")
    axs[1].set_ylabel("Moment Coefficient (Cm)")
    axs[1].set_title("Moment Coefficient Comparison")
    axs[1].grid(True)
    axs[1].legend()

    fig.suptitle("Neuralfoil vs Weisinger's Approximation", fontsize=12)
    fig.savefig("nf_comparison.png", dpi=300)
    plt.close(fig)
