import apps_header

from matplotlib import pyplot as plt
import numpy as np
import math

import models.WeisingersApprox as wa
import models.ThinAirfoilTheory as tat

import geometry.FlatPlate as fp

if __name__ == "__main__":
    # Load flat plate and plot at various angles of attack
    # Provides a simple sanity check of geometry rotation.
    chord = 1.0  # chord length
    alpha = np.linspace(0, math.pi/2.0, 5)  # angles of attack in radians
    flat_plate_geom = fp.FlatPlate(chord=chord)

    plt.figure(figsize=(8, 4))
    for a in alpha:
         # Load flate plate geometry with coordinates rotated
         # according to angle of attack (alpha_rad).
        #  flat_plate_geom.set_alpha(alpha_rad)
         flat_plate_geom.alpha['value'] = a
         flat_plate_geom.compute_endpoint_values()
         
         # Plit coordinates for comparison.
         plt.plot(flat_plate_geom.x, flat_plate_geom.y, label=f"Flat Plate at {a} deg")
         plt.title(f"Flat Plate Airfoil Geometry at {a} deg")
         plt.xlabel("x")
         plt.ylabel("y")
         plt.legend()
         plt.grid(True)
    plt.savefig(f"flat_plate_rotations.png")

    # # Run Thin Airfoil Theory model on flat plate using same geometry and alpha values.
    # for a in alpha:
    #     alpha_rad = a * 3.14159 / 180.0  # convert to radians
    #     flat_plate_geom = fp.FlatPlate.from_dimensions(chord=chord, alpha_rad=alpha_rad)

    #     tat_model = tat.ThinAirfoilTheory(
    #         geometry=flat_plate_geom,
    #         alpha_rad=alpha_rad,
    #         panels=1,
    #     )
    #     results = tat_model.solve()
    #     cl = results['cl']
    #     print(f"TAT: Alpha = {a} deg, CL = {cl:.4f}")

