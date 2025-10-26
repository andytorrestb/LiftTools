import apps_header

from matplotlib import pyplot as plt
import numpy as np
import math

import models.WeisingersApprox as wa
import models.ThinAirfoilTheory as tat

import geometry.FlatPlate as fp

if __name__ == "__main__":
     # Step 1) Load flat plate and plot at various angles of attack
     # Provides a simple sanity check of geometry rotation.
     chord = 1.0  # chord length
     alpha = np.linspace(0, math.pi / 12.0, 3)  # angles of attack in radians
     flat_plate_geom = fp.FlatPlate(chord=chord)

     plt.figure(figsize=(8, 4))
     for a in alpha:
          # Load flat plate geometry with coordinates rotated
          # according to angle of attack (alpha_rad).
          # flat_plate_geom.set_alpha(alpha_rad)
          flat_plate_geom.alpha['value'] = a
          flat_plate_geom.compute_endpoint_values()

          # Plot coordinates for comparison.
          plt.plot(
               flat_plate_geom.x,
               flat_plate_geom.y,
               label=f"Flat Plate at {math.degrees(a):.0f} deg",
          )
          plt.title(f"Flat Plate Airfoil Geometry at {math.degrees(a):.0f} deg")
          plt.xlabel("x")
          plt.ylabel("y")
          plt.legend()
          plt.grid(True)
     plt.savefig("flat_plate_rotations.png")
     plt.clf()

     # Step 2) Run Thin Airfoil Theory model on flat plate using same geometry and alpha values.
     for a in alpha:
          flat_plate_geom.alpha['value'] = a
          flat_plate_geom.compute_endpoint_values()

          tat_model = tat.ThinAirfoilTheory(
               geometry=flat_plate_geom,
               alpha_rad=a,
               panels=1,
          )

          tat_model.set_z_expr(foil_type="plate")
          results = tat_model.solve()
          cl = results['cl']
          print(f"TAT: Alpha = {math.degrees(a):.0f} deg, CL = {cl:.4f}")

     # Step 3) Plot profiles of flapped plate for visual verification.
     deflection_rad = math.radians(-15)  # 15 degree flap deflection
     flap_l_le = 0.7                    # hinge at 70% chord from LE
     plt.figure(figsize=(8, 4))
     for a in alpha:
          # Update geometry state for this angle
          flat_plate_geom.alpha['value'] = a
          flat_plate_geom.compute_endpoint_values()

          # Configure flap parameters on geometry
          flat_plate_geom.flap_deflection_rad = deflection_rad
          flat_plate_geom.flap_length_le = flap_l_le

          # Plot flapped plate
          wing, flap =flat_plate_geom.plot_flapped_plate(subplot=True)

          plt.plot(wing[0], wing[1], label=f"Wing at {math.degrees(a):.0f} deg")
          plt.plot(flap[0], flap[1], label=f"Flap at {math.degrees(a):.0f} deg")
     plt.title(f"Flapped Flat Plate at {math.degrees(a):.0f} deg AoA")
     plt.xlabel("x")
     plt.ylabel("y")
     plt.legend()
     plt.grid(True)
     plt.savefig("flapped_flat_plate.png")
     plt.clf()
     # Step 4) Run Weisinger's Approximation model on flat plate with a flap.
     for a in alpha:
          flat_plate_geom.alpha['value'] = a
          flat_plate_geom.compute_endpoint_values()

          wa_model = wa.WeisingersApprox(
               geometry=flat_plate_geom,
               alpha_rad=a,
               panels=2,
          )

          wa_model.set_flap(deflection_rad=deflection_rad, length_le=0.7)
          wa_model.set_flow_conditions(
               U_inf=10.0,  # m/s
               rho_inf=1.225,  # kg/m^3
          )
          results = wa_model.solve_plate()
          cl = results['cl']
          print(f"Weisinger's Approx: Alpha = {math.degrees(a):.0f} deg, CL = {cl:.4f}")
