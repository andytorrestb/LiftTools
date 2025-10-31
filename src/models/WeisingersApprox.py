import math
import numpy as np
import matplotlib.pyplot as plt

from models.AirFoilModel import AirfoilModel

def rotate_point(point: tuple[float, float], pivot: tuple[float, float], theta: float) -> tuple[float, float]:
    """Rotate a point (x, z) about a pivot by angle theta (radians, CCW-positive)."""
    # Flip sign of theta for CW-positive rotation
    theta = -theta
    
    # Translate point to pivot-relative coordinates
    x_rel = point[0] - pivot[0]
    z_rel = point[1] - pivot[1]

    # Apply standard 2D rotation without reusing the updated x
    x_new = x_rel * math.cos(theta) - z_rel * math.sin(theta)
    z_new = x_rel * math.sin(theta) + z_rel * math.cos(theta)

    # Translate back to world coordinates
    return (x_new + pivot[0], z_new + pivot[1])

class WeisingersApprox(AirfoilModel):
    def set_panels(self, n_panels: int) -> None:
        assert self.geometry.z is not None, "Camber line must be set before discretizing."
        # Target total number of panels (segments)
        self.n_panels = int(n_panels)
        # establish local name references for important data.
        z = self.geometry.z
        c = self.geometry.chord['value']
        k = self.geometry.k

        # Determine number of panels on main wing and flap
        # Use a clean partition so that:
        #   total_panels = n_wing_panels + n_flap_panels == self.n_panels
        n_wing_panels = int(round(k * self.n_panels))
        n_flap_panels = self.n_panels - n_wing_panels
        self.n_wing_panels = n_wing_panels
        self.n_flap_panels = n_flap_panels

        # Discretize and interpolate panels on main wing and flap
        # We create N_panels + 1 points to define N_panels segments.
        wing_panel_x = np.linspace(0, k*c, n_wing_panels + 1)
        wing_panel_z = np.interp(wing_panel_x, self.geometry.x_lower, z)

        flap_panel_x = np.linspace(k*c, c, n_flap_panels + 1)
        flap_panel_z = np.interp(flap_panel_x, self.geometry.x_lower, z)

        # Set to object properties for later use.
        self.wing_panel_x = wing_panel_x
        self.wing_panel_z = wing_panel_z    
        self.flap_panel_x = flap_panel_x
        self.flap_panel_z = flap_panel_z
        # Debug prints to verify discretization
        print("Set panels:")
        print("  Wing panel x:", self.wing_panel_x)
        print("  Wing panel z:", self.wing_panel_z)
        print("  Flap panel x:", self.flap_panel_x)
        print("  Flap panel z:", self.flap_panel_z)

    def plot_panels(self, file_suffix="") -> None:
        plt.figure(figsize=(10, 4))
        # plt.plot(self.geometry.x_lower, self.geometry.z, 'b-', label='Camber Line')
        plt.plot(self.wing_panel_x, self.wing_panel_z, 'ro--', label='Wing Panels')
        plt.plot(self.flap_panel_x, self.flap_panel_z, 'go--', label='Flap Panels')
        plt.xlabel('x (chordwise)')
        plt.ylabel('z (camber)')
        plt.ylim(-0.25, 0.25)
        plt.title('Airfoil Camber Line with Discretized Panels')
        plt.legend()
        plt.grid(True)
        plt.savefig(f"weisingers_approx_panels{file_suffix}_naca{self.geometry.naca_code}.png")         

    def set_points(self) -> None:
        # Placeholder for setting control points specific to Weisinger's Approximation
        # Set arrays for various points needed in WA calculations.
        N = self.n_panels
        C = np.zeros((N, 2))  # Control points (midpoint of panel endpoints)
        QC = np.zeros((N, 2))  # Vortex source points (quarter-chord locations)
        TC = np.zeros((N, 2))  # Tangency condition points (three-quarter-chord locations)

        # Add points for wing panels (indices reference segment i between i and i+1)
        for i in range(self.n_wing_panels):
            C_i = np.array([(self.wing_panel_x[i] + self.wing_panel_x[i+1]) / 2.0,
                    (self.wing_panel_z[i] + self.wing_panel_z[i+1]) / 2.0 ])
            C[i, :] = C_i

            Q_i = np.array([ self.wing_panel_x[i] + 0.25 * (self.wing_panel_x[i+1] - self.wing_panel_x[i]),
                             self.wing_panel_z[i] + 0.25 * (self.wing_panel_z[i+1] - self.wing_panel_z[i]) ])
            QC[i, :] = Q_i

            T_i = np.array([ self.wing_panel_x[i] + 0.75 * (self.wing_panel_x[i+1] - self.wing_panel_x[i]),
                             self.wing_panel_z[i] + 0.75 * (self.wing_panel_z[i+1] - self.wing_panel_z[i]) ])
            TC[i, :] = T_i

        # Add points for flap panels
        k = self.geometry.k
        c = self.geometry.chord['value']
        for i in range(self.n_wing_panels, self.n_panels):
            j = i - self.n_wing_panels  # index within flap segments
            C_i = np.array([(self.flap_panel_x[j] + self.flap_panel_x[j+1]) / 2.0,
                            (self.flap_panel_z[j] + self.flap_panel_z[j+1]) / 2.0 ])
            C[i, :] = C_i
            
            Q_i = np.array([ self.flap_panel_x[j] + 0.25 * (self.flap_panel_x[j+1] - self.flap_panel_x[j]),
                             self.flap_panel_z[j] + 0.25 * (self.flap_panel_z[j+1] - self.flap_panel_z[j]) ])
            QC[i, :] = Q_i

            T_i = np.array([ self.flap_panel_x[j] + 0.75 * (self.flap_panel_x[j+1] - self.flap_panel_x[j]),
                             self.flap_panel_z[j] + 0.75 * (self.flap_panel_z[j+1] - self.flap_panel_z[j]) ])
            TC[i, :] = T_i
        
        self.C = C
        self.QC = QC
        self.TC = TC
        print("Set points:")
        print("  Control points C:", self.C)
        print(self.C.shape)
        print("  Vortex points QC:", self.QC)
        print("  Tangency points TC:", self.TC)

    def compute_distances(self) -> None:
        # Calculate distances Rij from each vortex point to each control point
        R = np.zeros((len(self.QC), len(self.TC)))
        for i in range(len(self.QC)):
            for j in range(len(self.TC)):
                dx = self.QC[i][0] - self.TC[j][0]
                dz = self.QC[i][1] - self.TC[j][1]
                R_ij = math.sqrt(dx**2 + dz**2)
                R[i, j] = R_ij

        self.R = R

    def compute_panel_lengths(self) -> None:
        # There are self.n_panels segments defined by self.total_x/z (which have N+1 points)
        S = np.zeros(self.n_panels)
        for i in range(self.n_panels):
            S[i] = math.sqrt((self.total_x[i+1] - self.total_x[i])**2 + (self.total_z[i+1] - self.total_z[i])**2)
        self.S = S

    def compute_panel_normals(self) -> None:
        N = np.zeros((self.n_panels, 2))
        for i in range(self.n_panels):
            dx = self.total_x[i+1] - self.total_x[i]
            dz = self.total_z[i+1] - self.total_z[i]

            N[i, :] = np.array([-dz, dx]) / math.sqrt(dx**2 + dz**2)
        self.N = N



    def plot_points(self) -> None:
        # Placeholder for plotting control points, vortex points, etc.
        plt.figure(figsize=(10, 4))
        plt.plot(self.wing_panel_x, self.wing_panel_z, 'ro--', label='Wing Panels')
        plt.plot(self.flap_panel_x, self.flap_panel_z, 'go--', label='Flap Panels')
        plt.plot(self.C[:, 0], self.C[:, 1], 'kx', label='Control Points C')
        plt.plot(self.QC[:, 0], self.QC[:, 1], 'm+', label='Vortex Points QC')
        plt.plot(self.TC[:, 0], self.TC[:, 1], 'cs', label='Tangency Points TC')
        plt.xlabel('x (chordwise)')
        plt.ylabel('z (camber)')
        plt.ylim(-0.25, 0.25)
        plt.title('Airfoil Camber Line with Panels and Control Points')
        plt.legend()
        plt.grid(True)
        plt.savefig(f"weisingers_approx_points_naca{self.geometry.naca_code}.png")

        pass

    def solve(self) -> dict:
        # Placeholder for Weisinger's Approximation solution method
        print("Solving using Weisinger's Approximation...")

        # Geometry has already been rotated to alpha/delta in orient_panels().
        # Therefore, keep the freestream aligned with the x-axis here.
        alpha = float(self.geometry.alpha['value'])
        delta = float(self.geometry.delta)

        A = np.zeros((self.n_panels, self.n_panels))
        b = np.zeros(self.n_panels)

        # Build influence matrix A so that rows correspond to control/tangency
        # points (TC) and columns correspond to vortex/source points (QC).
        # self.R was built as R[q, t] = distance from QC[q] to TC[t].
        for i in range(self.n_panels):  # row: TC index
            for j in range(self.n_panels):  # col: QC index
                # WA sign convention: a vortex induces opposite-signed normal
                # velocity downstream vs upstream. Using panel index order as a
                # proxy for streamwise direction along the camber line.
                if i < j:  # TC upstream of QC
                    ijk = 1.0
                else:      # TC downstream of QC
                    ijk = -1.0

                # Use distance from QC[j] to TC[i]
                A[i, j] = ijk / (2.0 * math.pi * self.R[j, i])

        # Set flow velocity for wing and flap panels.
        # Use a single, global freestream vector. Geometry has already been
        # oriented by alpha/delta, so do NOT rotate the freestream again.
        U_inf = np.zeros((self.n_panels, 2))
        for i in range(self.n_panels):
            U_inf[i] = (self.U_inf, 0.0)
        
        
        # Set up right-hand side vector b based on flow tangency conditions.
        for i in range(self.n_panels):
            b[i] = -np.dot(U_inf[i, :], self.N[i, :])

        # Solve for vortex strengths G
        G = np.linalg.solve(A, b)
        print("Solved vortex strengths G:", G)
        L = self.rho_inf * self.U_inf * np.sum(G)
        M_cg = 0

        for i in range(self.n_panels):
            # pivot['value'] is already an absolute x-location (in chord units)
            r_x = self.C[i, 0] - (self.geometry.pivot['value'])
            M_cg += self.rho_inf * self.U_inf * G[i] * r_x
        
        c = self.geometry.chord['value']
        cl = L / (0.5 * self.rho_inf * self.U_inf**2 * c)
        cd = 0.0  # Cd is zero for TAT
        cm_cg = M_cg / (0.5 * self.rho_inf * self.U_inf**2 * c**2)
        results = {
            'L': L,
            'cl': cl,
            'cd': cd,
            'cm_cg': cm_cg
        }
        print("Weisinger's Approximation results:", results)
        return results

    def set_flap(self, deflection_rad: float, length_le: float) -> None:
        # Placeholder for setting flap deflection specific to Weisinger's Approximation
        self.geometry.flap_deflection_rad = float(deflection_rad)
        self.geometry.flap_length_le = float(length_le)
        pass


    def orient_panels(self) -> None:
        """Update panel coordinates for the current alpha if supported.

        Subclasses that provide `compute_endpoint_values()` will have this
        method call it to refresh x and y based on the current parameters.
        """
        assert self.geometry.alpha is not None, "Geometry must have alpha set before orienting panels."
        assert self.geometry.delta is not None, "Geometry must have flap deflection set before orienting panels."
        alpha = float(self.geometry.alpha['value'])
        delta = float(self.geometry.delta)
        p = float(self.geometry.pivot['value'])  # pivot for main wing rotation (e.g., quarter-chord)
        k = float(self.geometry.k)               # hinge fraction along chord
        c = float(self.geometry.chord['value'])

        # Rotate main-wing panels about pivot p by alpha
        wing_rot = [
            rotate_point((x, z), (p, 0.0), alpha)
            for x, z in zip(self.wing_panel_x, self.wing_panel_z)
        ]
        self.wing_panel_x = np.array([pt[0] for pt in wing_rot], dtype=float)
        self.wing_panel_z = np.array([pt[1] for pt in wing_rot], dtype=float)

        # Rotate flap panels around pivot point by alpha.
        flap_rot = [
            rotate_point((x, z), (p, 0.0), alpha)
            for x, z in zip(self.flap_panel_x, self.flap_panel_z)
        ]

        new_hinge = flap_rot[0]  # New hinge location after main wing rotation
        flap_rot = [
            rotate_point(point, new_hinge, delta)
            for point in flap_rot
        ]

        self.flap_panel_x = np.array([pt[0] for pt in flap_rot], dtype=float)
        self.flap_panel_z = np.array([pt[1] for pt in flap_rot], dtype=float)

        wing_panel_z = self.wing_panel_z[:-1]
        wing_panel_x = self.wing_panel_x[:-1]

        total_z = np.array(wing_panel_z.tolist() + self.flap_panel_z.tolist())
        total_x = np.array(wing_panel_x.tolist() + self.flap_panel_x.tolist())

        self.total_x = total_x
        self.total_z = total_z

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
    