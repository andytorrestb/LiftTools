import numpy as np

from .Airfoil import Airfoil

def parse_naca4_digits(digits: str) -> tuple[float, float, float]:
    """Parse NACA 4-digit code into parameters.

    Parameters
    - digits: string of 4 digits, e.g., "2412"

    Returns
    - (m, p, t): tuple where
      m = maximum camber as fraction of chord
      p = location of maximum camber as fraction of chord
      t = maximum thickness as fraction of chord

    Raises
    - ValueError if the input is not a valid 4-digit NACA code
    """
    if len(digits) != 4 or not digits.isdigit():
        raise ValueError("NACA code must be a string of 4 digits.")

    m = int(digits[0]) / 100.0
    p = int(digits[1]) / 10.0
    t = int(digits[2:4]) / 100.0

    return m, p, t

class NACA(Airfoil):
    def __init__(self, naca_code: str):
        # Parse NACA 4-digit code into parameters
        m, p, t = parse_naca4_digits(naca_code)
        self.naca_code = naca_code
        self.m = m
        self.p = p
        self.t = t
        
        # Initialize base Airfoil with default chord and alpha
        super().__init__(chord=1.0, alpha_rad=0.0, pivot_frac=0.25)


    def set_n_panels(self, n: int) -> None:
        """Set the number of panels for discretization."""
        assert n > 0, "Number of panels must be positive."
        assert n > 1, "Number of panels must be greater than 1."
        self.n_panels = int(n)

    def camber_line_naca4(self, m: float, p: float) -> None:
        """Set the camber line for a NACA 4-digit airfoil.

        Parameters
        - m: maximum camber as fraction of chord (e.g., 0.02 for 2%)
        - p: location of maximum camber as fraction of chord (e.g., 0.4 for 40%)
        """
        import sympy as sp

        x = sp.Symbol('x', real=True, positive=True)
        c = self.chord['symbol']

        # Define piecewise camber line
        z_expr = sp.Piecewise(
            ( (m / p**2) * (2 * p * (x / c) - (x / c)**2), (x / c) < p ),
            ( (m / (1 - p)**2) * ( (1 - 2 * p) + 2 * p * (x / c) - (x / c)**2 ), True )
        )

        # Set the camber line in the geometry
        # Note: call the correctly named setter
        self.set_camber(z_expr)
    
    def read_dat_file(self, filepath: str) -> None:
        """Read airfoil coordinates from a DAT file.

        Parameters
        - filepath: path to the DAT file containing airfoil coordinates
        """
        import numpy as np

        data = np.loadtxt(filepath, skiprows=1)
        x = data[:, 0].tolist()
        y = data[:, 1].tolist()

        # Split at the repeated leading-edge x≈0 so each series starts at x=0
        mid_index = max(1, len(x) // 2)
        self.x_upper = np.array([0.0] + x[:mid_index][::-1])
        self.y_upper = np.array([0.0] + y[:mid_index][::-1])
        self.x_lower = np.array(x[mid_index:])
        self.y_lower = np.array(y[mid_index:])

        # print("Airfoil coordinates loaded from DAT file.")
        # print(f"Upper surface points: {self.x_upper}")
        # print(f"Lower surface points: {self.x_lower}")

    pass

    def is_valid_dat(self) -> bool:
        """Validate loaded DAT coordinates.

        This verifies that the coordinate arrays exist, are finite, and that
        the x-coordinates for upper and lower surfaces are identical at each
        index (within a tiny numerical tolerance for floating-point data).

        Returns
        - True if validation passes, False otherwise.
        """
        import numpy as np

        # Ensure attributes exist and are non-empty
        required_attrs = ["x_upper", "y_upper", "x_lower", "y_lower"]
        for attr in required_attrs:
            if not hasattr(self, attr):
                return False
            val = getattr(self, attr)
            if val is None:
                return False

        x_upper = np.asarray(getattr(self, "x_upper"))
        x_lower = np.asarray(getattr(self, "x_lower"))

        # Basic shape and finiteness checks
        if x_upper.size == 0 or x_lower.size == 0:
            return False
        if x_upper.shape != x_lower.shape:
            return False
        if not (np.isfinite(x_upper).all() and np.isfinite(x_lower).all()):
            return False

        # Fast exact check first
        if np.array_equal(x_upper, x_lower):
            return True
        else:
            return False    

    def set_camber(self, z_expr = None) -> None:
        """Set the camber line expression for the airfoil.

        Parameters
        - z_expr: symbolic expression for camber line z(x)
        """
        if z_expr is None:
            assert self.x_upper is not None, "Airfoil coordinates (x_upper) must be set before setting camber."
            assert self.x_lower is not None, "Airfoil coordinates (x_lower) must be set before setting camber."
            assert self.y_upper is not None, "Airfoil coordinates (y_upper) must be set before setting camber."
            assert self.y_lower is not None, "Airfoil coordinates (y_lower) must be set before setting camber."
            assert self.is_valid_dat(), "Loaded DAT coordinates are not valid."

            self.z = 0.5*(self.y_upper + self.y_lower)  # Simplified assumption for camber line
        else:
            self.z = z_expr


    def plot_airfoil(self) -> None:
        """Plot the airfoil shape using stored coordinates."""
        import matplotlib.pyplot as plt

        plt.figure(figsize=(10, 4))
        plt.plot(self.x_upper, self.y_upper, label='Upper Surface')
        plt.plot(self.x_lower, self.y_lower, label='Lower Surface')
        plt.title(f'NACA {self.naca_code} Airfoil')
        plt.xlabel('x (chordwise)')
        plt.ylabel('y (thickness)')
        plt.axis('equal')
        plt.grid(True)

        if self.z is not None and isinstance(self.z, np.ndarray):
            
            plt.plot(self.x_upper, self.z, 'r--', label='Camber Line')

        plt.legend()
        plt.savefig(f"NACA_{self.naca_code}_airfoil.png")
