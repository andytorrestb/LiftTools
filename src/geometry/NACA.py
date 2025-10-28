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
        self.set_cambmer(z_expr)
    
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
        self.x_upper = [0.0] + x[:mid_index][::-1]
        self.y_upper = [0.0] + y[:mid_index][::-1]
        self.x_lower = x[mid_index:]
        self.y_lower = y[mid_index:]

        print("Airfoil coordinates loaded from DAT file.")
        print(f"Upper surface points: {self.x_upper}")
        print(f"Lower surface points: {self.x_lower}")

    pass


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
        plt.legend()
        plt.savefig(f"NACA_{self.naca_code}_airfoil.png")