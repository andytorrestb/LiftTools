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
    
    pass