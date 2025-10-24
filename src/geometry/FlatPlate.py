from .Airfoil import Airfoil
import math
from typing import List


class FlatPlate(Airfoil):
    def __init__(self, chord: float, alpha_rad: float = 0.0):
        """Simple FlatPlate geometry container."""
        self.chord = float(chord)
        self.alpha_rad = float(alpha_rad)
        self.x: List[float] = []
        self.y: List[float] = []

    def set_alpha(self, alpha_rad: float):
        """Set angle of attack in radians."""
        self.alpha_rad = float(alpha_rad)

    def orient_to_alpha(self):
        """Create a flat plate geometry given chord and angle of attack.

        Args:
            chord (float): chord length
            alpha_rad (float): angle of attack in radians

        Returns:
            FlatPlate: instance with computed x, y coordinates
        """

        assert self.chord > 0.0, "Chord must be positive."
        assert self.alpha_rad >= -math.pi/2 and self.alpha_rad <= math.pi/2, "Alpha must be between -90 and 90 degrees in radians."


        # original coordinates: leading (0) and trailing (chord) edges
        x = [0.0, self.chord]
        y = [0.0, 0.0]

        # rotate about quarter-chord point c/4 by alpha_rad
        c4 = self.chord / 4.0
        ca = math.cos(-1*self.alpha_rad)
        sa = math.sin(-1*self.alpha_rad)

        x_rot: List[float] = []
        y_rot: List[float] = []
        for xi, yi in zip(x, y):
            xt = xi - c4
            yt = yi
            xr = xt * ca - yt * sa
            yr = xt * sa + yt * ca
            x_rot.append(xr + c4)
            y_rot.append(yr)

        self.x = x_rot
        self.y = y_rot

        return