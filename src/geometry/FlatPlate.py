from .Airfoil import Airfoil
import math
from typing import List


class FlatPlate(Airfoil):
    def __init__(self, chord: float, thickness: float = 0.0, alpha_rad: float = 0.0):
        """Simple FlatPlate geometry container."""
        self.chord = float(chord)
        self.thickness = float(thickness)
        self.alpha_rad = float(alpha_rad)
        self.x: List[float] = []
        self.y: List[float] = []

    @classmethod
    def from_dimensions(cls, chord: float, alpha_rad: float = 0.0) -> "FlatPlate":
        """Create a flat plate geometry given chord and angle of attack.

        Args:
            chord (float): chord length
            alpha_rad (float): angle of attack in radians

        Returns:
            FlatPlate: instance with computed x, y coordinates
        """
        inst = cls(chord=chord, thickness=0.0, alpha_rad=alpha_rad)

        # original coordinates: leading (0) and trailing (chord) edges
        x = [0.0, inst.chord]
        y = [0.0, 0.0]

        # rotate about quarter-chord point c/4 by alpha_rad
        c4 = inst.chord / 4.0
        ca = math.cos(-1*inst.alpha_rad)
        sa = math.sin(inst.alpha_rad)

        x_rot: List[float] = []
        y_rot: List[float] = []
        for xi, yi in zip(x, y):
            xt = xi - c4
            yt = yi
            xr = xt * ca - yt * sa
            yr = xt * sa + yt * ca
            x_rot.append(xr + c4)
            y_rot.append(yr)

        inst.x = x_rot
        inst.y = y_rot

        return inst