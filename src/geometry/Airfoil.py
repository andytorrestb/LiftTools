import numpy as np


class Airfoil:
    @classmethod
    def from_flat_plate(cls, chord: float):
        """Create a flat plate airfoil geometry.

        Args:
            chord: Chord length of the flat plate.
            thickness: Thickness of the flat plate (usually small).

        Returns:
            An Airfoil instance representing the flat plate.
        """
        # Simple flat plate coordinates
        x = np.array([0.0, chord])
        y = np.array([0.0, 0.0])
        return cls(x, y)

