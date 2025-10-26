from __future__ import annotations

from typing import List
from unittest import case

from symbolic.utils import sym


class Airfoil:
    """Base class for airfoil geometries.

    Concrete geometries (e.g., FlatPlate, NACA) may inherit from this class.
    Provides common storage for chord, pivot, and angle-of-attack along with
    coordinate buffers x and y.
    """

    def __init__(self, chord: float, alpha_rad: float = 0.0, pivot_frac: float = 0.25) -> None:
        """Initialize base geometry values.

        Parameters
        - chord: chord length (must be > 0)
        - alpha_rad: initial angle of attack in radians
        - pivot_frac: pivot location as a fraction of chord (0..1), default 0.25
        """
        c_val = float(chord)
        p_abs = float(pivot_frac) * c_val

        self.chord = {
            'value': c_val,
            'symbol': sym.Symbol('c', real=True),
        }

        # Store absolute x-location of the pivot in the same units as chord
        self.pivot = {
            'value': p_abs,
            'symbol': sym.Symbol('p', real=True),
        }

        self.alpha = {
            'value': float(alpha_rad),
            'symbol': sym.Symbol('a', real=True),
        }

        self.x: List[float] = []
        self.y: List[float] = []

    def set_alpha(self, alpha_rad: float) -> None:
        """Set the numeric angle of attack in radians."""
        self.alpha['value'] = float(alpha_rad)

    def orient_to_alpha(self) -> None:
        """Update geometry coordinates for the current alpha if supported.

        Subclasses that provide `compute_endpoint_values()` will have this
        method call it to refresh x and y based on the current parameters.
        """
        if hasattr(self, 'compute_endpoint_values') and callable(getattr(self, 'compute_endpoint_values')):
            # type: ignore[attr-defined]
            self.compute_endpoint_values()  # noqa: F401