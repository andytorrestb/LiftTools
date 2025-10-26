from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Callable, Iterable, List, Optional, Tuple
import math
import numpy as np

class AirfoilModel(ABC):
    """Abstract base for airfoil models.

    Subclasses must implement "solve" which performs the computation and
    returns a result dictionary. Shared helpers for geometry sampling and
    angle conversion are provided here.
    """

    def __init__(
        self,
        geometry: Geometry | None = None,
        alpha_rad: float = 0.0,
        panels: int = 100,
    ) -> None:
        self.geometry = geometry or Geometry()
        self.alpha_rad = float(alpha_rad)
        self.panels = int(panels)

    # --- small contract -------------------------------------------------
    # solve() -> dict with keys at least: 'cl' (float) and 'alpha_rad' (float)

    @abstractmethod
    def solve(self) -> dict:
        """Run the model and return results.

        Implementations should return a dictionary with a minimum of:
        - 'cl': lift coefficient
        - 'alpha_rad': angle of attack in radians
        Additional keys (pressure distribution, gamma, etc.) are allowed.
        """

    def set_cambmer(self, z_expr) -> None:

        assert self.geometry is not None, "Geometry must be set before setting camber."
        assert hasattr(self.geometry.camber[''], 'set_camber'), (
            "Geometry must implement set_camber method."
        )
        assert callable(self.geometry.set_camber), (
            "Geometry's set_camber must be callable."
        )
        assert self.geom
        assert z_expr is not None, "Camber expression must be provided."
        self.geometry.set_camber(z_expr)
        return

    def set_flow_conditions(
        self,
        U_inf: float = 1.0, # m/s
        rho_inf: float = 1.225, # kg/m^3
    ) -> None:
        """Set flow conditions for the model.

        Parameters
        - alpha_rad: angle of attack in radians (optional)
        """
        self.U_inf = float(U_inf)
        self.rho_inf = float(rho_inf)
        return
