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
        alpha_deg: float = 0.0,
        panels: int = 100,
    ) -> None:
        self.geometry = geometry or Geometry()
        self.alpha_deg = float(alpha_deg)
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

    # --- helpers --------------------------------------------------------

    @property
    def alpha_rad(self) -> float:
        return math.radians(self.alpha_deg)

    def sample_camber(self, n: int | None = None) -> Tuple[np.ndarray, np.ndarray]:
        """Return arrays x (0..1) and z(x) sampled from camber_func or camber_points.

        If no camber is provided, returns zeros.
        """
        n = n or max(3, self.panels // 4)
        if self.geometry.camber_func is not None:
            xs = np.linspace(0.0, 1.0, n)
            zs = np.array([self.geometry.camber_func(x) for x in xs], dtype=float)
            return xs, zs

        if self.geometry.camber_points is not None:
            pts = sorted(self.geometry.camber_points, key=lambda p: p[0])
            xs, zs = zip(*pts)
            return np.array(xs, dtype=float), np.array(zs, dtype=float)

        # default: zero camber
        xs = np.linspace(0.0, 1.0, n)
        zs = np.zeros_like(xs)
        return xs, zs

    def camber_slope_mean(self, n: int | None = None) -> float:
        """Approximate the mean camber line slope d(z)/d(x) over 0..1.

        Used by simple analytical models (thin-airfoil effect of camber).
        """
        xs, zs = self.sample_camber(n)
        # central differences
        dx = np.diff(xs)
        dz = np.diff(zs)
        slopes = dz / dx
        # weight by dx
        mean = float(np.sum(slopes * dx) / np.sum(dx)) if slopes.size else 0.0
        return mean

    def chordwise_panels(self, clustering: str = "cosine") -> np.ndarray:
        """Return panel x-locations along chord [0,1].

        clustering: 'uniform' or 'cosine' (cosine clustering places more nodes near
        the leading/trailing edges). Returns `panels+1` node coordinates.
        """
        m = self.panels
        if clustering == "uniform":
            return np.linspace(0.0, 1.0, m + 1)
        # cosine clustering (standard for panel methods)
        beta = np.linspace(0.0, math.pi, m + 1)
        return 0.5 * (1.0 - np.cos(beta))

