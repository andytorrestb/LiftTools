import math

import numpy as np

from models.AirFoilModel import AirfoilModel
from models.WeisingersApprox import WeisingersApprox

class ThinAirfoilTheory(AirfoilModel):
    def set_naca_camber(self) -> None:
        assert self.geometry is not None, "Geometry must be set before setting camber."
        assert self.geometry.z is not None, "Geometry must have a camber line defined."
        self.c = self.geometry.chord['value']
        self.x =np.linspace(0, self.c, num=100)
        self.z = np.interp(self.x, self.geometry.x_lower, self.geometry.z)

    def alpha_n(self, n: int) -> float:
        assert self.x is not None, "x-coordinates must be set before computing alpha_n."
        assert self.c is not None, "Chord length must be set before computing alpha_n."
        
        theta =  (self.c/2.0) * (1 - np.cos(np.pi * self.x / self.c))
        dzdx_theta = np.gradient(self.z, theta)
        integrand = dzdx_theta * np.cos(n*theta)
        return (2.0 / np.pi) * np.trapz(integrand, theta)

    def solve(self):

        a = self.alpha_rad

        results = {
            "cl": 2*math.pi*(a - self.alpha_n(0)),
            "cm": (math.pi/4) * (self.alpha_n(2) - self.alpha_n(1)),
            "alpha_rad": self.alpha_rad,
        }
        return results
