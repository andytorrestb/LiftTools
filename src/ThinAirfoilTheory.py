from AirFoilModel import AirfoilModel

class ThinAirfoilTheory(AirfoilModel):

    def solve(self) -> dict:
        cl = 0.0
        results = {'cl': cl}
        return results
    pass