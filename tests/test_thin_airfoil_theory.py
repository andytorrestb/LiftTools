import sys
from pathlib import Path

# Ensure 'src' is on path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from ThinAirfoilTheory import ThinAirfoilTheory  # noqa: E402


def test_thin_airfoil_solve_returns_expected_keys():
    tat = ThinAirfoilTheory(alpha_rad=0.123)
    res = tat.solve()
    assert isinstance(res, dict)
    assert 'cl' in res
    assert 'alpha_rad' in res
    assert abs(res['alpha_rad'] - 0.123) < 1e-12
