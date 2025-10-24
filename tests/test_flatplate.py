import math
import sys
from pathlib import Path

# Ensure 'src' is on path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from geometry.FlatPlate import FlatPlate  # noqa: E402


def test_flatplate_endpoints_rotation():
    chord = 2.0
    p = chord / 4.0
    for alpha in [0.0, math.radians(10.0), math.radians(-10.0), math.radians(30.0)]:
        plate = FlatPlate(chord=chord, alpha_rad=alpha)
        plate.orient_to_alpha()
        assert len(plate.x) == 2 and len(plate.y) == 2

        x0_exp = p + (0.0 - p) * math.cos(alpha)
        y0_exp = (0.0 - p) * math.sin(alpha)
        x1_exp = p + (chord - p) * math.cos(alpha)
        y1_exp = (chord - p) * math.sin(alpha)

        assert math.isclose(plate.x[0], x0_exp, rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(plate.y[0], y0_exp, rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(plate.x[1], x1_exp, rel_tol=0.0, abs_tol=1e-12)
        assert math.isclose(plate.y[1], y1_exp, rel_tol=0.0, abs_tol=1e-12)


def test_flatplate_sampling_includes_endpoints():
    chord = 2.0
    alpha = math.radians(15.0)
    plate = FlatPlate(chord=chord, alpha_rad=alpha)
    plate.orient_to_alpha()

    xs, ys = plate.sample(5)
    assert len(xs) == 5 and len(ys) == 5
    assert math.isclose(xs[0], plate.x[0], rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(ys[0], plate.y[0], rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(xs[-1], plate.x[1], rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(ys[-1], plate.y[1], rel_tol=0.0, abs_tol=1e-12)
