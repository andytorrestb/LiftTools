# LiftTools

LiftTools is a collection of utilities and experimental solvers for analyzing two-dimensional airfoils. The project couples symbolic geometry descriptions with classic aerodynamic panel methods to explore lift, moment, and circulation characteristics for canonical airfoil families.

## Repository Structure

The top-level layout below highlights the major packages and scripts in the repository:

```
LiftTools/
├── apps/
│   ├── apps_header.py
│   ├── flat_plate.py
│   ├── naca_foil.py
│   └── tat_wa_compare.py
├── data/
│   ├── NACA0012.dat
│   └── NACA2412.dat
├── src/
│   ├── geometry/
│   │   ├── Airfoil.py
│   │   ├── FlatPlate.py
│   │   └── NACA.py
│   ├── models/
│   │   ├── AirFoilModel.py
│   │   ├── ThinAirfoilTheory.py
│   │   ├── ThinAirfoilTheoryOld.py
│   │   └── WeisingersApprox.py
│   └── symbolic/
│       └── utils.py
└── tests/
    ├── test_flatplate.py
    └── test_thin_airfoil_theory.py
```

## Key Components

- **`src/geometry/`** – Defines airfoil primitives. The base `Airfoil` class stores geometry parameters (chord, pivot, and angle of attack) as dictionaries containing both numeric values and SymPy symbols so they can feed symbolic derivations and numeric evaluations. `FlatPlate` and `NACA` extend the base with analytic camber constructions and DAT ingestion utilities.
- **`src/symbolic/`** – Houses helpers that normalize SymPy substitution and lambdification, letting geometry classes keep a single representation that can be sampled numerically on demand.
- **`src/models/`** – Implements aerodynamic solvers. `AirFoilModel` defines the solver protocol, while `WeisingersApprox` builds a discrete vortex panel system and the `ThinAirfoilTheory` variants evaluate analytic approximations of lift and moment coefficients.
- **`apps/`** – Collection of runnable scripts that combine geometry definitions and solvers for experiments (e.g., comparing thin-airfoil predictions with panel results or sweeping angle of attack).
- **`data/`** – Reference airfoil coordinate files for common NACA sections used by the applications and tests.
- **`tests/`** – Pytest-based regression tests that validate geometry setup and thin-airfoil predictions.

## Getting Started

1. Create a Python virtual environment and install project dependencies (see `requirements.txt` if present or inspect module imports for guidance).
2. Run the automated checks:
   ```bash
   pytest
   ```
3. Execute scripts under `apps/` to generate plots or perform comparative studies. For example:
   ```bash
   python apps/naca_foil.py
   ```

## Contributing

Please open an issue or submit a pull request with a clear description of the changes. Include unit tests or sample scripts that demonstrate new behavior whenever possible.
