"""The range separated RI-JK gradient, against four centres and against finite differences.

    python rs_scf_gradient_check.py [--cases water,hydroxyl] [--functionals ...]

Two comparisons are made of every case, and they answer different questions.

**Against the four-centre gradient.** Both sides are analytic gradients of their own
converged density, and the densities are not the same one: the resolution of the
identity fits the two-electron integrals and the four-centre way does not. So the
agreement here is bounded by the fitting error and lands around 1e-4, which is what
the plain (non range separated) rows show as well. It catches a gradient that is
wrong in shape or grossly wrong in a coefficient, and it is **not** sharp enough to
catch a coefficient wrong by a few percent.

**Against a finite difference of the SCF energy of the same method.** This is the
decisive check. The energy and the gradient come from one code path and one fitted
density, so nothing of the fitting cancels or hides: if the attenuated exchange
carried the wrong factor, the analytic gradient would not be the derivative of the
energy it was built from, whatever it agreed with elsewhere. A coefficient wrong by
a factor shows here as a difference of the size of the term, not of the fitting.

B3LYP is carried as a control. It is not range separated and goes through the plain
entry, so its rows say what the two comparisons cost when nothing new is exercised
-- the four-centre column reads the fitting error and the finite difference column
reads the noise floor of the grid and the convergence.

The four-centre finite difference is run for one range separated functional as well,
to show that the floor of that column is a property of the differencing and not of
the driver being checked.
"""
import argparse
import sys

import numpy as np
import veloxchem as vlx
from veloxchem.outputstream import OutputStream
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver

BASIS = "def2-svp"
AUX_BASIS = "def2-universal-jkfit"
# NOTE: 1e-9 and not tighter. At 1e-10 a plain B3LYP water does not converge in two
# hundred iterations -- the same gradient noise floor the c60 case ran into -- and a
# finite difference of a run which did not converge is not a finite difference. The
# floor this leaves is the convergence over the step, around 1e-7, which is well
# under what the grid contributes.
CONV_THRESH = 1.0e-9
GRID_LEVEL = 6
STEP = 0.005  # bohr

# NOTE: the geometries are deliberately off equilibrium. At a minimum the whole
# gradient is small and the terms of it cancel, so a coefficient wrong by a few
# percent would show as a small absolute difference; a stretched bond leaves a
# gradient of a few hundredths for the check to be read against.
CASES = {
    "water": {
        "symbols": ["O", "H", "H"],
        "coordinates": [[0.0, 0.0, 0.1173], [0.0, 0.8500, -0.5200],
                        [0.0, -0.7000, -0.4300]],
        "charge": 0,
        "multiplicity": 1,
    },
    "hydroxyl": {
        "symbols": ["O", "H"],
        "coordinates": [[0.0, 0.0, 0.0], [0.0, 0.0, 1.1000]],
        "charge": 0,
        "multiplicity": 2,
    },
}

FUNCTIONALS = ["B3LYP", "CAM-B3LYP", "WB97X-D4"]

# The functionals whose four-centre finite difference is taken as well.
FOUR_CENTRE_FD = ["CAM-B3LYP"]

BOHR_PER_ANGSTROM = 1.0 / 0.52917721092


def make_molecule(case, coordinates_bohr):
    molecule = vlx.Molecule(case["symbols"], coordinates_bohr, units="bohr")
    molecule.set_charge(case["charge"])
    molecule.set_multiplicity(case["multiplicity"])
    return molecule


def make_driver(case, functional, ri):
    scf_class = (ScfRestrictedDriver
                 if case["multiplicity"] == 1 else ScfUnrestrictedDriver)
    driver = scf_class(ostream=OutputStream(None))
    driver.conv_thresh = CONV_THRESH
    driver.max_iter = 200
    driver.grid_level = GRID_LEVEL
    driver.xcfun = functional
    if ri:
        driver.ri_jk = True
        driver.ri_jk_simd = True
        driver.ri_auxiliary_basis = AUX_BASIS.upper()
        driver.ri_mode = "in_memory"
    return driver


def energy(case, functional, ri, coordinates_bohr):
    molecule = make_molecule(case, coordinates_bohr)
    basis = vlx.MolecularBasis.read(molecule, BASIS.upper(), ostream=None)
    driver = make_driver(case, functional, ri)
    results = driver.compute(molecule, basis)
    if not driver.scf_results:
        raise SystemExit(f"{functional}: the SCF did not converge")
    return results["scf_energy"]


def gradient(case, functional, ri, coordinates_bohr):
    molecule = make_molecule(case, coordinates_bohr)
    basis = vlx.MolecularBasis.read(molecule, BASIS.upper(), ostream=None)
    driver = make_driver(case, functional, ri)
    driver.compute(molecule, basis)
    if not driver.scf_results:
        raise SystemExit(f"{functional}: the SCF did not converge")
    grad_driver = ScfGradientDriver(driver)
    grad_driver.ostream = OutputStream(None)
    grad_driver.compute(molecule, basis)
    return grad_driver.gradient.copy()


def finite_difference(case, functional, ri, coordinates_bohr):
    """The central difference of the SCF energy, one pair of runs a coordinate."""
    grad = np.zeros_like(coordinates_bohr)
    for iatom in range(coordinates_bohr.shape[0]):
        for icoord in range(3):
            plus = coordinates_bohr.copy()
            plus[iatom, icoord] += STEP
            minus = coordinates_bohr.copy()
            minus[iatom, icoord] -= STEP
            grad[iatom, icoord] = (energy(case, functional, ri, plus) -
                                   energy(case, functional, ri, minus)) / (
                                       2.0 * STEP)
    return grad


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", default=None)
    parser.add_argument("--functionals", default=None)
    args = parser.parse_args()

    names = ([c.strip() for c in args.cases.split(",")]
             if args.cases else list(CASES))
    functionals = ([f.strip() for f in args.functionals.split(",")]
                   if args.functionals else FUNCTIONALS)

    rows = []
    for name in names:
        case = CASES[name]
        coordinates = np.array(case["coordinates"]) * BOHR_PER_ANGSTROM
        shell = "closed" if case["multiplicity"] == 1 else "open"

        for functional in functionals:
            g_ri = gradient(case, functional, True, coordinates)
            g_4c = gradient(case, functional, False, coordinates)
            fd_ri = finite_difference(case, functional, True, coordinates)

            row = {
                "case": name,
                "shell": shell,
                "functional": functional,
                "ri_vs_4c": np.max(np.abs(g_ri - g_4c)),
                "ri_vs_fd": np.max(np.abs(g_ri - fd_ri)),
                "4c_vs_fd": None,
                "largest": np.max(np.abs(g_ri)),
            }
            if functional in FOUR_CENTRE_FD:
                fd_4c = finite_difference(case, functional, False, coordinates)
                row["4c_vs_fd"] = np.max(np.abs(g_4c - fd_4c))
            rows.append(row)

            print(f"  {name:9s} {functional:10s} "
                  f"ri-4c {row['ri_vs_4c']:.2e}  ri-fd {row['ri_vs_fd']:.2e}  "
                  f"4c-fd " +
                  ("--------" if row["4c_vs_fd"] is None else
                   f"{row['4c_vs_fd']:.2e}"),
                  flush=True)

    print()
    print(f"{'case':10s} {'shell':7s} {'functional':11s} {'|g|max':>10s} "
          f"{'ri vs 4c':>10s} {'ri vs fd':>10s} {'4c vs fd':>10s}")
    for row in rows:
        four = ("" if row["4c_vs_fd"] is None else f"{row['4c_vs_fd']:.2e}")
        print(f"{row['case']:10s} {row['shell']:7s} {row['functional']:11s} "
              f"{row['largest']:10.2e} {row['ri_vs_4c']:10.2e} "
              f"{row['ri_vs_fd']:10.2e} {four:>10s}")

    # NOTE: the finite difference is what settles it, so that is what decides the
    # exit status. The four-centre column is reported and not judged: its size is
    # the fitting error, which is a property of the approximation and not of this
    # code being right or wrong.
    floor = max((row["4c_vs_fd"] for row in rows if row["4c_vs_fd"]),
                default=1.0e-5)
    tolerance = max(10.0 * floor, 1.0e-5)
    bad = [row for row in rows if row["ri_vs_fd"] > tolerance]
    print()
    if bad:
        print(f"FAILED at a tolerance of {tolerance:.1e}: " +
              ", ".join(f"{r['case']}/{r['functional']}" for r in bad))
        return 1
    print(f"passed: every analytic gradient is the derivative of its own energy "
          f"to better than {tolerance:.1e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
