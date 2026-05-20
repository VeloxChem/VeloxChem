"""Full two-center charge-dipole benchmark — Tabula vs VeloxChem.

For each (molecule, basis) it reports, best-of-N (N adaptive to size):
  - Tabula dense `compute` wall time (the contracted matrix Sum_N d_N·field_N,
    explicit exact, bypassing screening),
  - VeloxChem `NuclearPotentialGeom010Driver` wall time (it returns the three
    field-component matrices),
  - VLX/Tab — dense speedup over VeloxChem,
  - Tabula `compute_sparse` wall time, its stored footprint, and VLX/Tab(sp).

The point dipoles sit at the molecule's own nuclei with unit moments — a
QM/MM-like load comparable to the nuclear benchmark's charges-at-nuclei. Note
the engines do different amounts of output work: Tabula returns the single
contracted charge-dipole matrix (what a QM/MM Fock build needs), VeloxChem
returns the three field components. Correctness is validated separately in
test_tabula_chargedipole.py (the VeloxChem multi-site CMatrices Python
extraction is unreliable, so it is not used for a max-diff here). Run inside the
VeloxChem environment.
"""

import gc
import time

import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem import NuclearPotentialGeom010Driver
from veloxchem.tabulalib import TabulaChargeDipoleDriver

GEOMETRY = {
    "c60":       "/Users/rinkevic/Development/VeloxChem/c60.xyz",
    "taxol":     "/Users/rinkevic/Development/VeloxChem/taxol.xyz",
    "olestra":   "/Users/rinkevic/Development/VeloxChem/olestra.xyz",
    "ubiquitin": "/Users/rinkevic/Downloads/ubiquitin.xyz",
}

PLAN = [
    ("c60",       ["def2-svp", "def2-tzvp", "cc-pvdz"]),
    ("taxol",     ["def2-svp", "def2-tzvp"]),
    ("olestra",   ["def2-svp", "def2-tzvp"]),
    ("ubiquitin", ["def2-svp", "cc-pvdz"]),
]

SCREEN_THRESHOLD = 1.0e-10


def read_molecule(path):
    with open(path) as handle:
        lines = handle.read().splitlines()
    return Molecule.read_str("\n".join(lines[2:]), "angstrom")


def reps_for(n):
    return 3 if n < 2000 else (2 if n < 6000 else 1)


def measure(make, reps):
    make(); gc.collect()  # warm-up
    best = float("inf")
    for _ in range(reps):
        gc.collect()
        start = time.perf_counter()
        result = make()
        best = min(best, time.perf_counter() - start)
        del result; gc.collect()
    return best


def main():
    print("| molecule  | basis       | n_AO | dense/ms | VLX/ms | VLX/Tab | sparse/ms | VLX/Tab(sp) | stored |")
    print("|-----------|-------------|-----:|--------:|-------:|--------:|----------:|------------:|-------:|")
    for name, bases in PLAN:
        mol = read_molecule(GEOMETRY[name])
        coords = [list(c) for c in np.array(mol.get_coordinates_in_bohr())]
        moments = [[1.0, 1.0, 1.0] for _ in coords]
        weights = [1.0 for _ in coords]
        for basis_label in bases:
            bas = MolecularBasis.read(mol, basis_label, ostream=None)
            drv = TabulaChargeDipoleDriver()
            g010 = NuclearPotentialGeom010Driver()

            sp = drv.compute_sparse(mol, bas, moments, coords, SCREEN_THRESHOLD)
            n = sp.dimension()
            stored = sp.stored_element_count() / (n * n)
            del sp; gc.collect()

            reps = reps_for(n)
            t_dense = measure(lambda: drv.compute(mol, bas, moments, coords, 0.0), reps)
            t_vlx = measure(lambda: g010.compute(mol, bas, weights, coords), reps)
            t_sparse = measure(lambda: drv.compute_sparse(mol, bas, moments, coords, SCREEN_THRESHOLD), reps)

            print(f"| {name:9s} | {basis_label:11s} | {n:5d} | {t_dense*1e3:6.1f} | {t_vlx*1e3:6.1f} | "
                  f"{t_vlx/t_dense:5.2f}x | {t_sparse*1e3:7.1f} | {t_vlx/t_sparse:8.2f}x | {stored*100:4.1f}% |",
                  flush=True)


if __name__ == "__main__":
    main()
