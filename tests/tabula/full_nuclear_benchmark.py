"""Full two-center nuclear-attraction benchmark — Tabula vs VeloxChem.

For each (molecule, basis) it reports, best-of-N (N adaptive to size):
  - Tabula dense `compute` wall time (explicit exact, bypassing the auto-screen),
  - VeloxChem `NuclearPotentialDriver` wall time,
  - VLX/Tab — dense speedup over VeloxChem,
  - Tabula `compute_sparse` wall time, its stored footprint, and VLX/Tab(sp) —
    the sparse speedup over VeloxChem,
  - the maximum element deviation of Tabula's dense matrix from VeloxChem's
    (skipped — "-" — when a matrix would exceed the memory budget).

Charges are the molecule's own nuclei. Both engines use VeloxChem's positive
convention. `compute_sparse` is run at a lossless screening threshold. Run
inside the VeloxChem environment.
"""

import gc
import time

import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem import NuclearPotentialDriver
from veloxchem.tabulalib import TabulaNuclearAttractionDriver

GEOMETRY = {
    "c60":       "/Users/rinkevic/Development/VeloxChem/c60.xyz",
    "taxol":     "/Users/rinkevic/Development/VeloxChem/taxol.xyz",
    "olestra":   "/Users/rinkevic/Development/VeloxChem/olestra.xyz",
    "ubiquitin": "/Users/rinkevic/Downloads/ubiquitin.xyz",
}

PLAN = [
    ("c60",       ["def2-svp", "def2-tzvp", "def2-qzvp", "cc-pvdz"]),
    ("taxol",     ["def2-svp", "def2-tzvp", "def2-qzvp"]),
    ("olestra",   ["def2-svp", "def2-tzvp"]),
    ("ubiquitin", ["def2-svp", "def2-svpd", "cc-pvdz", "aug-cc-pvdz"]),
]

SCREEN_THRESHOLD = 1.0e-10
MAXDIFF_BUDGET = 1.5e9  # bytes per matrix — skip the dense max-diff above this


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


def vlx_compute(mol, bas):
    drv = NuclearPotentialDriver()
    try:
        return drv.compute(mol, bas)
    except Exception:
        return drv.compute(bas, mol)


def main():
    print("| molecule  | basis       | n_AO | dense/ms | VLX/ms | VLX/Tab | sparse/ms | VLX/Tab(sp) | stored | max-diff |")
    print("|-----------|-------------|-----:|--------:|-------:|--------:|----------:|------------:|-------:|---------:|")
    for name, bases in PLAN:
        mol = read_molecule(GEOMETRY[name])
        for basis_label in bases:
            bas = MolecularBasis.read(mol, basis_label, ostream=None)
            drv = TabulaNuclearAttractionDriver()

            # the sparse compute is cheap on large cells and gives the dimension
            sp = drv.compute_sparse(mol, bas, SCREEN_THRESHOLD)
            n = sp.dimension()
            stored = sp.stored_element_count() / (n * n)
            del sp; gc.collect()

            reps = reps_for(n)
            t_dense = measure(lambda: drv.compute(mol, bas, 0.0), reps)
            t_vlx = measure(lambda: vlx_compute(mol, bas), reps)
            t_sparse = measure(lambda: drv.compute_sparse(mol, bas, SCREEN_THRESHOLD), reps)

            if n * n * 8 <= MAXDIFF_BUDGET:
                ref = vlx_compute(mol, bas).full_matrix().to_numpy()
                dense = drv.compute(mol, bas, 0.0)
                diff = f"{float(np.max(np.abs(dense.to_numpy() - ref))):.0e}"
                del ref, dense; gc.collect()
            else:
                diff = "-"

            print(f"| {name:9s} | {basis_label:11s} | {n:5d} | {t_dense*1e3:6.1f} | {t_vlx*1e3:6.1f} | "
                  f"{t_vlx/t_dense:5.2f}x | {t_sparse*1e3:7.1f} | {t_vlx/t_sparse:8.2f}x | {stored*100:4.1f}% | {diff} |",
                  flush=True)


if __name__ == "__main__":
    main()
