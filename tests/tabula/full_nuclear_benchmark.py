"""Full two-center nuclear-attraction benchmark — Tabula vs VeloxChem.

For each (molecule, basis) it reports, best-of-N:
  - Tabula dense `compute` wall time,
  - VeloxChem `NuclearPotentialDriver` wall time,
  - Tabula `compute_sparse` wall time and its stored footprint (vs dense),
  - the maximum element deviation of Tabula's dense matrix from VeloxChem's.

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
    "c60":     "/Users/rinkevic/Development/VeloxChem/c60.xyz",
    "taxol":   "/Users/rinkevic/Development/VeloxChem/taxol.xyz",
    "olestra": "/Users/rinkevic/Development/VeloxChem/olestra.xyz",
}

PLAN = [
    ("c60",     ["def2-svp", "def2-tzvp", "def2-qzvp", "cc-pvdz"]),
    ("taxol",   ["def2-svp", "def2-tzvp", "def2-qzvp"]),
    ("olestra", ["def2-svp", "def2-tzvp"]),
]

REPS = 3
SCREEN_THRESHOLD = 1.0e-10


def read_molecule(path):
    with open(path) as handle:
        lines = handle.read().splitlines()
    return Molecule.read_str("\n".join(lines[2:]), "angstrom")


def best_time(make):
    discard = make(); del discard; gc.collect()
    best = float("inf")
    for _ in range(REPS):
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
    print("| molecule | basis     | n_AO | dense/ms | VLX/ms | VLX/Tab | sparse/ms | stored | max-diff |")
    print("|----------|-----------|-----:|--------:|-------:|--------:|----------:|-------:|---------:|")
    for name, bases in PLAN:
        mol = read_molecule(GEOMETRY[name])
        for basis_label in bases:
            bas = MolecularBasis.read(mol, basis_label, ostream=None)
            drv = TabulaNuclearAttractionDriver()

            dense = drv.compute(mol, bas, 0.0)  # explicit exact (bypass the auto-screen default)
            n = dense.rows()
            t_dense = best_time(lambda: drv.compute(mol, bas, 0.0))
            t_vlx = best_time(lambda: vlx_compute(mol, bas))
            t_sparse = best_time(lambda: drv.compute_sparse(mol, bas, SCREEN_THRESHOLD))

            sp = drv.compute_sparse(mol, bas, SCREEN_THRESHOLD)
            stored = sp.stored_element_count() / (n * n)

            ref = vlx_compute(mol, bas).full_matrix().to_numpy()
            diff = float(np.max(np.abs(dense.to_numpy() - ref)))

            print(f"| {name:8s} | {basis_label:9s} | {n:5d} | {t_dense*1e3:6.1f} | {t_vlx*1e3:6.1f} | "
                  f"{t_vlx/t_dense:5.2f}x | {t_sparse*1e3:7.1f} | {stored*100:4.1f}% | {diff:.0e} |")


if __name__ == "__main__":
    main()
