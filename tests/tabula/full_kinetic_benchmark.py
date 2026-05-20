"""Full two-center kinetic benchmark — Tabula vs VeloxChem.

For each (molecule, basis) it reports, best-of-5:
  - Tabula dense `compute` wall time,
  - VeloxChem `KineticEnergyDriver` wall time,
  - Tabula `compute_sparse` wall time and its stored footprint (vs dense),
  - the maximum element deviation of Tabula's dense kinetic from VeloxChem's
    (skipped — shown "-" — when the two full matrices would not fit in RAM).

`compute_sparse` is run at a lossless screening threshold. Run inside the
VeloxChem environment.
"""

import gc
import time

import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import KineticEnergyDriver
from veloxchem.tabulalib import TabulaKineticDriver

GEOMETRY = {
    "c60":       "/Users/rinkevic/Development/VeloxChem/c60.xyz",
    "c240":      "/Users/rinkevic/Development/VeloxChem/c240.xyz",
    "taxol":     "/Users/rinkevic/Development/VeloxChem/taxol.xyz",
    "olestra":   "/Users/rinkevic/Development/VeloxChem/olestra.xyz",
    "ubiquitin": "/Users/rinkevic/Downloads/ubiquitin.xyz",
}

FULL_BASES = ["def2-SVP", "def2-SVPD", "def2-TZVP", "def2-TZVPD", "def2-QZVP", "def2-QZVPD",
              "cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ"]

PLAN = [
    ("c60",       FULL_BASES),
    ("c240",      FULL_BASES),
    ("taxol",     FULL_BASES),
    ("olestra",   FULL_BASES),
    ("ubiquitin", ["def2-SVP", "def2-SVPD", "cc-pVDZ", "aug-cc-pVDZ"]),
]

REPETITIONS = 5
SCREEN_THRESHOLD = 1.0e-12
VALIDATE_BYTE_BUDGET = 24.0e9   # skip the two-matrix diff above this


def read_molecule(path):
    with open(path) as handle:
        lines = handle.read().splitlines()
    return Molecule.read_str("\n".join(lines[2:]), "angstrom")


def best_time(make):
    """Best-of-REPETITIONS wall time; the result is freed each iteration."""
    discard = make()
    del discard
    gc.collect()
    best = float("inf")
    for _ in range(REPETITIONS):
        gc.collect()
        start = time.perf_counter()
        result = make()
        best = min(best, time.perf_counter() - start)
        del result
        gc.collect()
    return best


def max_abs_difference(a, b):
    step = max(1, (4096 * 4096) // max(1, a.shape[1]))
    worst = 0.0
    for start in range(0, a.shape[0], step):
        end = min(start + step, a.shape[0])
        worst = max(worst, float(np.max(np.abs(a[start:end] - b[start:end]))))
    return worst


print("| molecule  | basis       | n_AO | dense/ms | VLX/ms | VLX/Tab | sparse/ms | VLX/Tab(sp) | stored | max-diff |", flush=True)
print("|-----------|-------------|-----:|--------:|-------:|--------:|----------:|------------:|-------:|---------:|", flush=True)

for molecule_name, bases in PLAN:
    molecule = read_molecule(GEOMETRY[molecule_name])
    for basis_label in bases:
        try:
            basis = MolecularBasis.read(molecule, basis_label, ostream=None)

            tabula_time = best_time(lambda: TabulaKineticDriver().compute(molecule, basis))
            veloxchem_time = best_time(lambda: KineticEnergyDriver().compute(molecule, basis))
            sparse_time = best_time(lambda: TabulaKineticDriver().compute_sparse(molecule, basis, SCREEN_THRESHOLD))

            sparse = TabulaKineticDriver().compute_sparse(molecule, basis, SCREEN_THRESHOLD)
            dimension = sparse.dimension()
            stored = sparse.stored_element_count() / (dimension * dimension)
            del sparse
            gc.collect()

            ratio = veloxchem_time / tabula_time

            if 2.0 * dimension * dimension * 8.0 < VALIDATE_BYTE_BUDGET:
                tabula_matrix = TabulaKineticDriver().compute(molecule, basis)
                tabula_values = tabula_matrix.to_numpy()
                veloxchem_full = KineticEnergyDriver().compute(molecule, basis).full_matrix().to_numpy()
                deviation = max_abs_difference(tabula_values, veloxchem_full)
                del tabula_matrix, tabula_values, veloxchem_full
                gc.collect()
                deviation_text = f"{deviation:.0e}"
            else:
                deviation_text = "-"

            print(f"| {molecule_name:9s} | {basis_label:11s} | {dimension:5d} | "
                  f"{tabula_time*1e3:6.1f} | {veloxchem_time*1e3:6.1f} | {ratio:5.2f}x | "
                  f"{sparse_time*1e3:7.1f} | {veloxchem_time/sparse_time:8.2f}x | {stored*100:4.1f}% | {deviation_text} |",
                  flush=True)
        except Exception as error:
            print(f"{molecule_name:10s} {basis_label:12s}  FAILED: {error}", flush=True)
            gc.collect()
