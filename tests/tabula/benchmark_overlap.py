"""Overlap benchmark — tabula::OverlapDriver vs VeloxChem's OverlapDriver.

For each (molecule, basis) it best-of-5 times each engine sequentially
(never concurrent) and validates by the maximum absolute element difference.

Run inside the VeloxChem environment:
    python tests/tabula/benchmark_overlap.py
"""

import gc
import time

import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.tabulalib import TabulaOverlapDriver

GEOMETRY_DIR = "/Users/rinkevic/Development/VeloxChem"
MOLECULES = ["c60", "taxol", "c240"]
BASES = ["def2-SVP", "def2-SVPD", "def2-TZVP", "def2-TZVPD", "def2-QZVP", "def2-QZVPD"]
REPETITIONS = 5


def read_molecule(name):
    with open(f"{GEOMETRY_DIR}/{name}.xyz") as handle:
        lines = handle.read().splitlines()
    # standard xyz — line 0 is the atom count, line 1 a comment
    return Molecule.read_str("\n".join(lines[2:]), "angstrom")


def best_time(make):
    """Best-of-REPETITIONS wall time of `make`, keeping the last result."""
    discard = make()                      # warm-up
    del discard
    gc.collect()

    best = float("inf")
    result = None
    for _ in range(REPETITIONS):
        del result
        gc.collect()
        start = time.perf_counter()
        result = make()
        best = min(best, time.perf_counter() - start)
    return best, result


def max_abs_difference(a, b):
    """Maximum |a - b|, evaluated in row chunks to bound memory."""
    step = max(1, (4096 * 4096) // max(1, a.shape[1]))
    worst = 0.0
    for start in range(0, a.shape[0], step):
        end = min(start + step, a.shape[0])
        worst = max(worst, float(np.max(np.abs(a[start:end] - b[start:end]))))
    return worst


print(f"{'molecule':9s} {'basis':12s} {'n_AO':>7s} "
      f"{'tabula/s':>11s} {'veloxchem/s':>12s} {'ratio':>7s} {'max-diff':>11s}",
      flush=True)

for molecule_name in MOLECULES:
    molecule = read_molecule(molecule_name)
    for basis_label in BASES:
        try:
            basis = MolecularBasis.read(molecule, basis_label, ostream=None)

            tabula_time, tabula_matrix = best_time(
                lambda: TabulaOverlapDriver().compute(molecule, basis))
            veloxchem_time, veloxchem_matrix = best_time(
                lambda: OverlapDriver().compute(molecule, basis))

            # validation — Tabula vs VeloxChem
            veloxchem_full = veloxchem_matrix.full_matrix()
            del veloxchem_matrix
            gc.collect()

            tabula_values = tabula_matrix.to_numpy()
            veloxchem_values = veloxchem_full.to_numpy()
            n_ao = tabula_values.shape[0]
            difference = max_abs_difference(tabula_values, veloxchem_values)

            del tabula_values, veloxchem_values, tabula_matrix, veloxchem_full
            gc.collect()

            print(f"{molecule_name:9s} {basis_label:12s} {n_ao:7d} "
                  f"{tabula_time:11.4f} {veloxchem_time:12.4f} "
                  f"{veloxchem_time / tabula_time:7.2f} {difference:11.2e}",
                  flush=True)
        except Exception as error:
            print(f"{molecule_name:9s} {basis_label:12s}  FAILED: {error}",
                  flush=True)
            gc.collect()
