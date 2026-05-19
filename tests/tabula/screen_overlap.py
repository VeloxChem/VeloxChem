"""Screened overlap — tabula::OverlapDriver at a range of thresholds.

For each (molecule, basis) it times the overlap at several screening
thresholds and reports the maximum element deviation from VeloxChem's
(unscreened) overlap — which a sound screener bounds by the threshold.

Run inside the VeloxChem environment.
"""

import gc
import time

import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.tabulalib import TabulaOverlapDriver

GEOMETRY_DIR = "/Users/rinkevic/Development/VeloxChem"
MOLECULES = ["c240"]
BASES = ["def2-TZVP", "def2-QZVP"]
THRESHOLDS = [0.0, 1.0e-14, 1.0e-12, 1.0e-10, 1.0e-8]
REPETITIONS = 5


def read_molecule(name):
    with open(f"{GEOMETRY_DIR}/{name}.xyz") as handle:
        lines = handle.read().splitlines()
    return Molecule.read_str("\n".join(lines[2:]), "angstrom")


def best_time(make):
    discard = make()
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
    step = max(1, (4096 * 4096) // max(1, a.shape[1]))
    worst = 0.0
    for start in range(0, a.shape[0], step):
        end = min(start + step, a.shape[0])
        worst = max(worst, float(np.max(np.abs(a[start:end] - b[start:end]))))
    return worst


print(f"{'molecule':9s} {'basis':12s} {'threshold':>11s} "
      f"{'tabula/s':>11s} {'max-dev':>11s}", flush=True)

for molecule_name in MOLECULES:
    molecule = read_molecule(molecule_name)
    for basis_label in BASES:
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        reference = OverlapDriver().compute(molecule, basis).full_matrix().to_numpy()

        for threshold in THRESHOLDS:
            t, matrix = best_time(
                lambda: TabulaOverlapDriver().compute(molecule, basis, threshold))
            dev = max_abs_difference(matrix.to_numpy(), reference)
            del matrix
            gc.collect()
            print(f"{molecule_name:9s} {basis_label:12s} {threshold:11.0e} "
                  f"{t:11.4f} {dev:11.2e}", flush=True)
        del reference
        gc.collect()
