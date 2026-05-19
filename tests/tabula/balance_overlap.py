"""Per-thread load balance of tabula::OverlapDriver's task parallel loop —
for both the dense `compute` and the block-sparse `computeSparse` path.

The driver parallelizes one #pragma omp parallel for, schedule(dynamic),
over tasks — each task is one block pair restricted to a bra-CGTO range,
sized to a roughly uniform contracted-pair count. `overlap_thread_balance`
reports the most recent run.
"""

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.tabulalib import (tabula_overlap_profile,
                                 tabula_overlap_profile_sparse,
                                 overlap_thread_balance)

GEOMETRY_DIR = "/Users/rinkevic/Development/VeloxChem"
SPARSE_THRESHOLD = 1.0e-12

# (molecule, basis) cells — a compact span plus the large cells.
CELLS = [
    ("c240", "def2-SVP"),
    ("c240", "def2-TZVP"),
    ("c240", "def2-QZVP"),
    ("olestra", "def2-QZVPD"),
    ("olestra", "aug-cc-pVQZ"),
]


def read_molecule(name):
    with open(f"{GEOMETRY_DIR}/{name}.xyz") as handle:
        lines = handle.read().splitlines()
    return Molecule.read_str("\n".join(lines[2:]), "angstrom")


def report(label, balance):
    busy = balance["busy"]
    pairs = balance["pairs"]
    wall = balance["wall"]
    nthreads = len(busy)
    n_pairs = sum(pairs)

    busy_max = max(busy)
    busy_min = min(busy)
    busy_mean = sum(busy) / nthreads
    ideal_wall = sum(busy) / nthreads
    utilisation = sum(busy) / (wall * nthreads) if wall > 0 else 0.0

    print(f"  {label}  —  {n_pairs} tasks, {nthreads} threads")
    print(f"    region wall      {wall * 1e3:8.2f} ms")
    print(f"    busy  min/mean/max {busy_min * 1e3:7.2f} /"
          f"{busy_mean * 1e3:7.2f} /{busy_max * 1e3:7.2f} ms")
    print(f"    imbalance        max/mean = {busy_max / busy_mean:5.2f}x   "
          f"ideal wall {ideal_wall * 1e3:.2f} ms vs actual {wall * 1e3:.2f} ms")
    print(f"    utilisation      {utilisation * 100:5.1f} %  "
          f"(busy time / (wall x threads))")
    print("    per thread (busy ms : tasks):")
    for t in range(nthreads):
        bar = "#" * int(round(40 * busy[t] / busy_max)) if busy_max > 0 else ""
        print(f"      t{t:02d}  {busy[t] * 1e3:8.2f} : {pairs[t]:3d}  {bar}")


for molecule_name, basis_label in CELLS:
    molecule = read_molecule(molecule_name)
    basis = MolecularBasis.read(molecule, basis_label, ostream=None)

    print(f"\n{molecule_name} / {basis_label}")

    tabula_overlap_profile(molecule, basis, 0.0)                      # warm-up
    tabula_overlap_profile(molecule, basis, 0.0)                      # measured
    report("dense", overlap_thread_balance())

    tabula_overlap_profile_sparse(molecule, basis, SPARSE_THRESHOLD)  # warm-up
    tabula_overlap_profile_sparse(molecule, basis, SPARSE_THRESHOLD)  # measured
    report("sparse", overlap_thread_balance())
