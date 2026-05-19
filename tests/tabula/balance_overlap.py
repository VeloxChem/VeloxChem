"""Per-thread load balance of tabula::OverlapDriver's block-pair parallel
loop, across c240 and a span of def2 bases.

The driver parallelizes one #pragma omp parallel for, schedule(dynamic),
over the triangular list of basis-function-block pairs — a coarse unit of
work whose cost varies by orders of magnitude.
"""

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.tabulalib import tabula_overlap_profile, overlap_thread_balance

GEOMETRY_DIR = "/Users/rinkevic/Development/VeloxChem"


def read_molecule(name):
    with open(f"{GEOMETRY_DIR}/{name}.xyz") as handle:
        lines = handle.read().splitlines()
    return Molecule.read_str("\n".join(lines[2:]), "angstrom")


for molecule_name in ["c240"]:
    molecule = read_molecule(molecule_name)
    for basis_label in ["def2-SVP", "def2-TZVP", "def2-QZVP"]:
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        tabula_overlap_profile(molecule, basis, 0.0)          # warm-up
        tabula_overlap_profile(molecule, basis, 0.0)          # measured
        balance = overlap_thread_balance()

        busy = balance["busy"]
        pairs = balance["pairs"]
        wall = balance["wall"]
        nthreads = len(busy)
        n_pairs = sum(pairs)

        busy_max = max(busy)
        busy_min = min(busy)
        busy_mean = sum(busy) / nthreads
        # ideal wall = perfectly balanced; speedup vs the sum of busy time
        ideal_wall = sum(busy) / nthreads
        utilisation = sum(busy) / (wall * nthreads)

        print(f"\n{molecule_name} / {basis_label}   "
              f"{n_pairs} block pairs, {nthreads} threads")
        print(f"  region wall      {wall * 1e3:8.2f} ms")
        print(f"  busy  min/mean/max {busy_min * 1e3:7.2f} /"
              f"{busy_mean * 1e3:7.2f} /{busy_max * 1e3:7.2f} ms")
        print(f"  imbalance        max/mean = {busy_max / busy_mean:5.2f}x   "
              f"ideal wall {ideal_wall * 1e3:.2f} ms vs actual {wall * 1e3:.2f} ms")
        print(f"  utilisation      {utilisation * 100:5.1f} %  "
              f"(busy time / (wall x threads))")
        print("  per thread (busy ms : block pairs):")
        for t in range(nthreads):
            bar = "#" * int(round(40 * busy[t] / busy_max)) if busy_max > 0 else ""
            print(f"    t{t:02d}  {busy[t] * 1e3:8.2f} : {pairs[t]:3d}  {bar}")
