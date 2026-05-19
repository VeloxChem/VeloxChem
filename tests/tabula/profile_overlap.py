"""Per-phase wall-time breakdown of tabula::OverlapDriver::compute, plus the
per-section breakdown of the overlap seed.

Thread-summed seconds, across c240 and a span of def2 bases.
"""

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.tabulalib import (tabula_overlap_profile,
                                 seed_profile,
                                 reset_seed_profile)

GEOMETRY_DIR = "/Users/rinkevic/Development/VeloxChem"
PHASES = ["make_blocks", "pair_setup", "seed", "contract", "md", "transform", "scatter", "symmetrize"]


def read_molecule(name):
    with open(f"{GEOMETRY_DIR}/{name}.xyz") as handle:
        lines = handle.read().splitlines()
    return Molecule.read_str("\n".join(lines[2:]), "angstrom")


for molecule_name in ["c240"]:
    molecule = read_molecule(molecule_name)
    for basis_label in ["def2-SVP", "def2-TZVP", "def2-QZVP"]:
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        tabula_overlap_profile(molecule, basis, 0.0)            # warm-up

        reset_seed_profile()
        profile = tabula_overlap_profile(molecule, basis, 0.0)
        seed = seed_profile()

        total = sum(profile[p] for p in PHASES)
        print(f"\n{molecule_name} / {basis_label}   (thread-summed; total {total * 1e3:.1f} ms)")
        for phase in PHASES:
            ms = profile[phase] * 1e3
            print(f"  {phase:12s} {ms:9.2f} ms  ({ms / (total * 1e3) * 100:5.1f} %)")

        seed_total = profile["seed"] * 1e3
        print(f"  seed breakdown ({seed_total:.1f} ms):")
        for part in ("allocate", "row0", "ladder"):
            ms = seed[part] * 1e3
            print(f"    {part:10s}  {ms:9.2f} ms  ({ms / seed_total * 100:5.1f} %)")
