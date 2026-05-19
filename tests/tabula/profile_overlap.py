"""Per-phase wall-time breakdown of tabula::OverlapDriver::compute, across
c240 and a span of def2 bases.

Thread-summed seconds.
"""

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.tabulalib import tabula_overlap_profile

GEOMETRY_DIR = "/Users/rinkevic/Development/VeloxChem"
PHASES = ["make_blocks", "pair_setup", "kernel", "scatter", "symmetrize"]


def read_molecule(name):
    with open(f"{GEOMETRY_DIR}/{name}.xyz") as handle:
        lines = handle.read().splitlines()
    return Molecule.read_str("\n".join(lines[2:]), "angstrom")


for molecule_name in ["c240"]:
    molecule = read_molecule(molecule_name)
    for basis_label in ["def2-SVP", "def2-TZVP", "def2-QZVP"]:
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        tabula_overlap_profile(molecule, basis, 0.0)            # warm-up
        profile = tabula_overlap_profile(molecule, basis, 0.0)

        total = sum(profile[p] for p in PHASES)
        print(f"\n{molecule_name} / {basis_label}   (thread-summed; total {total * 1e3:.1f} ms)")
        for phase in PHASES:
            ms = profile[phase] * 1e3
            print(f"  {phase:12s} {ms:9.2f} ms  ({ms / (total * 1e3) * 100:5.1f} %)")
