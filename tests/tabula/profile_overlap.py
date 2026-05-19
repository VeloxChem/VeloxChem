"""Per-phase wall-time breakdown of tabula::OverlapDriver::compute, plus the
per-section breakdown of the single-centre MD recursion.

Thread-summed seconds, across c240 and a span of def2 bases.
"""

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.tabulalib import (tabula_overlap_profile,
                                 md_recursion_profile,
                                 reset_md_recursion_profile)

GEOMETRY_DIR = "/Users/rinkevic/Development/VeloxChem"
PHASES = ["make_blocks", "pair_setup", "seed", "contract", "md", "transform", "scatter"]


def read_molecule(name):
    with open(f"{GEOMETRY_DIR}/{name}.xyz") as handle:
        lines = handle.read().splitlines()
    return Molecule.read_str("\n".join(lines[2:]), "angstrom")


for molecule_name in ["c240"]:
    molecule = read_molecule(molecule_name)
    for basis_label in ["def2-SVP", "def2-TZVP", "def2-QZVP"]:
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        tabula_overlap_profile(molecule, basis, 0.0)            # warm-up

        reset_md_recursion_profile()
        profile = tabula_overlap_profile(molecule, basis, 0.0)
        md = md_recursion_profile()

        total = sum(profile[p] for p in PHASES)
        print(f"\n{molecule_name} / {basis_label}   (thread-summed; total {total * 1e3:.1f} ms)")
        for phase in PHASES:
            ms = profile[phase] * 1e3
            print(f"  {phase:12s} {ms:9.2f} ms  ({ms / (total * 1e3) * 100:5.1f} %)")

        md_total = profile["md"] * 1e3
        print(f"  md breakdown ({md_total:.1f} ms):")
        print(f"    allocate    {md['allocate'] * 1e3:9.2f} ms  "
              f"({md['allocate'] * 1e5 / md_total:5.1f} %)")
        print(f"    seed_copy   {md['seed_copy'] * 1e3:9.2f} ms  "
              f"({md['seed_copy'] * 1e5 / md_total:5.1f} %)")
        for degree, seconds in enumerate(md["build_by_degree"]):
            if seconds > 0.0:
                ms = seconds * 1e3
                print(f"    build d={degree}  {ms:9.2f} ms  ({ms / md_total * 100:5.1f} %)")
