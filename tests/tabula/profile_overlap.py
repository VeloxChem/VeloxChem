"""Per-phase wall-time breakdown of tabula::OverlapDriver — both the dense
`compute` and the block-sparse `computeSparse` path.

Dense is profiled unscreened (threshold 0); sparse at the lossless 1e-12
threshold the benchmark uses. Times are thread-summed seconds.
"""

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.tabulalib import tabula_overlap_profile, tabula_overlap_profile_sparse

GEOMETRY_DIR = "/Users/rinkevic/Development/VeloxChem"
PHASES = ["make_blocks", "pair_setup", "screen", "kernel", "scatter", "symmetrize"]
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


def report(label, profile):
    total = sum(profile[p] for p in PHASES)
    print(f"  {label:8s}  total {total * 1e3:9.1f} ms")
    for phase in PHASES:
        ms = profile[phase] * 1e3
        share = ms / (total * 1e3) * 100 if total > 0 else 0.0
        print(f"    {phase:12s} {ms:9.2f} ms  ({share:5.1f} %)")


for molecule_name, basis_label in CELLS:
    molecule = read_molecule(molecule_name)
    basis = MolecularBasis.read(molecule, basis_label, ostream=None)

    tabula_overlap_profile(molecule, basis, 0.0)                          # warm-up
    dense = tabula_overlap_profile(molecule, basis, 0.0)

    tabula_overlap_profile_sparse(molecule, basis, SPARSE_THRESHOLD)      # warm-up
    sparse = tabula_overlap_profile_sparse(molecule, basis, SPARSE_THRESHOLD)

    print(f"\n{molecule_name} / {basis_label}")
    report("dense", dense)
    report("sparse", sparse)
