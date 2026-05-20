"""Charge-dipole field: sparse-density vs dense-density path.

For each (molecule, basis) it reports, best-of-N:
  - dense/ms  — compute_field over a dense density (every atom-pair),
  - sparse/ms — compute_field_sparse over the same density in block-sparse
    storage (atom-pairs whose density block is dropped are skipped),
  - dense/sparse — the speedup,
  - stored — the sparse density's footprint vs the dense dimension^2,
  - max-diff — the two fields agree (a correctness check).

A nuclear-attraction matrix in block-sparse storage stands in for the AO density
(the field contraction is linear in D, and this carries the same bra-ket spatial
sparsity a real density has). The external points are the nuclei, capped so the
(shared) point count does not dominate — the speedup is the atom-pair skip and
is independent of the point count. Run inside the VeloxChem environment.
"""

import gc
import time

import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.tabulalib import TabulaNuclearAttractionDriver, TabulaChargeDipoleDriver

GEOMETRY = {
    "c60":       "/Users/rinkevic/Development/VeloxChem/c60.xyz",
    "taxol":     "/Users/rinkevic/Development/VeloxChem/taxol.xyz",
    "olestra":   "/Users/rinkevic/Development/VeloxChem/olestra.xyz",
    "ubiquitin": "/Users/rinkevic/Downloads/ubiquitin.xyz",
}

PLAN = [
    ("c60",       ["def2-svp", "def2-tzvp"]),
    ("taxol",     ["def2-svp", "def2-tzvp"]),
    ("olestra",   ["def2-svp"]),
    ("ubiquitin", ["def2-svp"]),
]

SCREEN_THRESHOLD = 1.0e-10  # lossless density screen
MAX_POINTS = 128


def read_molecule(path):
    with open(path) as handle:
        lines = handle.read().splitlines()
    return Molecule.read_str("\n".join(lines[2:]), "angstrom")


def reps_for(n):
    return 3 if n < 2000 else (2 if n < 6000 else 1)


def measure(make, reps):
    make(); gc.collect()
    best = float("inf")
    for _ in range(reps):
        gc.collect()
        start = time.perf_counter()
        result = make()
        best = min(best, time.perf_counter() - start)
        del result; gc.collect()
    return best


def main():
    nuc = TabulaNuclearAttractionDriver()
    cdp = TabulaChargeDipoleDriver()
    print("| molecule  | basis     | n_AO | n_pts | dense/ms | sparse/ms | dense/sparse | stored | max-diff |")
    print("|-----------|-----------|-----:|------:|---------:|----------:|-------------:|-------:|---------:|")
    for name, bases in PLAN:
        mol = read_molecule(GEOMETRY[name])
        coords_all = [list(c) for c in np.array(mol.get_coordinates_in_bohr())]
        stride = max(1, len(coords_all) // MAX_POINTS)
        coords = coords_all[::stride][:MAX_POINTS]
        npts = len(coords)
        for basis_label in bases:
            bas = MolecularBasis.read(mol, basis_label, ostream=None)

            d_sparse = nuc.compute_sparse(mol, bas, SCREEN_THRESHOLD)
            dim = d_sparse.dimension()
            stored = d_sparse.stored_element_count() / (dim * dim)
            d_dense = d_sparse.to_dense()

            reps = reps_for(dim)
            t_dense = measure(lambda: cdp.compute_field(mol, bas, d_dense, coords, 0.0), reps)
            t_sparse = measure(lambda: cdp.compute_field_sparse(mol, bas, d_sparse, coords, 0.0), reps)

            e_dense = cdp.compute_field(mol, bas, d_dense, coords, 0.0)
            e_sparse = cdp.compute_field_sparse(mol, bas, d_sparse, coords, 0.0)
            diff = float(np.max(np.abs(e_dense - e_sparse)))

            print(f"| {name:9s} | {basis_label:9s} | {dim:5d} | {npts:5d} | {t_dense*1e3:7.1f} | "
                  f"{t_sparse*1e3:8.1f} | {t_dense/t_sparse:11.2f}x | {stored*100:4.1f}% | {diff:.0e} |",
                  flush=True)


if __name__ == "__main__":
    main()
