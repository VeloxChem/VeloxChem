"""How tight the three-center screening is, against the values it keeps.

    python t3c_screening_check.py [molecule.xyz] [basis] [aux] [thresholds]

The two halves of every comparison in this file screen differently. The four-center
build uses Cauchy-Schwarz, `sqrt((mn|mn)) sqrt((ls|ls))`, which is built from computed
integral magnitudes and is tight by construction. The three-center driver uses an
analytic primitive bound, `ScreeningFunc.hpp:209`, built from the **smallest exponent**
of each contracted function and the distance between the two orbital atoms:

    fmu = a.smallest_exponent() * b.smallest_exponent() / fexp
    ... * exp(-fmu * r * r)

Taking the most diffuse primitive of a contracted function makes that bound loose:
it keeps combinations whose integrals are far below the threshold it was given. How
loose has never been measured, and everything downstream inherits it -- the density of
the B vectors, the padding of the gathered chunk, and the exponent the fitted path is
compared at.

**What this measures is the tightness of the bound, not a head to head with
Cauchy-Schwarz.** It builds the tensor at a threshold and asks what fraction of the
values that survived are actually above it. A tight bound keeps little that is
negligible; a loose one keeps a lot, and that fraction is work which is computed,
stored, and contracted for nothing.

It is deliberately not a timing, so it needs no quiet machine.
"""
import sys

import numpy as np
import veloxchem as vlx
from veloxchem.veloxchemlib import SimdThreeCenterElectronRepulsionDriver

DEFAULT_THRESHOLDS = [1.0e-10, 1.0e-12, 1.0e-14]


def magnitudes(tensor):
    """Every stored magnitude, as one array."""
    parts = []

    for iblock in range(tensor.number_of_blocks()):
        parts.append(np.abs(np.asarray(tensor.block_to_numpy(iblock))).ravel())

    return np.concatenate(parts) if parts else np.zeros(0)


def main():
    xyz = sys.argv[1] if len(sys.argv) > 1 else '/Users/rinkevic/Downloads/005.xyz'

    basis_name = sys.argv[2] if len(sys.argv) > 2 else 'def2-svp'

    aux_name = sys.argv[3] if len(sys.argv) > 3 else 'def2-universal-jkfit'

    thresholds = ([float(t) for t in sys.argv[4].split(',')]
                  if len(sys.argv) > 4 else DEFAULT_THRESHOLDS)

    molecule = vlx.Molecule.read_xyz_file(xyz)

    basis = vlx.MolecularBasis.read(molecule, basis_name.upper(), ostream=None)

    aux_basis = vlx.MolecularBasis.read(molecule, aux_name.upper(), ostream=None)

    print(f"# {xyz}  {basis_name} / {aux_name}")
    print(f"# {molecule.number_of_atoms()} atoms, nao {basis.get_dimensions_of_basis()}, "
          f"naux {aux_basis.get_dimensions_of_basis()}")
    print()
    print(f"{'threshold':>10s} {'values kept':>13s} {'above it':>9s} {'above /10':>10s} "
          f"{'above /100':>11s} {'median':>10s}")

    drv = SimdThreeCenterElectronRepulsionDriver()

    previous = None

    for threshold in thresholds:
        tensor = drv.compute(molecule, basis, aux_basis, threshold)

        values = magnitudes(tensor)

        if values.size == 0:
            print(f"{threshold:10.1e}  nothing kept")
            continue

        # NOTE: the fraction of what was kept which is worth keeping. A bound which
        # is tight keeps little below the threshold it was given; this one is being
        # asked how much of its own output it need not have produced.
        share = lambda cut: 100.0 * float((values >= cut).sum()) / values.size

        print(f"{threshold:10.1e} {values.size:13d} {share(threshold):8.1f} % "
              f"{share(threshold / 10):9.1f} % {share(threshold / 100):10.1f} % "
              f"{np.median(values):10.1e}")

        previous = (threshold, values.size)

    print()
    print("'above it' is the share of the stored values which reach the threshold the")
    print("pattern was built at. The rest is computed, stored and contracted for")
    print("nothing: it is the looseness of the bound, and it is the same fraction the")
    print("B vectors, the padding and the contraction all carry.")

    if previous is not None:
        print()
        print("NOTE: the screening decides per combination of basis functions and atom")
        print("pair, not per value, so a block which is worth keeping carries small")
        print("values with it. This is therefore an upper bound on what a tighter")
        print("bound of the same granularity could remove, not a promise.")


if __name__ == '__main__':
    sys.exit(main())
