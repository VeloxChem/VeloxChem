"""The range separated Coulomb integrals against pyscf.

    python rs_integral_check.py

Checks the two-center driver and the three-center one, each in both of its blocks,
against pyscf's `int2c2e` and `int3c2e` under `with_range_coulomb`. Prints a line
per case and exits non-zero if anything is outside tolerance.

Nothing in the repository computes attenuated Coulomb integrals other than the
drivers this checks, so pyscf is the only outside reference there is for them. It
is not a benchmark and takes seconds; it belongs here rather than in `tests`
because it needs pyscf, which the test suite does not require.

Three things about comparing with pyscf which are not obvious and which this file
exists to encode:

**The plain block is checked first.** It must be exactly what the unattenuated
driver gives -- same kernels, same path -- so any difference there is a mistake in
the range separated driver and not an accuracy question. Only then does the
attenuated block against pyscf mean anything.

**Magnitudes, not signed values.** For at least lithium's def2-universal-jkfit the
two codes differ by a sign convention, a diagonal plus or minus one similarity. The
spectrum, the trace and every magnitude agree to 1e-13 while the signed elements
differ by 20. Hydrogen, carbon and fluorine agree signed, so a case being clean
says nothing about the next one.

**Two different elements, for the three-center driver.** With one atom the product
center sits on the auxiliary center, so the Boys function's argument is zero and
the scaling of that argument goes untested -- which is the same blind spot that
sending omega to infinity has, theta going to one in that limit too.
"""
import sys

import numpy as np
import veloxchem as vlx
from pyscf import df, gto

from veloxchem.veloxchemlib import (SimdThreeCenterElectronRepulsionDriver,
                                    SimdThreeCenterElectronRepulsionRsDriver,
                                    SimdTwoCenterElectronRepulsionDriver,
                                    SimdTwoCenterElectronRepulsionRsDriver)

OMEGAS = (0.2, 0.8, 3.0)
FITTING = "def2-universal-jkfit"
THRESHOLD = 1.0e-14

# NOTE: the three-center cases are molecules of distinct elements alone, so a block
# is either a diagonal atom pair or an off-diagonal one and the two are told apart
# by the atom basis indices. A molecule with two atoms of one element would need the
# pairs within a block separated, which the block does not offer.
PAIRS = (("H", "F", 0.917, "def2-svp"),
         ("Li", "H", 1.595, "def2-svp"),
         ("H", "F", 0.917, "def2-tzvp"),
         ("C", "O", 1.128, "def2-tzvp"),
         ("B", "N", 1.281, "cc-pvtz"))

TWO_CENTER = ("def2-svp", "def2-universal-jkfit", "def2-tzvp")

WATER = [["O", (0.0, 0.0, 0.1173)], ["H", (0.0, 0.7572, -0.4692)],
         ["H", (0.0, -0.7572, -0.4692)]]

TOLERANCE = 1.0e-9


def _as_xyz(geometry):
    return f"{len(geometry)}\n\n" + "\n".join(
        f"{a}  {x} {y} {z}" for a, (x, y, z) in geometry)


def check_two_center(failures):
    """The two-center matrices, compared through quantities free of the ordering."""
    mol = vlx.Molecule.read_xyz_string(_as_xyz(WATER))
    pm = gto.M(atom=WATER, basis="def2-svp", unit="Angstrom")

    for basis in TWO_CENTER:
        bas = vlx.MolecularBasis.read(mol, basis.upper(), ostream=None)
        plain_only = SimdTwoCenterElectronRepulsionDriver().compute(mol, bas).to_numpy()

        # NOTE: pyscf is asked for the same basis, which for a fitting set means
        # asking the auxiliary molecule for its own two-center integrals.
        if "fit" in basis:
            ref_mol = df.addons.make_auxmol(pm, basis)
        else:
            ref_mol = gto.M(atom=WATER, basis=basis, unit="Angstrom")

        for omega in OMEGAS:
            coulomb, attenuated = SimdTwoCenterElectronRepulsionRsDriver().compute(
                mol, bas, omega)
            coulomb, attenuated = coulomb.to_numpy(), attenuated.to_numpy()

            exact = float(np.abs(coulomb - plain_only).max())

            with ref_mol.with_range_coulomb(omega):
                reference = ref_mol.intor("int2c2e")

            ours = np.sort(np.abs(attenuated).ravel())
            theirs = np.sort(np.abs(reference).ravel())
            scale = max(float(np.abs(theirs).max()), 1.0)
            against = float(np.abs(ours - theirs).max()) / scale

            ok = exact == 0.0 and against < TOLERANCE
            if not ok:
                failures.append(f"two-center {basis} omega {omega}")
            print(f"  {'ok  ' if ok else 'FAIL'} two-center  {basis:22s} omega {omega:<4g} "
                  f" plain block exact {exact:.1e}   attenuated vs pyscf {against:.1e}")


def _magnitudes(tensor, want):
    """The magnitudes of a sparse tensor, padded to the dense count and sorted.

    The dense array holds an off-diagonal atom pair both ways round while the tensor
    holds it once, so those blocks are repeated. The screening drops elements the
    dense array still carries; their values are below the threshold, so counting
    them as zero costs at most that.
    """
    parts = []
    for iblock in range(tensor.number_of_blocks()):
        block = tensor.block(iblock)
        npairs = block.number_of_pairs()
        values = np.abs(np.asarray(tensor.block_to_numpy(iblock)))

        # NOTE: the pattern keeps the upper triangle of the **atom pairs**, which is
        # not the upper triangle of the atom basis pairs. A block whose two sides are
        # the same atom basis can still hold pairs of two different atoms: water's
        # hydrogens give one block with H1-H1, H2-H2 and H1-H2 in it, and only the
        # last of those stands for two. Repeating whole blocks is right for a
        # molecule of distinct elements, which is all this file measures, and wrong
        # for anything else. The atom pairs are the fastest varying index, so the
        # columns can be repeated individually.
        columns = values.reshape(-1, npairs)

        for k, (a, b) in enumerate(zip(block.a_atoms(), block.b_atoms())):
            parts.append(columns[:, k].ravel())
            if a != b:
                parts.append(columns[:, k].ravel())

    values = np.concatenate(parts)

    if values.size < want:
        values = np.concatenate([values, np.zeros(want - values.size)])

    return np.sort(values)


def check_three_center(failures):
    """The two tensors, compared as magnitudes against pyscf's dense array."""
    plain_drv = SimdThreeCenterElectronRepulsionDriver()
    rs_drv = SimdThreeCenterElectronRepulsionRsDriver()

    for a_elem, b_elem, distance, orbital in PAIRS:
        geometry = [[a_elem, (0.0, 0.0, 0.0)], [b_elem, (0.0, 0.0, distance)]]

        mol = vlx.Molecule.read_xyz_string(_as_xyz(geometry))
        bas = vlx.MolecularBasis.read(mol, orbital.upper(), ostream=None)
        aux = vlx.MolecularBasis.read(mol, FITTING.upper(), ostream=None)

        pm = gto.M(atom=geometry, basis=orbital, unit="Angstrom")
        am = df.addons.make_auxmol(pm, FITTING)

        plain_only = plain_drv.compute(mol, bas, aux, THRESHOLD)
        reference = np.sort(np.abs(df.incore.aux_e2(pm, am, intor="int3c2e").ravel()))

        for omega in OMEGAS:
            coulomb, attenuated = rs_drv.compute(mol, bas, aux, THRESHOLD, omega)

            exact = max(
                float(np.abs(np.asarray(coulomb.block_to_numpy(i))
                             - np.asarray(plain_only.block_to_numpy(i))).max())
                for i in range(coulomb.number_of_blocks()))

            with pm.with_range_coulomb(omega):
                att_ref = np.sort(np.abs(df.incore.aux_e2(pm, am, intor="int3c2e").ravel()))

            scale = max(float(reference.max()), 1.0)
            plain_vs = float(np.abs(_magnitudes(coulomb, reference.size) - reference).max()) / scale
            att_vs = float(np.abs(_magnitudes(attenuated, att_ref.size) - att_ref).max()) / scale

            ok = exact == 0.0 and plain_vs < TOLERANCE and att_vs < TOLERANCE
            if not ok:
                failures.append(f"three-center {a_elem}{b_elem} {orbital} omega {omega}")
            print(f"  {'ok  ' if ok else 'FAIL'} three-center {a_elem}{b_elem:2s} {orbital:9s} "
                  f"omega {omega:<4g}  plain block exact {exact:.1e}   "
                  f"plain vs pyscf {plain_vs:.1e}   attenuated vs pyscf {att_vs:.1e}")


def check_limits(failures):
    """What the operator does at the two ends, which no reference is needed for."""
    geometry = [["C", (0.0, 0.0, 0.0)], ["O", (0.0, 0.0, 1.128)]]
    mol = vlx.Molecule.read_xyz_string(_as_xyz(geometry))
    bas = vlx.MolecularBasis.read(mol, "DEF2-SVP", ostream=None)
    aux = vlx.MolecularBasis.read(mol, FITTING.upper(), ostream=None)

    rs_drv = SimdThreeCenterElectronRepulsionRsDriver()

    _, at_zero = rs_drv.compute(mol, bas, aux, THRESHOLD, 0.0)
    empty = max(float(np.abs(np.asarray(at_zero.block_to_numpy(i))).max())
                for i in range(at_zero.number_of_blocks()))

    coulomb, at_large = rs_drv.compute(mol, bas, aux, THRESHOLD, 1.0e6)
    apart = max(float(np.abs(np.asarray(coulomb.block_to_numpy(i))
                             - np.asarray(at_large.block_to_numpy(i))).max())
                for i in range(coulomb.number_of_blocks()))

    ok = empty == 0.0 and apart < 1.0e-9
    if not ok:
        failures.append("the limits of omega")
    print(f"  {'ok  ' if ok else 'FAIL'} limits      omega zero leaves {empty:.1e}   "
          f"omega large differs from the Coulomb tensor by {apart:.1e}")


def main():
    failures = []

    check_two_center(failures)
    check_three_center(failures)
    check_limits(failures)

    print("\nall checks passed" if not failures else "\nFAILURES:")
    for failure in failures:
        print("  " + failure)

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
