import subprocess
import sys
import textwrap

import numpy as np
import pytest

from veloxchem.veloxchemlib import AtomBasis, BasisFunction, MolecularBasis
from veloxchem.veloxchemlib import NuclearPotentialDriver
from veloxchem.veloxchemlib import SimdNuclearPotentialDriver
from veloxchem.molecule import Molecule

# NOTE: the kernels reach angular momentum six. A combination above that stops with
# an error rather than returning the zeros of a kernel which does not exist, which a
# caller could not tell from integrals which are genuinely zero.
#
# NOTE: the plain CNuclearPotentialDriver, which is what these are checked against,
# returns **all zeros for h blocks** -- measured against PySCF on 2026-09-12, water
# in cc-pV5Z. It is therefore a reference only up to g, and the tests which compare
# stop there. The h and i kernels are checked for what can be checked without a
# reference: that they run, that the matrix is symmetric, and that they are not
# themselves zero.
#
# NOTE: the s only case has one angular component, so it exercises neither the
# transformation nor the angular coupling of the diagonal blocks -- and the nuclear
# potential, unlike the overlap, has no closed form there and computes them with the
# same kernels. Every case above s is what puts the driver through its paces.


def atom_basis_of(exponents, coefficients, identifier, momenta):
    """An atom basis with one contracted function of each angular momentum given."""

    atom_basis = AtomBasis()
    atom_basis.set_identifier(identifier)
    atom_basis.set_name("TEST")

    for momentum in momenta:
        basis_function = BasisFunction()
        basis_function.set_angular_momentum(momentum)
        basis_function.set_primitives(exponents, coefficients)
        basis_function.normalize()
        atom_basis.add(basis_function)

    return atom_basis


class TestSimdNuclearPotentialDriver:

    @pytest.fixture
    def molecule(self):

        return Molecule.read_str(
            """O  0.00 0.00 0.00
               H  0.00 0.00 0.95
               N  2.60 0.30 0.10
               C  1.40 1.70 0.20""", "angstrom")

    @pytest.fixture(params=[(0,), (1,), (2,), (3,), (4,), (0, 1), (2, 3), (0, 2, 4)],
                    ids=["s", "p", "d", "f", "g", "sp", "df", "sdg"])
    def basis(self, request):
        """Bases over the kernels the plain driver can be a reference for, singly and
        in combination. The mixed ones are what reach the off-diagonal kernels: an
        s and g basis dispatches sg and gs as well as ss and gg."""

        basis = MolecularBasis()

        for identifier, exponents in ((8, [7.5, 1.3]), (1, [3.4, 0.7]),
                                      (7, [6.1, 1.1]), (6, [4.9, 0.9])):
            basis.add(atom_basis_of(exponents, [1.0, 0.6], identifier,
                                    request.param))

        return basis

    def reference(self, molecule, basis, charges=None, points=None):
        """The matrix the plain driver gives, which is the one to match."""

        driver = NuclearPotentialDriver()

        if charges is None:
            return driver.compute(molecule, basis).to_numpy()

        coords = [list(points[3 * i:3 * i + 3]) for i in range(len(charges))]

        return driver.compute(molecule, basis, charges, coords).to_numpy()

    def test_the_nuclei_of_the_molecule(self, molecule, basis):
        """The form which takes no charges uses the nuclei of the molecule, and
        has to give what the plain driver gives for the same."""

        computed = SimdNuclearPotentialDriver().compute(molecule, basis)

        expected = self.reference(molecule, basis)

        scale = float(np.max(np.abs(expected)))

        assert scale > 0.0
        assert np.max(np.abs(computed.to_numpy(basis) - expected)) / scale < 1.0e-11

    def test_charges_which_are_not_the_nuclei(self, molecule, basis):
        """The form which takes charges is what an embedding or a set of external
        charges needs. The positions are a flat array of three coordinates each."""

        charges = [0.7, -0.4, 1.3]

        points = [0.3, -0.2, 1.1, 2.0, 1.0, -0.5, -1.4, 0.8, 0.25]

        computed = SimdNuclearPotentialDriver().compute(molecule, basis, charges,
                                                        points)

        expected = self.reference(molecule, basis, charges, points)

        scale = float(np.max(np.abs(expected)))

        assert scale > 0.0
        assert np.max(np.abs(computed.to_numpy(basis) - expected)) / scale < 1.0e-11

    def test_the_matrix_is_symmetric(self, molecule, basis):

        computed = SimdNuclearPotentialDriver().compute(molecule, basis).to_numpy(basis)

        assert np.array_equal(computed, computed.T)

    def test_no_charges_gives_no_potential(self, molecule, basis):
        """An empty set of charges is not an error: there is simply no operator,
        and the integrals are zero rather than undefined."""

        computed = SimdNuclearPotentialDriver().compute(molecule, basis, [], [])

        assert np.max(np.abs(computed.to_numpy(basis))) == 0.0

    def test_a_kernel_which_is_not_written_is_refused(self):
        """The kernels reach angular momentum six. A combination above it must stop
        rather than return the zeros of a kernel which does not exist -- a caller
        cannot tell those from integrals which are genuinely zero. This test moves up
        as the kernels are added; k functions are the first which are not written.

        It is checked in a process of its own because a critical error terminates
        the interpreter rather than raising, so pytest.raises cannot see it."""

        script = textwrap.dedent("""
            from veloxchem.veloxchemlib import AtomBasis, BasisFunction
            from veloxchem.veloxchemlib import MolecularBasis
            from veloxchem.veloxchemlib import SimdNuclearPotentialDriver
            from veloxchem.molecule import Molecule

            mol = Molecule.read_str("O 0.0 0.0 0.0\\nH 0.0 0.0 0.95", "angstrom")

            basis = MolecularBasis()
            for identifier in (8, 1):
                atom_basis = AtomBasis()
                atom_basis.set_identifier(identifier)
                atom_basis.set_name("TEST")
                function = BasisFunction()
                function.set_angular_momentum(7)
                function.set_primitives([2.1, 0.5], [1.0, 0.5])
                function.normalize()
                atom_basis.add(function)
                basis.add(atom_basis)

            SimdNuclearPotentialDriver().compute(mol, basis)
            print("NOT REFUSED")
            """)

        outcome = subprocess.run([sys.executable, "-c", script],
                                 capture_output=True, text=True)

        assert outcome.returncode != 0, "a kernel which is not written was not refused"

        assert "NOT REFUSED" not in outcome.stdout

        # NOTE: either the table of buffer rows or the dispatch refuses, whichever
        # the caller reaches first -- the arena is sized before a kernel is chosen,
        # so it is usually the table. Both are the refusal being asked for here.

        said = outcome.stdout + outcome.stderr

        assert ("Angular momentum is out of range" in said) or (
            "No kernel for the combination of angular momenta" in said), said[:400]

    @pytest.mark.parametrize("momentum, label", [(5, "h"), (6, "i")])
    def test_above_the_reference(self, molecule, momentum, label):
        """The plain driver returns all zeros for h blocks, so there is nothing to
        compare against at h and i. What can be asserted without a reference is
        asserted: that the kernels run, that the matrix is symmetric, and that it is
        not itself zero -- which is what a missing kernel would give.

        This is the gap worth closing: the SCF takes its core Hamiltonian from the
        plain drivers, so a basis with h functions gets zeros there today and
        converges quietly to an energy below the true one."""

        basis = MolecularBasis()

        for identifier, exponents in ((8, [7.5, 1.3]), (1, [3.4, 0.7]),
                                      (7, [6.1, 1.1]), (6, [4.9, 0.9])):
            basis.add(atom_basis_of(exponents, [1.0, 0.6], identifier, (momentum,)))

        computed = SimdNuclearPotentialDriver().compute(molecule,
                                                        basis).to_numpy(basis)

        ncomps = 2 * momentum + 1

        assert computed.shape == (4 * ncomps, 4 * ncomps)

        assert np.array_equal(computed, computed.T)

        assert np.max(np.abs(computed)) > 0.0, (
            f"the {label} kernels returned nothing but zeros")

        # NOTE: and the plain driver really does give zeros here, which is what makes
        # the comparison above impossible rather than merely omitted.

        if momentum == 5:
            reference = NuclearPotentialDriver().compute(molecule,
                                                         basis).to_numpy()

            assert np.max(np.abs(reference)) == 0.0, (
                "the plain driver no longer returns zeros at h, so these blocks can "
                "and should be compared against it")
