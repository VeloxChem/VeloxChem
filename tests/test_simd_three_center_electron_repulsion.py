import re

import numpy as np
import pytest

from veloxchem.veloxchemlib import AtomBasis, BasisFunction, MolecularBasis
from veloxchem.veloxchemlib import SimdThreeCenterElectronRepulsionDriver
from veloxchem.veloxchemlib import ThreeCenterElectronRepulsionDriver
from veloxchem.molecule import Molecule

# NOTE: the reference driver dispatches to angular momentum four on the bra and six
# on the auxiliary side, and returns zeros without an error above that. Comparing
# against it there would pass while checking nothing, so the sweep stops where it
# stops.
MAX_BRA_MOMENTUM = 4
MAX_AUX_MOMENTUM = 6

# NOTE: the sparse tensor and the flat buffer of the reference hold different sets of
# elements: the buffer keeps the upper triangle of an atom pair, the tensor its own
# blocks. Their sums therefore disagree by tens of per cent even when every integral
# matches, so this compares element by element and never totals.

LABELS = "spdfghikl"


class TestSimdThreeCenterElectronRepulsion:

    def one_function_basis(self, angular_momentum, identifier):
        """An atom basis of a single shell, so that (atom, momentum) fixes it."""

        atom_basis = AtomBasis()
        atom_basis.set_identifier(identifier)
        atom_basis.set_name("TEST")

        basis_function = BasisFunction()
        basis_function.set_angular_momentum(angular_momentum)
        basis_function.set_primitives([1.7 / (2.6**p) for p in range(2)],
                                      [1.0 / (p + 1) for p in range(2)])
        basis_function.normalize()
        atom_basis.add(basis_function)

        return atom_basis

    def dense_maps(self, basis, molecule):
        """(atom, momentum, component) -> row, and the components in dense order."""

        rows, per_shell = {}, {}

        for row, label in enumerate(basis.get_ao_basis_map(molecule)):
            fields = label.split()
            atom = int(fields[0]) - 1
            shell = re.match(r"\d+([spdfghikl])(.*)", fields[2])
            momentum = LABELS.index(shell.group(1))
            component = shell.group(2).strip()

            rows[(atom, momentum, component)] = row
            per_shell.setdefault((atom, momentum), []).append((row, component))

        order = {}
        for (atom, momentum), entries in per_shell.items():
            if momentum not in order:
                order[momentum] = [c for _, c in sorted(entries)]

        return rows, order

    def reference_slices(self, buffer, nao, naux):
        """The reference keeps the upper triangle, row major. Expand it."""

        upper = np.triu_indices(nao)
        slices = []

        for index in range(naux):
            matrix = np.zeros((nao, nao))
            matrix[upper] = np.asarray(buffer.values(index))
            slices.append(matrix + matrix.T - np.diag(np.diag(matrix)))

        return slices

    def worst_difference(self, molecule, identifiers, bra_momenta, aux_momentum):

        basis, aux_basis = MolecularBasis(), MolecularBasis()

        for identifier, momentum in zip(identifiers, bra_momenta):
            basis.add(self.one_function_basis(momentum, identifier))

        for identifier in identifiers:
            aux_basis.add(self.one_function_basis(aux_momentum, identifier))

        tensor = SimdThreeCenterElectronRepulsionDriver().compute(
            molecule, basis, aux_basis, 0.0)

        reference = ThreeCenterElectronRepulsionDriver().compute(
            molecule, basis, aux_basis)

        rows, order = self.dense_maps(basis, molecule)
        aux_rows, aux_order = self.dense_maps(aux_basis, molecule)

        slices = self.reference_slices(reference,
                                       basis.get_dimensions_of_basis(),
                                       aux_basis.get_dimensions_of_basis())

        worst, compared = 0.0, 0

        for index in range(tensor.number_of_blocks()):
            block = tensor.block(index)

            a_atoms, b_atoms = block.a_atoms(), block.b_atoms()
            c_atoms = block.c_atoms()

            if len(a_atoms) == 0 or len(c_atoms) == 0:
                continue

            values = tensor.block_to_numpy(index)

            la = basis.basis_sets()[block.a_index()].get_basis_functions()[
                0].get_angular_momentum()
            lb = basis.basis_sets()[block.b_index()].get_basis_functions()[
                0].get_angular_momentum()
            lc = aux_basis.basis_sets()[block.c_index()].get_basis_functions()[
                0].get_angular_momentum()

            npairs, natoms = len(a_atoms), len(c_atoms)
            nb, nc = 2 * lb + 1, 2 * lc + 1

            # NOTE: a combination is laid out as component, then atom on c side, then
            # atom pair, from its offset in the values of the block.
            base = block.element_offset(la, 0, lb, 0, lc, 0)

            bra = [[rows[(a, la, c)] for a in a_atoms] for c in order[la]]
            ket = [[rows[(b, lb, c)] for b in b_atoms] for c in order[lb]]
            aux = [[aux_rows[(c, lc, cc)] for c in c_atoms]
                   for cc in aux_order[lc]]

            for ma in range(2 * la + 1):
                for mb in range(nb):
                    for mc in range(nc):
                        offset = base + ((ma * nb + mb) * nc + mc) * natoms * npairs

                        for n in range(natoms):
                            got = values[offset + n * npairs:offset +
                                         (n + 1) * npairs]
                            want = slices[aux[mc][n]][bra[ma], ket[mb]]
                            worst = max(worst, float(np.max(np.abs(got - want))))
                            compared += npairs

        # NOTE: a comparison which visits nothing would pass, so the count is
        # returned and asserted on. Getting the layout wrong is the likely way to
        # visit nothing.
        return worst, compared

    @pytest.fixture
    def molecule(self):

        return Molecule.read_str(
            """O  0.00 0.00 0.00
               H  0.00 0.00 0.95
               N  2.60 0.30 0.10
               C  1.40 1.70 0.20""", "angstrom")

    def test_index_mapping_is_valid(self, molecule):
        """The sweep below is meaningless if the mapping is wrong, so prove it on the
        two simplest combinations first."""

        for bra in ((0, 0), (1, 0)):
            worst, compared = self.worst_difference(
                molecule, (8, 1, 7, 6), (bra[0], bra[0], bra[1], bra[1]), 0)

            assert compared > 0
            assert worst < 1.0e-10

    def test_against_reference(self, molecule):

        total = 0

        for la in range(MAX_BRA_MOMENTUM + 1):
            for lb in range(MAX_BRA_MOMENTUM + 1):
                for lc in range(MAX_AUX_MOMENTUM + 1):
                    worst, compared = self.worst_difference(
                        molecule, (8, 1, 7, 6), (la, la, lb, lb), lc)
                    total += compared

                    assert compared > 0, (
                        f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) compared nothing")

                    assert worst < 1.0e-10, (
                        f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) differs by {worst:.2e}")

        # NOTE: the sweep covered a million and a half integrals when it was written.
        # A large fall means the layout moved and most of them are no longer reached.
        assert total > 1000000
