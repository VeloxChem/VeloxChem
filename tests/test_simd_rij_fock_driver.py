import re

import numpy as np
import pytest

from veloxchem.veloxchemlib import AtomBasis, BasisFunction, MolecularBasis
from veloxchem.veloxchemlib import SimdRIJFockDriver
from veloxchem.veloxchemlib import SimdThreeCenterElectronRepulsionDriver
from veloxchem.veloxchemlib import SimdTwoCenterElectronRepulsionDriver
from veloxchem.molecule import Molecule

# NOTE: the B vectors are B(q)_ij = sum over p of (ij|p) Linv_pq, so they are checked
# against that contraction formed in numpy from the three-center integrals and the
# inverted metric. The layout of the sparse tensor is the one proven in
# test_simd_three_center_electron_repulsion.py, and the expansion below is shared by
# the integrals and the B vectors, which carry q where the integrals carry p.

LABELS = "spdfghikl"


class TestSimdRIJFockDriver:

    def one_function_basis(self, angular_momentum, identifier):

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

    def expand(self, tensor, basis, aux_basis, maps, aux_maps, nao, naux):
        """The sparse tensor as a dense (nao, nao, naux) array."""

        rows, order = maps
        aux_rows, aux_order = aux_maps

        dense = np.zeros((nao, nao, naux))
        visited = 0

        for index in range(tensor.number_of_blocks()):
            block = tensor.block(index)

            a_atoms, b_atoms, c_atoms = block.a_atoms(), block.b_atoms(), block.c_atoms()

            if len(a_atoms) == 0 or len(c_atoms) == 0:
                continue

            values = tensor.block_to_numpy(index)

            la = basis.basis_sets()[block.a_index()].get_basis_functions()[0].get_angular_momentum()
            lb = basis.basis_sets()[block.b_index()].get_basis_functions()[0].get_angular_momentum()
            lc = aux_basis.basis_sets()[block.c_index()].get_basis_functions()[0].get_angular_momentum()

            npairs = block.number_of_pairs(la, 0, lb, 0, lc, 0)
            natoms = len(c_atoms)
            nb, nc = 2 * lb + 1, 2 * lc + 1

            base = block.element_offset(la, 0, lb, 0, lc, 0)

            bra = [[rows[(a, la, c)] for a in a_atoms] for c in order[la]]
            ket = [[rows[(b, lb, c)] for b in b_atoms] for c in order[lb]]
            aux = [[aux_rows[(c, lc, cc)] for c in c_atoms] for cc in aux_order[lc]]

            for ma in range(2 * la + 1):
                for mb in range(nb):
                    for mc in range(nc):
                        offset = base + ((ma * nb + mb) * nc + mc) * natoms * npairs

                        for n in range(natoms):
                            vals = values[offset + n * npairs:offset + (n + 1) * npairs]
                            dense[bra[ma], ket[mb], aux[mc][n]] = vals
                            dense[ket[mb], bra[ma], aux[mc][n]] = vals
                            visited += npairs

        return dense, visited

    def run_case(self, molecule, identifiers, bra_momenta, aux_momentum):

        basis, aux_basis = MolecularBasis(), MolecularBasis()

        for identifier, momentum in zip(identifiers, bra_momenta):
            basis.add(self.one_function_basis(momentum, identifier))

        for identifier in identifiers:
            aux_basis.add(self.one_function_basis(aux_momentum, identifier))

        nao = basis.get_dimensions_of_basis()
        naux = aux_basis.get_dimensions_of_basis()

        maps = self.dense_maps(basis, molecule)
        aux_maps = self.dense_maps(aux_basis, molecule)

        # the inverted metric, in the packed format the driver takes

        metric = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis)
        inverse = metric.invert()

        # the three-center integrals, and the reference contraction in numpy

        integrals = SimdThreeCenterElectronRepulsionDriver().compute(
            molecule, basis, aux_basis, 0.0)

        dense_ints, _ = self.expand(integrals, basis, aux_basis, maps, aux_maps, nao, naux)

        expected = np.einsum('ijp,pq->ijq', dense_ints, inverse.to_numpy())

        # the B vectors

        bq = SimdRIJFockDriver().compute_bq_vectors(
            molecule, basis, aux_basis, inverse, 0.0)

        computed, visited = self.expand(bq, basis, aux_basis, maps, aux_maps, nao, naux)

        scale = max(float(np.max(np.abs(expected))), 1.0)

        return float(np.max(np.abs(computed - expected))) / scale, visited

    @pytest.fixture
    def molecule(self):

        return Molecule.read_str(
            """O  0.00 0.00 0.00
               H  0.00 0.00 0.95
               N  2.60 0.30 0.10
               C  1.40 1.70 0.20""", "angstrom")

    def test_simplest_combinations(self, molecule):
        """The sweep is meaningless if the mapping is wrong, so prove it first."""

        for bra in ((0, 0), (1, 0)):
            worst, visited = self.run_case(molecule, (8, 1, 7, 6),
                                           (bra[0], bra[0], bra[1], bra[1]), 0)
            assert visited > 0
            assert worst < 1.0e-10

    def test_against_contraction(self, molecule):

        total = 0

        for la in range(4):
            for lb in range(4):
                for lc in range(5):
                    worst, visited = self.run_case(molecule, (8, 1, 7, 6),
                                                   (la, la, lb, lb), lc)
                    total += visited

                    assert visited > 0, (
                        f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) visited nothing")

                    assert worst < 1.0e-10, (
                        f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) differs by {worst:.2e}")

        # NOTE: the sweep visited close to a quarter of a million elements when it
        # was written. A large fall means the layout moved and most of the B vectors
        # are no longer reached, which would let the comparison pass while checking
        # almost nothing.
        assert total > 200000
