import re

import numpy as np
import pytest

from veloxchem.veloxchemlib import AtomBasis, BasisFunction, MolecularBasis
from veloxchem.veloxchemlib import PackedMatrix, mat_t
from veloxchem.veloxchemlib import SimdRIFockDriver
from veloxchem.veloxchemlib import SimdThreeCenterElectronRepulsionDriver
from veloxchem.veloxchemlib import SimdTwoCenterElectronRepulsionDriver
from veloxchem.molecule import Molecule

# NOTE: the B vectors are B(q)_ij = sum over p of (ij|p) Linv_pq, so they are checked
# against that contraction formed in numpy from the three-center integrals and the
# inverted metric. The layout of the sparse tensor is the one proven in
# test_simd_three_center_electron_repulsion.py, and the expansion below is shared by
# the integrals and the B vectors, which carry q where the integrals carry p.

LABELS = "spdfghikl"


class TestSimdRIFockDriver:

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

        bq = SimdRIFockDriver().compute_bq_vectors(
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

    def y_case(self, molecule, identifiers, bra_momenta, aux_momentum, kind):
        """Y(q) = sum over i and j of B(q)_ij D_ij, against the same in numpy."""

        basis, aux_basis = MolecularBasis(), MolecularBasis()

        for identifier, momentum in zip(identifiers, bra_momenta):
            basis.add(self.one_function_basis(momentum, identifier))

        for identifier in identifiers:
            aux_basis.add(self.one_function_basis(aux_momentum, identifier))

        nao = basis.get_dimensions_of_basis()
        naux = aux_basis.get_dimensions_of_basis()

        maps = self.dense_maps(basis, molecule)
        aux_maps = self.dense_maps(aux_basis, molecule)

        inverse = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis).invert()

        drv = SimdRIFockDriver()
        bq = drv.compute_bq_vectors(molecule, basis, aux_basis, inverse, 0.0)

        dense_bq, visited = self.expand(bq, basis, aux_basis, maps, aux_maps, nao, naux)

        # the density: symmetric for the field, genuinely asymmetric for response

        rng = np.random.default_rng(17 + nao)
        raw = rng.standard_normal((nao, nao))

        if kind == 'symmetric':
            dmat = raw + raw.T
            packed = PackedMatrix(nao, nao, mat_t.symmetric)
        else:
            dmat = raw
            packed = PackedMatrix(nao, nao, mat_t.general)

        packed.from_numpy(np.ascontiguousarray(dmat))

        expected = np.einsum('ijq,ij->q', dense_bq, dmat)

        computed = np.asarray(drv.compute_y_vector(bq, basis, aux_basis, packed))

        scale = max(float(np.max(np.abs(expected))), 1.0)

        return float(np.max(np.abs(computed - expected))) / scale, visited, computed

    def test_y_vector_symmetric_density(self, molecule):
        """The density of a self consistent field calculation."""

        total = 0

        for la in range(4):
            for lb in range(4):
                for lc in range(4):
                    worst, visited, _ = self.y_case(molecule, (8, 1, 7, 6),
                                                    (la, la, lb, lb), lc, 'symmetric')
                    total += visited

                    assert visited > 0
                    assert worst < 1.0e-10, (
                        f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) differs by {worst:.2e}")

        # NOTE: the sweep visited close to two hundred thousand elements of the B
        # vectors when it was written. A large fall means the layout moved and the
        # contraction no longer reaches most of them.
        assert total > 150000

    def test_y_vector_general_density(self, molecule):
        """The density of a response calculation, which is not symmetric."""

        total = 0

        for la in range(4):
            for lb in range(4):
                for lc in range(4):
                    worst, visited, _ = self.y_case(molecule, (8, 1, 7, 6),
                                                    (la, la, lb, lb), lc, 'general')
                    total += visited

                    assert visited > 0
                    assert worst < 1.0e-10, (
                        f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) differs by {worst:.2e}")

        # NOTE: the sweep visited close to two hundred thousand elements of the B
        # vectors when it was written. A large fall means the layout moved and the
        # contraction no longer reaches most of them.
        assert total > 150000

    def test_general_path_reproduces_symmetric(self, molecule):
        """A symmetric density must give the same Y down either path."""

        basis, aux_basis = MolecularBasis(), MolecularBasis()

        for identifier, momentum in zip((8, 1, 7, 6), (1, 1, 0, 0)):
            basis.add(self.one_function_basis(momentum, identifier))

        for identifier in (8, 1, 7, 6):
            aux_basis.add(self.one_function_basis(1, identifier))

        nao = basis.get_dimensions_of_basis()

        inverse = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis).invert()

        drv = SimdRIFockDriver()
        bq = drv.compute_bq_vectors(molecule, basis, aux_basis, inverse, 0.0)

        rng = np.random.default_rng(5)
        raw = rng.standard_normal((nao, nao))
        dmat = np.ascontiguousarray(raw + raw.T)

        sym = PackedMatrix(nao, nao, mat_t.symmetric)
        sym.from_numpy(dmat)

        gen = PackedMatrix(nao, nao, mat_t.general)
        gen.from_numpy(dmat)

        y_sym = np.asarray(drv.compute_y_vector(bq, basis, aux_basis, sym))
        y_gen = np.asarray(drv.compute_y_vector(bq, basis, aux_basis, gen))

        assert np.max(np.abs(y_sym)) > 0.0
        assert np.allclose(y_sym, y_gen, rtol=1.0e-12, atol=1.0e-12)

    def test_resolution_of_identity_closes(self, molecule):
        """The property the whole chain rests on: with L inverted, where the metric
        is L L transposed, the B vectors satisfy

            sum over q of B(q)_ij B(q)_kl = sum over p, t of (ij|p) Jinv_pt (kl|t)

        which is what makes F = sum over q of B(q) Y(q) the Coulomb matrix. The full
        inverse of the metric does not have this property, and the second assertion
        below records that, so the two are never confused again.
        """

        compared = 0

        for la, lb, lc in ((0, 0, 0), (1, 0, 1), (1, 1, 2), (2, 1, 1)):

            basis, aux_basis = MolecularBasis(), MolecularBasis()

            for identifier, momentum in zip((8, 1, 7, 6), (la, la, lb, lb)):
                basis.add(self.one_function_basis(momentum, identifier))

            for identifier in (8, 1, 7, 6):
                aux_basis.add(self.one_function_basis(lc, identifier))

            nao = basis.get_dimensions_of_basis()
            naux = aux_basis.get_dimensions_of_basis()

            maps = self.dense_maps(basis, molecule)
            aux_maps = self.dense_maps(aux_basis, molecule)

            metric = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis)

            integrals = SimdThreeCenterElectronRepulsionDriver().compute(
                molecule, basis, aux_basis, 0.0)

            dense_ints, _ = self.expand(integrals, basis, aux_basis, maps, aux_maps, nao, naux)

            target = np.einsum('ijp,pt,klt->ijkl', dense_ints,
                               np.linalg.inv(metric.to_numpy(max_memory=8.0)), dense_ints)

            scale = float(np.max(np.abs(target)))

            drv = SimdRIFockDriver()

            # both forms of the metric must close: the inverted Cholesky factor,
            # and the inverted square root which drops the directions of the small
            # eigenvalues and is what a nearly linearly dependent basis needs

            for metric_matrix in (metric.cholesky_inverse(),
                                  metric.inverse_square_root(1.0e-12)):

                bq = drv.compute_bq_vectors(molecule, basis, aux_basis,
                                            metric_matrix, 0.0)

                dense_bq, visited = self.expand(bq, basis, aux_basis, maps, aux_maps,
                                                nao, naux)

                closes = np.einsum('ijq,klq->ijkl', dense_bq, dense_bq)

                compared += visited

                assert np.max(np.abs(closes - target)) / scale < 1.0e-11, (
                    f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) does not close")

            # the full inverse is the wrong matrix here, and wrong by order one
            # rather than subtly, which is what makes this worth asserting

            wrong = drv.compute_bq_vectors(molecule, basis, aux_basis, metric.invert(), 0.0)

            dense_wrong, _ = self.expand(wrong, basis, aux_basis, maps, aux_maps, nao, naux)

            missed = np.einsum('ijq,klq->ijkl', dense_wrong, dense_wrong)

            assert np.max(np.abs(missed - target)) / scale > 1.0e-2

        assert compared > 0

    def fock_case(self, molecule, identifiers, bra_momenta, aux_momentum, kind):
        """The whole chain: L inverted, B vectors, Y vector, Coulomb matrix."""

        basis, aux_basis = MolecularBasis(), MolecularBasis()

        for identifier, momentum in zip(identifiers, bra_momenta):
            basis.add(self.one_function_basis(momentum, identifier))

        for identifier in identifiers:
            aux_basis.add(self.one_function_basis(aux_momentum, identifier))

        nao = basis.get_dimensions_of_basis()
        naux = aux_basis.get_dimensions_of_basis()

        maps = self.dense_maps(basis, molecule)
        aux_maps = self.dense_maps(aux_basis, molecule)

        metric = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis)

        integrals = SimdThreeCenterElectronRepulsionDriver().compute(
            molecule, basis, aux_basis, 0.0)

        dense_ints, _ = self.expand(integrals, basis, aux_basis, maps, aux_maps, nao, naux)

        rng = np.random.default_rng(23 + nao)
        raw = rng.standard_normal((nao, nao))

        if kind == 'symmetric':
            dmat = np.ascontiguousarray(raw + raw.T)
            packed = PackedMatrix(nao, nao, mat_t.symmetric)
        else:
            dmat = np.ascontiguousarray(raw)
            packed = PackedMatrix(nao, nao, mat_t.general)

        packed.from_numpy(dmat)

        drv = SimdRIFockDriver()

        bq = drv.compute_bq_vectors(molecule, basis, aux_basis,
                                    metric.cholesky_inverse(), 0.0)

        dense_bq, visited = self.expand(bq, basis, aux_basis, maps, aux_maps, nao, naux)

        fock = drv.compute_fock_matrix(bq, basis, aux_basis, packed)

        assert fock.get_type() == mat_t.symmetric

        computed = fock.to_numpy(max_memory=8.0)

        # the Coulomb matrix is symmetric to the last bit, as one triangle is stored

        assert np.array_equal(computed, computed.T)

        # the reference, built independently from the raw integrals and the inverse

        reference = np.einsum('ijp,pt,klt,kl->ij', dense_ints,
                              np.linalg.inv(metric.to_numpy(max_memory=8.0)),
                              dense_ints, dmat)

        scale = max(float(np.max(np.abs(reference))), 1.0)

        return float(np.max(np.abs(computed - reference))) / scale, visited

    def test_fock_matrix_symmetric_density(self, molecule):
        """The Coulomb matrix against the resolution of the identity written out."""

        total = 0

        for la in range(4):
            for lb in range(4):
                for lc in range(4):
                    worst, visited = self.fock_case(molecule, (8, 1, 7, 6),
                                                    (la, la, lb, lb), lc, 'symmetric')
                    total += visited

                    assert visited > 0
                    assert worst < 1.0e-11, (
                        f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) differs by {worst:.2e}")

        # NOTE: the sweep visited close to two hundred thousand elements of the B
        # vectors when it was written, and a large fall means most of them are no
        # longer reached.
        assert total > 150000

    def test_fock_matrix_general_density(self, molecule):
        """A response density, which is not symmetric, still gives a symmetric matrix."""

        total = 0

        for la in range(4):
            for lb in range(4):
                for lc in range(4):
                    worst, visited = self.fock_case(molecule, (8, 1, 7, 6),
                                                    (la, la, lb, lb), lc, 'general')
                    total += visited

                    assert visited > 0
                    assert worst < 1.0e-11, (
                        f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) differs by {worst:.2e}")

        # NOTE: the sweep visited close to two hundred thousand elements of the B
        # vectors when it was written, and a large fall means most of them are no
        # longer reached.
        assert total > 150000

    def test_fock_matrix_from_y_vector(self, molecule):
        """The two forms must agree, as one is written in terms of the other."""

        basis, aux_basis = MolecularBasis(), MolecularBasis()

        for identifier, momentum in zip((8, 1, 7, 6), (1, 1, 0, 0)):
            basis.add(self.one_function_basis(momentum, identifier))

        for identifier in (8, 1, 7, 6):
            aux_basis.add(self.one_function_basis(1, identifier))

        nao = basis.get_dimensions_of_basis()

        metric = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis)

        drv = SimdRIFockDriver()
        bq = drv.compute_bq_vectors(molecule, basis, aux_basis,
                                    metric.cholesky_inverse(), 0.0)

        rng = np.random.default_rng(31)
        raw = rng.standard_normal((nao, nao))
        dmat = np.ascontiguousarray(raw + raw.T)

        packed = PackedMatrix(nao, nao, mat_t.symmetric)
        packed.from_numpy(dmat)

        yvec = drv.compute_y_vector(bq, basis, aux_basis, packed)

        from_density = drv.compute_fock_matrix(bq, basis, aux_basis, packed).to_numpy()
        from_yvector = drv.compute_fock_matrix(bq, basis, aux_basis, yvec).to_numpy()

        assert np.max(np.abs(from_density)) > 0.0

        # NOTE: not bit for bit. The blocks are handed to the threads dynamically,
        # so which of them a thread sums varies between runs and the last bit of the
        # total varies with it. The spread is one unit in the last place, which the
        # tolerance below is far tighter than while still catching a real difference.

        assert np.allclose(from_density, from_yvector, rtol=1.0e-13, atol=1.0e-14)

    def w_case(self, molecule, identifiers, bra_momenta, aux_momentum, nocc):
        """W(q)_is = sum over r of B(q)_ir C_rs, against the same in numpy."""

        basis, aux_basis = MolecularBasis(), MolecularBasis()

        for identifier, momentum in zip(identifiers, bra_momenta):
            basis.add(self.one_function_basis(momentum, identifier))

        for identifier in identifiers:
            aux_basis.add(self.one_function_basis(aux_momentum, identifier))

        nao = basis.get_dimensions_of_basis()
        naux = aux_basis.get_dimensions_of_basis()

        maps = self.dense_maps(basis, molecule)
        aux_maps = self.dense_maps(aux_basis, molecule)

        metric = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis)

        drv = SimdRIFockDriver()

        bq = drv.compute_bq_vectors(molecule, basis, aux_basis,
                                    metric.cholesky_inverse(), 0.0)

        dense_bq, _ = self.expand(bq, basis, aux_basis, maps, aux_maps, nao, naux)

        # NOTE: the columns of the coefficients are different from one another, so
        # that a contribution which is dropped cannot be hidden by another.

        rng = np.random.default_rng(41 + nao)
        cmat = np.ascontiguousarray(rng.standard_normal((nao, nocc)))

        packed = PackedMatrix(nao, nocc, mat_t.general)
        packed.from_numpy(cmat)

        wvecs = drv.compute_w_vectors(bq, basis, aux_basis, packed, 0, naux)

        assert len(wvecs) == naux

        computed = np.stack([w.to_numpy(max_memory=8.0) for w in wvecs])

        expected = np.einsum('irq,rs->qis', dense_bq, cmat)

        # every row of W must be reached. The B vectors keep one of the two orders
        # of an off-diagonal pair of atoms, and the transformation has to add into
        # the rows of both sides of it, so a row left at zero means half of the
        # contributions were dropped.

        rows_touched = int(np.sum(np.max(np.abs(computed), axis=(0, 2)) > 0.0))

        scale = max(float(np.max(np.abs(expected))), 1.0)

        return (float(np.max(np.abs(computed - expected))) / scale,
                computed.size, rows_touched, nao)

    def test_w_vectors_against_transformation(self, molecule):

        total = 0

        for la in range(4):
            for lb in range(4):
                for lc in range(4):
                    worst, visited, touched, nao = self.w_case(
                        molecule, (8, 1, 7, 6), (la, la, lb, lb), lc, 3)
                    total += visited

                    assert touched == nao, (
                        f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) reached "
                        f"{touched} of {nao} rows")

                    assert worst < 1.0e-12, (
                        f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) differs by {worst:.2e}")

        # NOTE: the sweep visited close to fifty thousand elements of W when it was
        # written, and a large fall means most of them are no longer reached.
        assert total > 40000

    def test_w_vectors_range(self, molecule):
        """A range of the auxiliary basis must be the same as that slice of all."""

        basis, aux_basis = MolecularBasis(), MolecularBasis()

        for identifier, momentum in zip((8, 1, 7, 6), (1, 1, 0, 0)):
            basis.add(self.one_function_basis(momentum, identifier))

        for identifier in (8, 1, 7, 6):
            aux_basis.add(self.one_function_basis(1, identifier))

        nao = basis.get_dimensions_of_basis()
        naux = aux_basis.get_dimensions_of_basis()

        metric = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis)

        drv = SimdRIFockDriver()
        bq = drv.compute_bq_vectors(molecule, basis, aux_basis,
                                    metric.cholesky_inverse(), 0.0)

        rng = np.random.default_rng(7)
        cmat = np.ascontiguousarray(rng.standard_normal((nao, 4)))

        packed = PackedMatrix(nao, 4, mat_t.general)
        packed.from_numpy(cmat)

        whole = drv.compute_w_vectors(bq, basis, aux_basis, packed, 0, naux)

        first, last = 2, min(7, naux)

        part = drv.compute_w_vectors(bq, basis, aux_basis, packed, first, last)

        assert len(part) == last - first
        assert np.max(np.abs(whole[first].to_numpy())) > 0.0

        for j in range(last - first):
            assert np.array_equal(part[j].to_numpy(), whole[first + j].to_numpy())

    def exchange_case(self, molecule, identifiers, bra_momenta, aux_momentum, nocc):
        """K_ij = sum over q and s of W(q)_is W(q)_js, against the same in numpy."""

        basis, aux_basis = MolecularBasis(), MolecularBasis()

        for identifier, momentum in zip(identifiers, bra_momenta):
            basis.add(self.one_function_basis(momentum, identifier))

        for identifier in identifiers:
            aux_basis.add(self.one_function_basis(aux_momentum, identifier))

        nao = basis.get_dimensions_of_basis()
        naux = aux_basis.get_dimensions_of_basis()

        metric = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis)

        drv = SimdRIFockDriver()

        bq = drv.compute_bq_vectors(molecule, basis, aux_basis,
                                    metric.cholesky_inverse(), 0.0)

        rng = np.random.default_rng(53 + nao)
        cmat = np.ascontiguousarray(rng.standard_normal((nao, nocc)))

        packed = PackedMatrix(nao, nocc, mat_t.general)
        packed.from_numpy(cmat)

        wvecs = drv.compute_w_vectors(bq, basis, aux_basis, packed, 0, naux)

        dense_w = np.stack([w.to_numpy(max_memory=8.0) for w in wvecs])

        expected = np.einsum('qis,qjs->ij', dense_w, dense_w)

        matrix = PackedMatrix(nao, nao, mat_t.symmetric)
        matrix.zero()

        drv.compute_exchange_matrix(wvecs, matrix)

        computed = matrix.to_numpy(max_memory=8.0)

        assert np.array_equal(computed, computed.T)

        scale = max(float(np.max(np.abs(expected))), 1.0)

        return float(np.max(np.abs(computed - expected))) / scale, computed.size

    def test_exchange_against_contraction(self, molecule):

        total = 0

        for la in range(3):
            for lb in range(3):
                for lc in range(3):
                    worst, visited = self.exchange_case(
                        molecule, (8, 1, 7, 6), (la, la, lb, lb), lc, 3)
                    total += visited

                    assert visited > 0
                    assert worst < 1.0e-12, (
                        f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) differs by {worst:.2e}")

        assert total > 2000

    def test_exchange_accumulates(self, molecule):
        """It adds to the matrix rather than replacing it, applies the factor, and
        gives the same result whether the auxiliary basis is taken in one range or
        in several. All three are what a calculation relies on."""

        basis, aux_basis = MolecularBasis(), MolecularBasis()

        for identifier, momentum in zip((8, 1, 7, 6), (1, 1, 0, 0)):
            basis.add(self.one_function_basis(momentum, identifier))

        for identifier in (8, 1, 7, 6):
            aux_basis.add(self.one_function_basis(1, identifier))

        nao = basis.get_dimensions_of_basis()
        naux = aux_basis.get_dimensions_of_basis()

        metric = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis)

        drv = SimdRIFockDriver()
        bq = drv.compute_bq_vectors(molecule, basis, aux_basis,
                                    metric.cholesky_inverse(), 0.0)

        rng = np.random.default_rng(61)
        cmat = np.ascontiguousarray(rng.standard_normal((nao, 3)))

        packed = PackedMatrix(nao, 3, mat_t.general)
        packed.from_numpy(cmat)

        whole = drv.compute_w_vectors(bq, basis, aux_basis, packed, 0, naux)

        once = PackedMatrix(nao, nao, mat_t.symmetric)
        once.zero()
        drv.compute_exchange_matrix(whole, once)
        single = once.to_numpy(max_memory=8.0)

        assert np.max(np.abs(single)) > 0.0

        twice = PackedMatrix(nao, nao, mat_t.symmetric)
        twice.zero()
        drv.compute_exchange_matrix(whole, twice)
        drv.compute_exchange_matrix(whole, twice)

        assert np.allclose(twice.to_numpy(), 2.0 * single, rtol=1.0e-13, atol=1.0e-14)

        scaled = PackedMatrix(nao, nao, mat_t.symmetric)
        scaled.zero()
        drv.compute_exchange_matrix(whole, scaled, -0.5)

        assert np.allclose(scaled.to_numpy(), -0.5 * single, rtol=1.0e-13, atol=1.0e-14)

        # the same auxiliary basis, taken in two ranges

        split = PackedMatrix(nao, nao, mat_t.symmetric)
        split.zero()
        middle = naux // 2
        drv.compute_exchange_matrix(
            drv.compute_w_vectors(bq, basis, aux_basis, packed, 0, middle), split)
        drv.compute_exchange_matrix(
            drv.compute_w_vectors(bq, basis, aux_basis, packed, middle, naux), split)

        assert np.allclose(split.to_numpy(), single, rtol=1.0e-12, atol=1.0e-13)

    def test_coulomb_and_exchange_compose(self, molecule):
        """The exchange is added onto the matrix the Coulomb build produced."""

        basis, aux_basis = MolecularBasis(), MolecularBasis()

        for identifier, momentum in zip((8, 1, 7, 6), (1, 1, 0, 0)):
            basis.add(self.one_function_basis(momentum, identifier))

        for identifier in (8, 1, 7, 6):
            aux_basis.add(self.one_function_basis(1, identifier))

        nao = basis.get_dimensions_of_basis()
        naux = aux_basis.get_dimensions_of_basis()

        metric = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis)

        drv = SimdRIFockDriver()
        bq = drv.compute_bq_vectors(molecule, basis, aux_basis,
                                    metric.cholesky_inverse(), 0.0)

        nocc = 3
        rng = np.random.default_rng(67)
        cmat = np.ascontiguousarray(rng.standard_normal((nao, nocc)))

        coeffs = PackedMatrix(nao, nocc, mat_t.general)
        coeffs.from_numpy(cmat)

        # the density of those orbitals, which is what a field calculation carries

        dmat = np.ascontiguousarray(cmat @ cmat.T)
        density = PackedMatrix(nao, nao, mat_t.symmetric)
        density.from_numpy(dmat)

        fock = drv.compute_fock_matrix(bq, basis, aux_basis, density)
        coulomb = fock.to_numpy(max_memory=8.0).copy()

        wvecs = drv.compute_w_vectors(bq, basis, aux_basis, coeffs, 0, naux)
        drv.compute_exchange_matrix(wvecs, fock, -1.0)

        combined = fock.to_numpy(max_memory=8.0)

        dense_w = np.stack([w.to_numpy(max_memory=8.0) for w in wvecs])
        exchange = np.einsum('qis,qjs->ij', dense_w, dense_w)

        assert np.max(np.abs(coulomb)) > 0.0
        assert np.max(np.abs(exchange)) > 0.0
        assert np.array_equal(combined, combined.T)

        scale = float(np.max(np.abs(combined)))

        assert np.max(np.abs(combined - (coulomb - exchange))) / scale < 1.0e-12

    def test_both_forms_of_the_transformation_agree(self, molecule):
        """W is formed either by walking the values of the B vectors or by expanding
        them into a square and handing that to a matrix product. Which is taken
        follows from how dense they are, and the two must give the same matrices.
        The product is the default, so without this the sum would go unexercised."""

        total = 0

        for la, lb, lc in ((0, 0, 0), (1, 0, 1), (1, 1, 1), (2, 1, 2)):

            basis, aux_basis = MolecularBasis(), MolecularBasis()

            for identifier, momentum in zip((8, 1, 7, 6), (la, la, lb, lb)):
                basis.add(self.one_function_basis(momentum, identifier))

            for identifier in (8, 1, 7, 6):
                aux_basis.add(self.one_function_basis(lc, identifier))

            nao = basis.get_dimensions_of_basis()
            naux = aux_basis.get_dimensions_of_basis()

            metric = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis)

            drv = SimdRIFockDriver()

            # the default takes the product, so the sum is what needs asking for

            assert drv.get_dense_threshold() == 0.0

            bq = drv.compute_bq_vectors(molecule, basis, aux_basis,
                                        metric.cholesky_inverse(), 0.0)

            rng = np.random.default_rng(89 + nao)
            cmat = np.ascontiguousarray(rng.standard_normal((nao, 4)))

            coeffs = PackedMatrix(nao, 4, mat_t.general)
            coeffs.from_numpy(cmat)

            matrices = {}

            for tag, threshold in (('product', 0.0), ('sum', 2.0)):
                drv.set_dense_threshold(threshold)
                wvecs = drv.compute_w_vectors(bq, basis, aux_basis, coeffs, 0, naux)
                matrices[tag] = np.stack([w.to_numpy(max_memory=8.0) for w in wvecs])

            scale = max(float(np.max(np.abs(matrices['sum']))), 1.0)

            assert np.max(np.abs(matrices['sum'])) > 0.0

            assert np.max(np.abs(matrices['product'] - matrices['sum'])) / scale < 1.0e-12, (
                f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) the two forms disagree")

            total += matrices['sum'].size

        # NOTE: the four combinations visited a couple of thousand elements of W
        # when this was written, over both forms of the transformation.
        assert total > 2000
