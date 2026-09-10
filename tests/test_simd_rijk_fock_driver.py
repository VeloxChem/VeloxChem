import numpy as np
import pytest

from veloxchem.veloxchemlib import AtomBasis, BasisFunction, MolecularBasis
from veloxchem.veloxchemlib import PackedMatrix, mat_t
from veloxchem.veloxchemlib import SimdRIJFockDriver, SimdRIJKFockDriver
from veloxchem.veloxchemlib import SimdTwoCenterElectronRepulsionDriver
from veloxchem.molecule import Molecule

# NOTE: the driver is the free standing routines of CSimdRIJFockDriver held
# together, so it is checked against those composed by hand rather than against a
# reference of its own. What it adds is the state it keeps between the calls, the
# memory check, and the ranges the W matrices are formed in, and those are what the
# tests below are about.

LABELS = "spdfghikl"


def one_function_basis(angular_momentum, identifier):

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


class TestSimdRIJKFockDriver:

    @pytest.fixture
    def molecule(self):

        return Molecule.read_str(
            """O  0.00 0.00 0.00
               H  0.00 0.00 0.95
               N  2.60 0.30 0.10
               C  1.40 1.70 0.20""", "angstrom")

    @pytest.fixture
    def chain(self):
        """Enough atoms that the auxiliary basis outruns one range of W matrices."""

        lines = [f"H  0.00 0.00 {1.4 * i:.2f}" for i in range(10)]

        return Molecule.read_str("\n".join(lines), "angstrom")

    def bases(self, identifiers, bra_momenta, aux_momentum):

        basis, aux_basis = MolecularBasis(), MolecularBasis()

        for identifier, momentum in zip(identifiers, bra_momenta):
            basis.add(one_function_basis(momentum, identifier))

        for identifier in identifiers:
            aux_basis.add(one_function_basis(aux_momentum, identifier))

        return basis, aux_basis

    def orbitals(self, nao, norb, seed):

        rng = np.random.default_rng(seed)

        cmat = np.ascontiguousarray(rng.standard_normal((nao, norb)))

        coeffs = PackedMatrix(nao, norb, mat_t.general)
        coeffs.from_numpy(cmat)

        density = PackedMatrix(nao, nao, mat_t.symmetric)
        density.from_numpy(np.ascontiguousarray(cmat @ cmat.T))

        return coeffs, density

    def by_hand(self, molecule, basis, aux_basis, density, coeffs, factor, threshold):
        """The same Fock matrix from the free standing routines, with the W matrices
        formed in one range rather than in several."""

        naux = aux_basis.get_dimensions_of_basis()

        metric = SimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis)

        drv = SimdRIJFockDriver()

        bq = drv.compute_bq_vectors(molecule, basis, aux_basis,
                                    metric.cholesky_inverse(), threshold)

        fock = drv.compute_fock_matrix(bq, basis, aux_basis, density)

        if factor != 0.0:
            wvecs = drv.compute_w_vectors(bq, basis, aux_basis, coeffs, 0, naux)
            drv.compute_exchange_matrix(wvecs, fock, factor)

        return fock.to_numpy(max_memory=8.0)

    def test_against_the_routines_it_holds(self, molecule):

        for la, lb, lc in ((0, 0, 0), (1, 0, 1), (1, 1, 1), (2, 1, 2)):

            basis, aux_basis = self.bases((8, 1, 7, 6), (la, la, lb, lb), lc)

            nao = basis.get_dimensions_of_basis()

            coeffs, density = self.orbitals(nao, 3, 13 + nao)

            driver = SimdRIJKFockDriver()
            driver.prepare(molecule, basis, aux_basis, 0.0, 1 << 30)

            assert driver.is_prepared()

            for factor in (-1.0, 0.0, -0.25):

                computed = driver.compute(density, coeffs, factor).to_numpy(max_memory=8.0)

                expected = self.by_hand(molecule, basis, aux_basis, density, coeffs,
                                        factor, 0.0)

                assert np.array_equal(computed, computed.T)

                scale = max(float(np.max(np.abs(expected))), 1.0)

                assert np.max(np.abs(computed - expected)) / scale < 1.0e-12, (
                    f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) at factor {factor}")

    def test_exchange_factor_scales(self, molecule):
        """A hybrid functional scales the exchange, so the factor has to be a factor
        and not a sign."""

        basis, aux_basis = self.bases((8, 1, 7, 6), (1, 1, 0, 0), 1)

        nao = basis.get_dimensions_of_basis()

        coeffs, density = self.orbitals(nao, 3, 29)

        driver = SimdRIJKFockDriver()
        driver.prepare(molecule, basis, aux_basis, 0.0, 1 << 30)

        coulomb = driver.compute(density, coeffs, 0.0).to_numpy(max_memory=8.0)
        full = driver.compute(density, coeffs, -1.0).to_numpy(max_memory=8.0)
        quarter = driver.compute(density, coeffs, -0.25).to_numpy(max_memory=8.0)

        exchange = coulomb - full

        assert np.max(np.abs(coulomb)) > 0.0
        assert np.max(np.abs(exchange)) > 0.0

        assert np.allclose(quarter, coulomb - 0.25 * exchange, rtol=1.0e-12, atol=1.0e-13)

    def test_ranges_of_w_do_not_change_the_answer(self, chain):
        """The W matrices are formed a range at a time. With ten atoms and an f shell
        on the auxiliary side there are seventy of them, so the ranges wrap and the
        last one is shorter than the others."""

        basis, aux_basis = self.bases(tuple([1] * 10), tuple([0] * 10), 3)

        nao = basis.get_dimensions_of_basis()
        naux = aux_basis.get_dimensions_of_basis()

        assert naux > 64, f"the ranges do not wrap at {naux} auxiliary functions"

        coeffs, density = self.orbitals(nao, 4, 37)

        driver = SimdRIJKFockDriver()
        driver.prepare(chain, basis, aux_basis, 0.0, 1 << 30)

        computed = driver.compute(density, coeffs, -1.0).to_numpy(max_memory=8.0)

        expected = self.by_hand(chain, basis, aux_basis, density, coeffs, -1.0, 0.0)

        scale = float(np.max(np.abs(expected)))

        assert scale > 0.0
        assert np.max(np.abs(computed - expected)) / scale < 1.0e-12

    def test_required_memory_is_what_is_held(self, molecule):
        """The check is worth nothing if it is an estimate, so it is the count of the
        values the B vectors will hold."""

        basis, aux_basis = self.bases((8, 1, 7, 6), (1, 1, 0, 0), 2)

        driver = SimdRIJKFockDriver()

        needed = driver.required_memory(molecule, basis, aux_basis, 0.0)

        driver.prepare(molecule, basis, aux_basis, 0.0, 1 << 30)

        assert needed == driver.get_bq_vectors().number_of_elements() * 8

    def test_a_budget_which_is_too_small_is_refused(self, molecule):
        """It has to be an exception the caller can act on. A critical error would
        end the interpreter, which is what the check exists to avoid."""

        basis, aux_basis = self.bases((8, 1, 7, 6), (1, 1, 0, 0), 2)

        driver = SimdRIJKFockDriver()

        with pytest.raises(RuntimeError, match="B vectors need"):
            driver.prepare(molecule, basis, aux_basis, 0.0, 8)

        assert not driver.is_prepared()

    def test_compute_before_prepare_is_refused(self):

        driver = SimdRIJKFockDriver()

        assert not driver.is_prepared()
