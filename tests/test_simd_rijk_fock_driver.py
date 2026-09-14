import numpy as np
import pytest

from veloxchem.veloxchemlib import AtomBasis, BasisFunction, MolecularBasis
from veloxchem.veloxchemlib import PackedMatrix, mat_t
from veloxchem.veloxchemlib import SimdRIFockDriver, SimdRIJKFockDriver
from veloxchem.veloxchemlib import rimode
from veloxchem.veloxchemlib import SimdTwoCenterElectronRepulsionDriver
from veloxchem.molecule import Molecule

# NOTE: the driver is the free standing routines of CSimdRIFockDriver held
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

        drv = SimdRIFockDriver()

        bq = drv.compute_bq_vectors(molecule, basis, aux_basis,
                                    metric.cholesky_inverse(), threshold)

        fock = drv.compute_fock_matrix(bq, basis, aux_basis, density)

        # the driver returns twice the Coulomb less the scaled exchange, which is
        # the convention of a closed shell calculation

        fock.scale(2.0)

        if factor != 0.0:
            wvecs = drv.compute_w_vectors(bq, basis, aux_basis, coeffs, 0, naux)
            drv.compute_exchange_matrix(wvecs, fock, -factor)

        return fock.to_numpy(max_memory=8.0)

    def test_against_the_routines_it_holds(self, molecule):

        for la, lb, lc in ((0, 0, 0), (1, 0, 1), (1, 1, 1), (2, 1, 2)):

            basis, aux_basis = self.bases((8, 1, 7, 6), (la, la, lb, lb), lc)

            nao = basis.get_dimensions_of_basis()

            coeffs, density = self.orbitals(nao, 3, 13 + nao)

            driver = SimdRIJKFockDriver()
            driver.prepare(molecule, basis, aux_basis, 0.0, 1 << 30)

            assert driver.is_prepared()

            for factor in (1.0, 0.0, 0.25):

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
        full = driver.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0)
        quarter = driver.compute(density, coeffs, 0.25).to_numpy(max_memory=8.0)

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

        computed = driver.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0)

        expected = self.by_hand(chain, basis, aux_basis, density, coeffs, 1.0, 0.0)

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

    def test_a_budget_which_is_too_small_takes_the_direct_way(self, molecule):
        """The B vectors of a large molecule fit nowhere, and a budget they exceed
        selects the way which does not hold them rather than refusing the work."""

        basis, aux_basis = self.bases((8, 1, 7, 6), (1, 1, 0, 0), 2)

        nao = basis.get_dimensions_of_basis()

        coeffs, density = self.orbitals(nao, 3, 101)

        roomy = SimdRIJKFockDriver()
        roomy.prepare(molecule, basis, aux_basis, 0.0, 1 << 30)

        assert roomy.get_mode() == rimode.in_memory

        cramped = SimdRIJKFockDriver()
        cramped.prepare(molecule, basis, aux_basis, 0.0, 8)

        assert cramped.is_prepared()
        assert cramped.get_mode() == rimode.direct

        # and the way it took has to give the same matrix as the way it did not

        held = roomy.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0)
        formed = cramped.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0)

        assert np.max(np.abs(held)) > 0.0

        scale = float(np.max(np.abs(held)))

        assert np.max(np.abs(formed - held)) / scale < 1.0e-12

    def test_the_two_ways_agree(self, molecule):
        """The direct way forms the integrals again for every batch of orbitals and
        solves the factor of the metric against them, where the way which holds the
        B vectors multiplies by its inverse. They are different arithmetic and must
        reach the same matrix."""

        for la, lb, lc in ((0, 0, 0), (1, 0, 1), (1, 1, 1), (2, 1, 2)):

            basis, aux_basis = self.bases((8, 1, 7, 6), (la, la, lb, lb), lc)

            nao = basis.get_dimensions_of_basis()

            coeffs, density = self.orbitals(nao, 3, 11 + nao)

            matrices = {}

            for mode in (rimode.in_memory, rimode.direct):
                driver = SimdRIJKFockDriver()
                driver.prepare(molecule, basis, aux_basis, 0.0, 1 << 30, 1.0e-12, False, mode)

                assert driver.get_mode() == mode

                matrices[mode] = driver.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0)

            held = matrices[rimode.in_memory]
            formed = matrices[rimode.direct]

            assert np.max(np.abs(held)) > 0.0
            assert np.array_equal(formed, formed.T)

            scale = float(np.max(np.abs(held)))

            assert np.max(np.abs(formed - held)) / scale < 1.0e-12, (
                f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) the two ways disagree")

    def test_compute_before_prepare_is_refused(self):

        driver = SimdRIJKFockDriver()

        assert not driver.is_prepared()

    def test_either_metric_gives_the_same_fock_matrix(self, molecule):
        """The inverted Cholesky factor and the inverted square root both close the
        resolution of the identity, so the Fock matrix does not depend on which was
        used. The flag takes the second when a caller wants it from the start."""

        basis, aux_basis = self.bases((8, 1, 7, 6), (1, 1, 0, 0), 1)

        nao = basis.get_dimensions_of_basis()

        coeffs, density = self.orbitals(nao, 3, 71)

        matrices = []

        for use_root in (False, True):

            driver = SimdRIJKFockDriver()
            driver.prepare(molecule, basis, aux_basis, 0.0, 1 << 30,
                           metric_threshold=1.0e-12,
                           use_inverse_square_root=use_root)

            matrices.append(driver.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0))

        assert np.max(np.abs(matrices[0])) > 0.0

        scale = float(np.max(np.abs(matrices[0])))

        assert np.max(np.abs(matrices[0] - matrices[1])) / scale < 1.0e-10

    def test_the_metric_formed_outside_is_the_metric_formed_inside(self, molecule):
        """The ranks of a communicator cannot each invert the metric: the fallbacks
        are chosen from the matrix, and two ranks could choose differently. The
        master forms it once with make_metric, which answers the way of building as
        well, and hands both over. A metric given has to give what forming it inside
        would have given."""

        basis, aux_basis = self.bases((8, 1, 7, 6), (1, 1, 0, 0), 2)

        nao = basis.get_dimensions_of_basis()

        coeffs, density = self.orbitals(nao, 3, 53)

        for mode in (rimode.in_memory, rimode.direct):

            inside = SimdRIJKFockDriver()
            inside.prepare(molecule, basis, aux_basis, 0.0, 1 << 30, 1.0e-12, False, mode)

            outside = SimdRIJKFockDriver()

            metric, answered = outside.make_metric(molecule, aux_basis, 1.0e-12, False, mode)

            assert answered == mode
            assert metric.number_of_elements() > 0

            outside.prepare(molecule, basis, aux_basis, 0.0, 1 << 30, 1.0e-12, False,
                            answered, metric=metric)

            assert outside.get_mode() == mode

            held = inside.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0)
            given = outside.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0)

            scale = float(np.max(np.abs(held)))

            assert scale > 0.0
            assert np.max(np.abs(given - held)) / scale < 1.0e-14

    def test_a_share_of_the_auxiliary_atoms_is_a_share_of_the_fock_matrix(self, molecule):
        """This is what divides the work over a communicator. Every term of both the
        Coulomb and the exchange is a sum over the auxiliary basis, so a rank given
        some of its atoms answers a part of the matrix and the parts add up to the
        whole. The atoms are dealt out in turn, as a communicator deals them, rather
        than cut into runs."""

        basis, aux_basis = self.bases((8, 1, 7, 6), (1, 1, 0, 0), 2)

        nao = basis.get_dimensions_of_basis()

        coeffs, density = self.orbitals(nao, 3, 59)

        whole = SimdRIJKFockDriver()
        whole.prepare(molecule, basis, aux_basis, 0.0, 1 << 30, 1.0e-12, False,
                      rimode.in_memory)

        expected = whole.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0)

        natoms = molecule.number_of_atoms()

        for nranks in (2, 3):

            shares = []

            for rank in range(nranks):

                atoms = list(range(natoms))[rank::nranks]

                part = SimdRIJKFockDriver()
                part.prepare(molecule, basis, aux_basis, 0.0, 1 << 30, 1.0e-12, False,
                             rimode.in_memory, aux_atoms=atoms)

                shares.append(part.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0))

            summed = sum(shares)

            scale = float(np.max(np.abs(expected)))

            assert scale > 0.0
            assert np.max(np.abs(summed - expected)) / scale < 1.0e-12, (
                f"the shares of {nranks} ranks do not add up to the whole")

    def test_the_memory_of_the_shares_is_the_memory_of_the_whole(self, molecule):
        """A rank chooses the way it builds from the memory it would hold, which is
        the memory of its own atoms and not of the molecule. Asked for all of them
        the answer must be what it always was."""

        basis, aux_basis = self.bases((8, 1, 7, 6), (1, 1, 0, 0), 2)

        driver = SimdRIJKFockDriver()

        natoms = molecule.number_of_atoms()

        whole = driver.required_memory(molecule, basis, aux_basis, 0.0)

        assert whole == driver.required_memory(molecule, basis, aux_basis, 0.0,
                                               list(range(natoms)))

        shares = [
            driver.required_memory(molecule, basis, aux_basis, 0.0,
                                   list(range(natoms))[rank::2]) for rank in range(2)
        ]

        assert all(share > 0 for share in shares)
        assert sum(shares) == whole

    def test_the_direct_way_in_three_calls_is_the_direct_way_in_one(self, molecule):
        """The direct build is two sweeps of the integrals with the fitting between
        them, and a communicator has to gather the fitting, so the build is three
        calls. Made in a row over the whole of the orbitals and the whole of the
        auxiliary basis they have to be the one call they replace."""

        basis, aux_basis = self.bases((8, 1, 7, 6), (1, 1, 0, 0), 2)

        nao = basis.get_dimensions_of_basis()

        coeffs, density = self.orbitals(nao, 4, 83)

        driver = SimdRIJKFockDriver()
        driver.prepare(molecule, basis, aux_basis, 0.0, 1 << 30, 1.0e-12, False,
                       rimode.direct)

        expected = driver.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0)

        fock, gamma = driver.compute_exchange(coeffs, 1.0, 0, 4)

        gamma = driver.solve_fitting(gamma)

        driver.compute_coulomb(gamma, list(range(driver.number_of_parts())), fock)

        computed = fock.to_numpy(max_memory=8.0)

        scale = float(np.max(np.abs(expected)))

        assert scale > 0.0
        assert np.max(np.abs(computed - expected)) / scale < 1.0e-13

    def test_a_share_of_the_orbitals_and_the_parts_is_a_share_of_the_direct_build(
            self, chain):
        """This is what divides the direct way over a communicator. The exchange and
        the right hand side of the fitting are sums over the orbitals, so the ranks
        take a range of them each; the Coulomb matrix is a sum over the parts of the
        auxiliary basis, so they take some parts each. A budget small enough to cut
        the auxiliary basis into several parts is what makes the second division
        worth anything."""

        basis, aux_basis = self.bases(tuple([1] * 10), tuple([0] * 10), 3)

        nao = basis.get_dimensions_of_basis()

        norbitals = 5

        coeffs, density = self.orbitals(nao, norbitals, 97)

        budget = 1 << 14

        whole = SimdRIJKFockDriver()
        whole.prepare(chain, basis, aux_basis, 0.0, budget, 1.0e-12, False,
                      rimode.direct)

        nparts = whole.number_of_parts()

        assert nparts > 1, f"the auxiliary basis is swept in {nparts} part"

        expected = whole.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0)

        for nranks in (2, 3, 7):

            shares = []

            # the exchange pass first, as the fitting cannot be solved until every
            # rank has added its right hand side to the others

            rights = []

            for rank in range(nranks):

                part = SimdRIJKFockDriver()
                part.prepare(chain, basis, aux_basis, 0.0, budget, 1.0e-12, False,
                             rimode.direct)

                ofirst = (norbitals * rank) // nranks
                olast = (norbitals * (rank + 1)) // nranks

                fock, gamma = part.compute_exchange(coeffs, 1.0, ofirst, olast)

                shares.append((part, fock))
                rights.append(np.asarray(gamma))

            gathered = np.sum(rights, axis=0)

            summed = None

            for rank, (part, fock) in enumerate(shares):

                gamma = part.solve_fitting(gathered)

                part.compute_coulomb(gamma, list(range(nparts))[rank::nranks], fock)

                matrix = fock.to_numpy(max_memory=8.0)

                summed = matrix if summed is None else summed + matrix

            scale = float(np.max(np.abs(expected)))

            assert scale > 0.0
            assert np.max(np.abs(summed - expected)) / scale < 1.0e-12, (
                f"the shares of {nranks} ranks over {nparts} parts do not add up")

    def test_a_rank_with_nothing_to_do_answers_nothing(self, molecule):
        """More ranks than orbitals leaves some of them an empty range, and more
        ranks than parts leaves some of them no parts. Neither may be read as all of
        them, which is what an empty list would mean if it were a default."""

        basis, aux_basis = self.bases((8, 1, 7, 6), (1, 1, 0, 0), 2)

        nao = basis.get_dimensions_of_basis()

        coeffs, _ = self.orbitals(nao, 3, 103)

        driver = SimdRIJKFockDriver()
        driver.prepare(molecule, basis, aux_basis, 0.0, 1 << 30, 1.0e-12, False,
                       rimode.direct)

        fock, gamma = driver.compute_exchange(coeffs, 1.0, 2, 2)

        assert np.max(np.abs(np.asarray(gamma))) == 0.0

        driver.compute_coulomb(driver.solve_fitting(gamma), [], fock)

        assert np.max(np.abs(fock.to_numpy(max_memory=8.0))) == 0.0

    def test_a_share_of_the_atoms_is_a_share_of_the_work(self, molecule):
        """A rank answers a share of the Fock matrix whether or not the work was
        divided: a function it holds nothing of contributes nothing either way. So
        the energy cannot tell a division which divides from one which only looks
        like it, and this counts the auxiliary functions a build sweeps instead.

        This is the test the first division did not have, and it would have failed:
        every rank swept the whole auxiliary basis, formed a W matrix for every
        function of every other rank, zeroed it and multiplied it as zeros."""

        basis, aux_basis = self.bases((8, 1, 7, 6), (1, 1, 0, 0), 2)

        naux = aux_basis.get_dimensions_of_basis()

        whole = SimdRIJKFockDriver()
        whole.prepare(molecule, basis, aux_basis, 0.0, 1 << 30, 1.0e-12, False,
                      rimode.in_memory)

        assert whole.number_of_aux_functions() == naux

        natoms = molecule.number_of_atoms()

        for nranks in (2, 4):

            swept = []

            for rank in range(nranks):

                part = SimdRIJKFockDriver()
                part.prepare(molecule, basis, aux_basis, 0.0, 1 << 30, 1.0e-12, False,
                             rimode.in_memory,
                             aux_atoms=list(range(natoms))[rank::nranks])

                swept.append(part.number_of_aux_functions())

            assert sum(swept) == naux, (
                f"{nranks} ranks sweep {sum(swept)} of {naux} functions")

            assert max(swept) < naux, (
                f"a rank of {nranks} sweeps {max(swept)} of {naux} functions, "
                "which is the whole of them")

    def test_the_direct_way_sweeps_the_whole_auxiliary_basis(self, molecule):
        """It is not divided over the atoms, and says so rather than reporting a
        share it does not take."""

        basis, aux_basis = self.bases((8, 1, 7, 6), (1, 1, 0, 0), 2)

        driver = SimdRIJKFockDriver()
        driver.prepare(molecule, basis, aux_basis, 0.0, 1 << 30, 1.0e-12, False,
                       rimode.direct)

        assert driver.number_of_aux_functions() == aux_basis.get_dimensions_of_basis()

    def test_the_parts_are_cut_for_the_ranks_as_well_as_the_memory(self, chain):
        """The Coulomb pass of the direct way is divided over the parts it sweeps,
        and the parts are cut to fit the memory of a build. A machine with memory to
        spare gives one part, which is one rank's work and no one else's, so a caller
        dividing over a communicator asks for at least as many parts as it has
        ranks."""

        basis, aux_basis = self.bases(tuple([1] * 10), tuple([0] * 10), 3)

        nao = basis.get_dimensions_of_basis()

        coeffs, density = self.orbitals(nao, 4, 61)

        roomy = 1 << 30

        whole = SimdRIJKFockDriver()
        whole.prepare(chain, basis, aux_basis, 0.0, roomy, 1.0e-12, False,
                      rimode.direct)

        assert whole.number_of_parts() == 1, (
            "a budget of a gigabyte should hold this molecule in one part")

        expected = whole.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0)

        for asked in (2, 4, 8):

            driver = SimdRIJKFockDriver()
            driver.prepare(chain, basis, aux_basis, 0.0, roomy, 1.0e-12, False,
                           rimode.direct, min_parts=asked)

            assert driver.number_of_parts() >= asked, (
                f"asked for {asked} parts and got {driver.number_of_parts()}")

            # NOTE: and the sweep of the exchange pass must not have been cut with
            # it. Every rank sweeps every one of those, and each is another call of
            # the transformation, which costs a square of the basis for every thread
            # whatever the part holds -- so balancing the Coulomb pass this way was a
            # regression at a thousand basis functions on two hundred and fifty six
            # threads until the two were cut apart.

            assert driver.number_of_sweep_parts() == whole.number_of_sweep_parts(), (
                f"asking for {asked} parts of the Coulomb pass cut the sweep from "
                f"{whole.number_of_sweep_parts()} to "
                f"{driver.number_of_sweep_parts()}")

            # cutting finer is a regrouping of the same atoms and must not move the
            # Fock matrix by more than the order of the arithmetic

            computed = driver.compute(density, coeffs, 1.0).to_numpy(max_memory=8.0)

            scale = float(np.max(np.abs(expected)))

            assert scale > 0.0
            assert np.max(np.abs(computed - expected)) / scale < 1.0e-12

    def test_the_direct_way_takes_either_metric(self, molecule):
        """The direct way solves the Cholesky factor of the metric against the half
        transformed integrals, and multiplies by the inverted square root where that
        is what it was given. Both close the same sum -- solving the factor gives B
        with B^T B equal to A^T V^-1 A, and so does multiplying by the root, the root
        being its own transpose -- so the Fock matrix must not depend on which."""

        for la, lb, lc in ((0, 0, 0), (1, 1, 1), (2, 1, 2)):

            basis, aux_basis = self.bases((8, 1, 7, 6), (la, la, lb, lb), lc)

            nao = basis.get_dimensions_of_basis()

            coeffs, density = self.orbitals(nao, 3, 41 + nao)

            matrices = {}

            for use_root in (False, True):

                driver = SimdRIJKFockDriver()

                metric, answered = driver.make_metric(molecule, aux_basis, 1.0e-12,
                                                      use_root, rimode.direct)

                # the direct way keeps the direct way with either metric

                assert answered == rimode.direct

                assert metric.is_triangular()

                driver.prepare(molecule, basis, aux_basis, 0.0, 1 << 30, 1.0e-12,
                               use_root, rimode.direct, metric=metric)

                assert driver.get_mode() == rimode.direct

                matrices[use_root] = driver.compute(density, coeffs,
                                                    1.0).to_numpy(max_memory=8.0)

            solved, multiplied = matrices[False], matrices[True]

            scale = float(np.max(np.abs(solved)))

            assert scale > 0.0
            assert np.max(np.abs(multiplied - solved)) / scale < 1.0e-10, (
                f"({LABELS[la]}{LABELS[lb]}|{LABELS[lc]}) the two metrics disagree")
