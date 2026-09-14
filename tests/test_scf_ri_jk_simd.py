import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver

# NOTE: the SIMD RI-JK driver is an alternative path through the closed shell Fock
# build, selected by ri_jk_simd. It is the same approximation as the conventional
# RI-JK path, with the same auxiliary basis and the same metric threshold, so the
# two must converge to the same energy. That one number exercises the whole chain,
# the inverted metric, the B vectors, the Y vector, the Coulomb matrix, the W
# matrices and the exchange, through a calculation rather than through constructed
# input.

# NOTE: run under mpirun the same comparisons cover the division of the work over
# the ranks, as each of them then forms a share of every Fock matrix and the shares
# are reduced. The two ways are divided over different things -- the way which holds
# the B vectors over the atoms of the auxiliary basis, the direct way over the
# orbitals and over the parts it sweeps -- so both are worth running there.


class TestScfRiJkSimd:

    @pytest.fixture
    def molecule(self):

        return Molecule.read_str(
            """O  0.000  0.000  0.117
               H  0.000  0.757 -0.467
               H  0.000 -0.757 -0.467""", 'angstrom')

    @pytest.fixture
    def basis(self, molecule):

        return MolecularBasis.read(molecule, 'def2-svp', ostream=None)

    def run_scf(self, molecule, basis, **settings):

        driver = ScfRestrictedDriver(ostream=OutputStream(None))

        driver.conv_thresh = 1.0e-8

        for key, value in settings.items():
            setattr(driver, key, value)

        results = driver.compute(molecule, basis)

        # NOTE: the results are returned on the master alone, and every rank has to
        # reach the comparison, so the energy is handed to all of them.

        energy = (results['scf_energy']
                  if driver.rank == mpi_master() else None)

        return driver.comm.bcast(energy, root=mpi_master())

    @pytest.mark.parametrize('xcfun', [None, 'PBE0', 'B3LYP'])
    def test_simd_matches_the_conventional_driver(self, molecule, basis, xcfun):
        """Hartree-Fock scales the exchange by one and the hybrids by a fraction of
        it, so the functionals cover the scaling as well as the path."""

        settings = {
            'ri_jk': True,
            'ri_auxiliary_basis': 'def2-universal-jkfit',
        }

        if xcfun is not None:
            settings['xcfun'] = xcfun

        conventional = self.run_scf(molecule, basis, **settings)

        simd = self.run_scf(molecule, basis, ri_jk_simd=True, **settings)

        assert abs(simd - conventional) < 1.0e-10, (
            f"{xcfun or 'HF'}: {simd:.12f} against {conventional:.12f}")

    def test_the_approximation_is_the_one_expected(self, molecule, basis):
        """The energy has to differ from the one without the approximation, or the
        test above would pass with the path never taken."""

        plain = self.run_scf(molecule, basis)

        simd = self.run_scf(molecule,
                            basis,
                            ri_jk=True,
                            ri_jk_simd=True,
                            ri_auxiliary_basis='def2-universal-jkfit')

        assert abs(simd - plain) > 1.0e-6
        assert abs(simd - plain) < 1.0e-3

    def test_the_flag_is_off_by_default(self):

        driver = ScfRestrictedDriver(ostream=OutputStream(None))

        assert not driver.ri_jk_simd
        assert driver.ri_memory_budget is None

    def test_the_mode_can_be_chosen(self, molecule, basis):
        """The way the Fock matrices are formed is an input setting, and both ways
        must reach the same energy."""

        settings = {
            'ri_jk': True,
            'ri_jk_simd': True,
            'ri_auxiliary_basis': 'def2-universal-jkfit',
        }

        energies = {}

        for mode in ('in_memory', 'direct'):
            energies[mode] = self.run_scf(molecule, basis, ri_mode=mode, **settings)

        assert abs(energies['in_memory'] - energies['direct']) < 1.0e-10

    def test_a_budget_which_does_not_hold_the_b_vectors_goes_direct(self, molecule,
                                                                    basis):
        """A budget below what the B vectors need selects the way which does not
        hold them, and the energy is the same either way. This is what makes a
        molecule too large to hold reachable."""

        settings = {
            'ri_jk': True,
            'ri_jk_simd': True,
            'ri_auxiliary_basis': 'def2-universal-jkfit',
        }

        roomy = self.run_scf(molecule, basis, ri_memory_budget=1.0, **settings)

        # ten kilobytes, which the B vectors of any molecule exceed

        cramped = self.run_scf(molecule, basis, ri_memory_budget=1.0e-5, **settings)

        assert abs(roomy - cramped) < 1.0e-10

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason="pytest.raises only valid in serial")
    def test_an_unknown_mode_is_refused(self, molecule, basis):
        """A name which is not one of the two ways has to be caught, not quietly
        treated as one of them."""

        from veloxchem.errorhandler import VeloxChemError

        with pytest.raises(VeloxChemError, match='ri_mode'):
            self.run_scf(molecule,
                         basis,
                         ri_jk=True,
                         ri_jk_simd=True,
                         ri_auxiliary_basis='def2-universal-jkfit',
                         ri_mode='sideways')

    def test_the_mode_defaults_to_automatic(self):

        driver = ScfRestrictedDriver(ostream=OutputStream(None))

        assert driver.ri_mode == 'automatic'

    def test_the_metric_route_is_a_setting(self, molecule, basis):
        """The metric may be inverted through its Cholesky factor or through its
        eigenvalues. The conventional RI-JK driver does the second, and convergence
        trouble with the resolution of the identity is known to come from the first
        on a poorly conditioned fitting basis, so the SIMD driver has to be able to
        do it too. Both close the same sum, so a well conditioned metric gives the
        same energy either way."""

        settings = {
            'ri_jk': True,
            'ri_jk_simd': True,
            'ri_auxiliary_basis': 'def2-universal-jkfit',
            'ri_mode': 'in_memory',
        }

        through_cholesky = self.run_scf(molecule, basis,
                                        ri_metric_route='cholesky', **settings)

        through_eigenvalues = self.run_scf(molecule, basis,
                                           ri_metric_route='eigenvalues',
                                           **settings)

        assert abs(through_eigenvalues - through_cholesky) < 1.0e-9, (
            f"{through_eigenvalues:.12f} against {through_cholesky:.12f}")

    def test_the_metric_route_defaults_to_cholesky(self):
        """It is an order of magnitude cheaper, and it is the right default where
        the metric is well conditioned, which is most of the time."""

        driver = ScfRestrictedDriver(ostream=OutputStream(None))

        assert driver.ri_metric_route == 'cholesky'

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason="pytest.raises only valid in serial")
    def test_an_unknown_metric_route_is_refused(self, molecule, basis):

        from veloxchem.errorhandler import VeloxChemError

        with pytest.raises(VeloxChemError, match='ri_metric_route'):
            self.run_scf(molecule,
                         basis,
                         ri_jk=True,
                         ri_jk_simd=True,
                         ri_auxiliary_basis='def2-universal-jkfit',
                         ri_metric_route='sideways')
