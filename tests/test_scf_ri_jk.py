import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver

# NOTE: the SIMD RI-JK driver is the path through the closed shell Fock
# build, selected by ri. It is the same approximation as the conventional
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


class TestScfRiJk:

    @pytest.fixture
    def molecule(self):

        return Molecule.read_str(
            """O  0.000  0.000  0.117
               H  0.000  0.757 -0.467
               H  0.000 -0.757 -0.467""", 'angstrom')

    @pytest.fixture
    def basis(self, molecule):

        return MolecularBasis.read(molecule, 'def2-svp', ostream=None)

    def run_scf(self, molecule, basis, ri=True, **settings):

        driver = ScfRestrictedDriver(ostream=OutputStream(None))
        driver.ri = ri

        driver.conv_thresh = 1.0e-8

        for key, value in settings.items():
            setattr(driver, key, value)

        results = driver.compute(molecule, basis)

        # NOTE: the results are returned on the master alone, and every rank has to
        # reach the comparison, so the energy is handed to all of them.

        energy = (results['scf_energy']
                  if driver.rank == mpi_master() else None)

        return driver.comm.bcast(energy, root=mpi_master())

    @pytest.mark.parametrize('xcfun, ref_energy', [
        (None, -75.961015181386),
        ('PBE0', -76.276223087084),
        ('B3LYP', -76.358082329970),
    ])
    def test_simd_matches_the_reference_energy(self, molecule, basis, xcfun,
                                               ref_energy):
        """Hartree-Fock scales the exchange by one and the hybrids by a fraction of
        it, so the functionals cover the scaling as well as the path. The reference
        energies are the conventional RI-JK energies, which the SIMD path has to
        reproduce."""

        settings = {
            'ri_auxiliary_basis': 'def2-universal-jkfit',
        }

        if xcfun is not None:
            settings['xcfun'] = xcfun

        energy = self.run_scf(molecule, basis, **settings)

        assert abs(energy - ref_energy) < 1.0e-9, (
            f"{xcfun or 'HF'}: {energy:.12f} against {ref_energy:.12f}")

    def test_the_approximation_is_the_one_expected(self, molecule, basis):
        """The energy has to differ from the one without the approximation, or the
        test above would pass with the path never taken."""

        plain = self.run_scf(molecule, basis, ri=False)

        simd = self.run_scf(molecule,
                            basis,
                            ri_auxiliary_basis='def2-universal-jkfit')

        assert abs(simd - plain) > 1.0e-6
        assert abs(simd - plain) < 1.0e-3

    def test_ri_is_enabled_by_default(self):

        driver = ScfRestrictedDriver(ostream=OutputStream(None))

        assert driver.ri
        assert driver.ri_memory_budget is None

    @pytest.mark.parametrize('xcfun', [None, 'BLYP'])
    def test_the_mode_can_be_chosen(self, molecule, basis, xcfun):
        """The way the Fock matrices are formed is an input setting, and both ways
        must reach the same energy. A pure functional takes RI-J and the others
        RI-JK, and each has both ways."""

        settings = {
            'ri_auxiliary_basis': 'def2-universal-jkfit',
        }

        if xcfun is not None:
            settings['xcfun'] = xcfun

        energies = {}

        for mode in ('in_memory', 'direct'):
            energies[mode] = self.run_scf(molecule, basis, ri_mode=mode, **settings)

        assert abs(energies['in_memory'] - energies['direct']) < 1.0e-10

    @pytest.mark.parametrize('unrestricted, xcfun', [
        (False, None),
        (False, 'CAM-B3LYP'),
        (True, None),
        (False, 'BLYP'),
        (True, 'BLYP'),
    ])
    def test_a_budget_below_the_peak_falls_back_to_four_center(
            self, molecule, basis, unrestricted, xcfun):
        """A budget below the peak RI is predicted to reach has the Fock
        matrices formed from the four-center integrals, which gives the energy
        without the approximation rather than the direct way. The settings which
        asked for the approximation are left as they were set. The closed shell,
        the range separated and the open shell builds each take the fallback
        apart, and a pure functional takes it through RI-J rather than RI-JK, so
        each is run."""

        if unrestricted:
            molecule.set_charge(1)
            molecule.set_multiplicity(2)

        def run(ri, budget=None):
            driver_type = (ScfUnrestrictedDriver
                           if unrestricted else ScfRestrictedDriver)
            driver = driver_type(ostream=OutputStream(None))
            driver.ri = ri
            driver.conv_thresh = 1.0e-8
            driver.ri_auxiliary_basis = 'def2-universal-jkfit'
            driver.ri_memory_budget = budget
            if xcfun is not None:
                driver.xcfun = xcfun
            results = driver.compute(molecule, basis)
            energy = (results['scf_energy']
                      if driver.rank == mpi_master() else None)
            return driver, driver.comm.bcast(energy, root=mpi_master())

        _, plain = run(ri=False)

        # ten kilobytes, which the peak of any molecule exceeds

        driver, cramped = run(ri=True, budget=1.0e-5)

        assert abs(cramped - plain) < 1.0e-9, (
            f"{cramped:.12f} against {plain:.12f}")

        pure = (xcfun == 'BLYP')

        assert driver.ri
        assert driver.ri_jk == (not pure)
        assert driver.ri_coulomb == pure
        assert driver.ri_mode == 'in_memory'
        assert driver._ri_drv is None

        _, roomy = run(ri=True, budget=1.0)

        assert abs(roomy - plain) > 1.0e-6

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason="pytest.raises only valid in serial")
    @pytest.mark.parametrize('mode', ['sideways', 'automatic'])
    def test_an_unknown_mode_is_refused(self, molecule, basis, mode):
        """A name which is not one of the two ways has to be caught, not quietly
        treated as one of them. The way is no longer chosen from the memory, so
        automatic is not one of them either."""

        from veloxchem.errorhandler import VeloxChemError

        with pytest.raises(VeloxChemError, match='ri_mode'):
            self.run_scf(molecule,
                         basis,
                         ri_auxiliary_basis='def2-universal-jkfit',
                         ri_mode=mode)

    def test_the_mode_defaults_to_in_memory(self):

        driver = ScfRestrictedDriver(ostream=OutputStream(None))

        assert driver.ri_mode == 'in_memory'

    def test_the_metric_route_is_a_setting(self, molecule, basis):
        """The metric may be inverted through its Cholesky factor or through its
        eigenvalues. The conventional RI-JK driver does the second, and convergence
        trouble with the resolution of the identity is known to come from the first
        on a poorly conditioned fitting basis, so the SIMD driver has to be able to
        do it too. Both close the same sum, so a well conditioned metric gives the
        same energy either way."""

        settings = {
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
        driver.ri = False

        assert driver.ri_metric_route == 'cholesky'

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason="pytest.raises only valid in serial")
    def test_an_unknown_metric_route_is_refused(self, molecule, basis):

        from veloxchem.errorhandler import VeloxChemError

        with pytest.raises(VeloxChemError, match='ri_metric_route'):
            self.run_scf(molecule,
                         basis,
                         ri_auxiliary_basis='def2-universal-jkfit',
                         ri_metric_route='sideways')
