from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.sanitychecks import get_ri_auxiliary_basis


@pytest.mark.solvers
class TestRIAuxiliaryBasis:
    """ri_auxiliary_basis is left as set. None, the default, takes
    def2-universal-jfit for RI-J and def2-universal-jkfit for RI-JK."""

    @pytest.fixture
    def molecule(self):

        return Molecule.read_molecule_string(
            """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
            """,
            units='angstrom')

    @pytest.fixture
    def basis(self, molecule):

        return MolecularBasis.read(molecule, 'def2-svp', ostream=None)

    @staticmethod
    def make_driver(tmp_path=None, name='output.txt'):
        """Makes an SCF driver, which writes its output to a file on the master
        where tmp_path is given."""

        if tmp_path is not None and MPI.COMM_WORLD.Get_rank() == mpi_master():
            outfile = tmp_path / name
            driver = ScfRestrictedDriver(ostream=OutputStream(outfile))
        else:
            outfile = None
            driver = ScfRestrictedDriver(ostream=OutputStream(None))

        return driver, outfile

    @staticmethod
    def run_scf(driver, molecule, basis, outfile=None, **settings):
        """Runs the SCF and gives the energy, the auxiliary basis set in the
        results and the output on every rank."""

        driver.conv_thresh = 1.0e-8

        for key, value in settings.items():
            setattr(driver, key, value)

        results = driver.compute(molecule, basis)

        if driver.rank == mpi_master():
            energy = results['scf_energy']
            ri_aux_basis = results['ri_auxiliary_basis']
            text = (outfile.read_text(encoding='utf-8')
                    if outfile is not None else '')
        else:
            energy, ri_aux_basis, text = None, None, ''

        return driver.comm.bcast((energy, ri_aux_basis, text),
                                 root=mpi_master())

    def test_none_takes_the_basis_set_of_the_ri_kind(self):

        assert get_ri_auxiliary_basis(None, 'RI-J') == 'def2-universal-jfit'
        assert get_ri_auxiliary_basis(None, 'RI-JK') == 'def2-universal-jkfit'
        assert get_ri_auxiliary_basis('def2-universal-jkfit',
                                      'RI-J') == 'def2-universal-jkfit'

    def test_the_default_is_none(self):

        driver = ScfRestrictedDriver(ostream=OutputStream(None))

        assert driver.ri_auxiliary_basis is None

    def test_jkfit_is_kept_for_ri_j(self, molecule, basis, tmp_path):
        """A pure functional with def2-universal-jkfit set runs RI-J with it,
        and the output says so; the default runs RI-J with
        def2-universal-jfit."""

        driver, outfile = self.make_driver(tmp_path, 'jkfit.txt')
        jkfit, jkfit_setting, text = self.run_scf(
            driver,
            molecule,
            basis,
            outfile,
            xcfun='blyp',
            ri_auxiliary_basis='def2-universal-jkfit')

        assert driver.ri_auxiliary_basis == 'def2-universal-jkfit'
        assert jkfit_setting == 'def2-universal-jkfit'
        assert 'RI-J uses def2-universal-jkfit as set' in text

        driver, outfile = self.make_driver(tmp_path, 'default.txt')
        default, default_setting, text = self.run_scf(driver,
                                                      molecule,
                                                      basis,
                                                      outfile,
                                                      xcfun='blyp')

        assert default_setting is None
        assert 'RI-J uses def2-universal-jkfit' not in text

        jfit, _, _ = self.run_scf(self.make_driver()[0],
                                  molecule,
                                  basis,
                                  xcfun='blyp',
                                  ri_auxiliary_basis='def2-universal-jfit')

        assert abs(default - jfit) < 1.0e-10
        assert abs(jkfit - jfit) > 1.0e-7

    def test_a_driver_used_again_chooses_again(self, molecule, basis):
        """A driver that ran a hybrid (RI-JK, def2-universal-jkfit) and then a
        pure functional runs RI-J with def2-universal-jfit, as a new driver
        does."""

        driver = self.make_driver()[0]
        self.run_scf(driver, molecule, basis, xcfun='b3lyp')
        again, _, _ = self.run_scf(driver, molecule, basis, xcfun='blyp')

        assert driver.ri_auxiliary_basis is None

        fresh, _, _ = self.run_scf(self.make_driver()[0],
                                   molecule,
                                   basis,
                                   xcfun='blyp')

        assert abs(again - fresh) < 1.0e-10

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason="pytest.raises only valid in serial")
    def test_jfit_is_refused_for_ri_jk(self, molecule, basis):

        from veloxchem.errorhandler import VeloxChemError

        with pytest.raises(VeloxChemError, match='jkfit'):
            self.run_scf(self.make_driver()[0],
                         molecule,
                         basis,
                         xcfun='b3lyp',
                         ri_auxiliary_basis='def2-universal-jfit')
