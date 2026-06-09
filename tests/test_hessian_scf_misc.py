from pathlib import Path
from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import parse_xc_func
from veloxchem.firstorderprop import FirstOrderProperties
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfhessiandriver import ScfHessianDriver
from veloxchem.dftutils import get_default_grid_level
from veloxchem.errorhandler import VeloxChemError


@pytest.mark.solvers
class TestScfHessianDriverMiscellaneous:

    @staticmethod
    def get_h2_molecule_and_basis():

        molecule = Molecule.read_xyz_string("""2
        h2
        H    0.000000    0.000000   -0.350000
        H    0.000000    0.000000    0.350000
        """)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        return molecule, basis

    @staticmethod
    def get_h2o_molecule_and_basis():

        molecule = Molecule.read_xyz_string("""3
        h2o
        O    0.000000    0.000000    0.000000
        H    0.000000   -0.757160    0.586260
        H    0.000000    0.757160    0.586260
        """)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        return molecule, basis

    @staticmethod
    def get_hcl_molecule_and_basis():

        molecule = Molecule.read_xyz_string("""2
        hcl
        H    0.000000    0.000000    0.000000
        Cl   0.000000    0.000000    1.275000
        """)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        return molecule, basis

    @staticmethod
    def get_agcl_molecule_and_basis():

        molecule = Molecule.read_xyz_string("""2
        agcl
        Ag   0.000000    0.000000    0.000000
        Cl   0.000000    2.800000    0.000000
        """)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        return molecule, basis

    @staticmethod
    def get_au2_molecule_and_basis():

        molecule = Molecule.read_xyz_string("""2
        au2
        Au   0.000000    0.000000    0.000000
        Au   0.000000    0.000000    2.800000
        """)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        return molecule, basis

    @staticmethod
    def get_embedded_water_molecule_and_basis():

        molecule = Molecule.read_xyz_string("""3
        xyz
        O    1.2361419   1.0137761  -0.0612424
        H    0.5104418   0.8944555   0.5514190
        H    1.9926927   1.1973129   0.4956931
        """)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        return molecule, basis

    @staticmethod
    def get_fake_molecule_and_basis(elem_id):

        class FakeMolecule:

            @staticmethod
            def get_identifiers():
                return [elem_id]

        class FakeBasis:

            @staticmethod
            def get_number_of_ecp_core_electrons():
                return [0]

            @staticmethod
            def has_ecp():
                return False

        return FakeMolecule(), FakeBasis()

    @staticmethod
    def run_restricted_scf(molecule, basis, xcfun='hf'):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun
        scf_drv.compute(molecule, basis)
        return scf_drv

    def test_update_settings_defaults_cphf_dict(self):

        molecule, basis = self.get_h2_molecule_and_basis()
        scf_drv = self.run_restricted_scf(molecule, basis)

        hess_drv = ScfHessianDriver(scf_drv)
        hess_drv.update_settings({'xcfun': 'hf'}, hess_dict={'numerical': 'no'})

        assert hess_drv.cphf_dict == {}

    def test_compute_runs_internal_scf_when_results_missing(self):

        molecule, basis = self.get_h2o_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()

        hess_drv = ScfHessianDriver(scf_drv)
        hess_drv.ostream.mute()
        hess_drv.compute(molecule, basis)

        assert scf_drv.scf_results is not None
        if scf_drv.rank == mpi_master():
            assert scf_drv.scf_results['scf_type'] == 'restricted'
            assert hess_drv.hessian.shape == (9, 9)
            assert np.max(np.abs(hess_drv.hessian)) > 0.0

    def test_compute_energy_gradient_and_dipole_helpers(self):

        molecule, basis = self.get_h2o_molecule_and_basis()
        scf_drv = self.run_restricted_scf(molecule, basis)
        hess_drv = ScfHessianDriver(scf_drv)

        original_compute = scf_drv.compute
        restart_states = []

        def wrapped_compute(*args, **kwargs):
            restart_states.append(scf_drv.restart)
            return original_compute(*args, **kwargs)

        scf_drv.compute = wrapped_compute

        energy = hess_drv.compute_energy(molecule, basis)

        ref_grad_drv = ScfGradientDriver(scf_drv)
        ref_grad_drv.ostream.mute()
        ref_grad_drv.compute(molecule, basis)

        gradient = hess_drv.compute_gradient(molecule, basis)

        prop = FirstOrderProperties(hess_drv.comm, hess_drv.ostream)
        prop.compute_scf_prop(molecule, basis, scf_drv.scf_results)
        dipole_moment = hess_drv.compute_electric_dipole_moment(molecule, basis)

        hess_drv.numerical = True
        numerical_energy = hess_drv.compute_energy(molecule, basis)

        assert restart_states == [True, True, False]
        assert np.isclose(energy, scf_drv.get_scf_energy())
        assert np.isclose(numerical_energy, scf_drv.get_scf_energy())
        if scf_drv.rank == mpi_master():
            np.testing.assert_allclose(gradient,
                                       ref_grad_drv.gradient,
                                       rtol=1.0e-12,
                                       atol=1.0e-12)
            np.testing.assert_allclose(dipole_moment,
                                       prop.get_property('dipole moment'),
                                       rtol=1.0e-12,
                                       atol=1.0e-12)
        else:
            assert dipole_moment is None

    def test_analytical_hessian_with_point_charges_and_vdw(self):

        molecule, basis = self.get_embedded_water_molecule_and_basis()
        here = Path(__file__).parent
        potfile = str(here / 'data' / 'pe_water.pot')
        vdwfile = str(here / 'data' / 'pe_water.qm_vdw_params.txt')
        reffile = str(here / 'data' /
                      'water_analytical_hessian_pointcharges_scf.txt')

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.point_charges = potfile
        scf_drv.qm_vdw_params = vdwfile
        scf_drv.compute(molecule, basis)

        hess_drv = ScfHessianDriver(scf_drv)
        hess_drv.ostream.mute()
        hess_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            ref_hessian = np.loadtxt(reffile)
            np.testing.assert_allclose(hess_drv.hessian,
                                       ref_hessian,
                                       rtol=1.0e-8,
                                       atol=1.0e-10)

    def test_determine_xc_hessian_grid_level_promotes_supported_cases(self):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        hess_drv = ScfHessianDriver(scf_drv)

        h2o_molecule, h2o_basis = self.get_h2o_molecule_and_basis()
        scf_drv.xcfun = parse_xc_func('PBE')
        assert hess_drv._determine_xc_hessian_grid_level(
            h2o_molecule, h2o_basis, 4) == 6

        scf_drv.xcfun = parse_xc_func('M06')
        assert hess_drv._determine_xc_hessian_grid_level(
            h2o_molecule, h2o_basis, 5) == 6

        hcl_molecule, hcl_basis = self.get_hcl_molecule_and_basis()
        assert hess_drv._determine_xc_hessian_grid_level(
            hcl_molecule, hcl_basis, 5) == 7

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize(
        ('xcfun_label', 'molecule_basis_getter', 'grid_level',
         'expected_message'),
        [('M06', 'get_agcl_molecule_and_basis', 6,
          r'Hessian calculation with M06 functional and effective core '
          r'potential is not supported'),
         ('M06', 'get_au2_molecule_and_basis', 6,
          r'Hessian calculation with M06 functional and effective core '
          r'potential is not supported'),
         ('SCAN', 'get_h2o_molecule_and_basis', 7,
          r'Hessian calculation with SCAN functional and max element id 8 '
          r'is not supported'),
         ('PBE', 'get_fake_molecule_and_basis_37', 4,
          r'Hessian calculation with PBE functional and max element id 37 '
          r'is not supported'),
         ('M06', 'get_fake_molecule_and_basis_19', 5,
          r'Hessian calculation with M06 functional and max element id 19 '
          r'is not supported')])
    def test_determine_xc_hessian_grid_level_rejects_unsupported_cases(
            self, xcfun_label, molecule_basis_getter, grid_level,
            expected_message):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = parse_xc_func(xcfun_label)
        hess_drv = ScfHessianDriver(scf_drv)

        if molecule_basis_getter == 'get_fake_molecule_and_basis_37':
            molecule, basis = self.get_fake_molecule_and_basis(37)
        elif molecule_basis_getter == 'get_fake_molecule_and_basis_19':
            molecule, basis = self.get_fake_molecule_and_basis(19)
        else:
            molecule, basis = getattr(self, molecule_basis_getter)()

        with pytest.raises(VeloxChemError, match=expected_message):
            hess_drv._determine_xc_hessian_grid_level(molecule, basis,
                                                      grid_level)

    def test_determine_xc_hessian_grid_level_with_all_ecp_case(self):

        au2_molecule, au2_basis = self.get_au2_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = parse_xc_func('PBE')

        hess_drv = ScfHessianDriver(scf_drv)
        hess_drv.ostream.mute()
        default_grid_level = get_default_grid_level(scf_drv.xcfun)

        assert hess_drv._determine_xc_hessian_grid_level(
            au2_molecule, au2_basis, default_grid_level) == 6
