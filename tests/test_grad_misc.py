from copy import deepcopy

from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.dispersionmodel import DispersionModel
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver


@pytest.mark.solvers
class TestScfGradientDriverCoverage:

    @staticmethod
    def get_h2_molecule_and_basis():

        xyz_string = """
        2
        h2
        H    0.000000    0.000000   -0.350000
        H    0.000000    0.000000    0.350000
        """
        molecule = Molecule.read_xyz_string(xyz_string)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        return molecule, basis

    @staticmethod
    def get_h2o_molecule_and_basis():

        xyz_string = """
        3
        h2o
        O    0.000000    0.000000    0.000000
        H    0.000000   -0.757160    0.586260
        H    0.000000    0.757160    0.586260
        """
        molecule = Molecule.read_xyz_string(xyz_string)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        return molecule, basis

    @staticmethod
    def get_ch3_molecule_and_basis():

        xyz_string = """
        4
        ch3
        C   -1.85334300   -0.63945100    1.29623300
        H   -2.40884500   -1.56570200    1.04276400
        H   -2.24160900   -0.22442700    2.24900500
        H   -1.98830700    0.10613700    0.48589200
        """
        molecule = Molecule.read_xyz_string(xyz_string)
        molecule.set_multiplicity(2)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        return molecule, basis

    @staticmethod
    def run_restricted_scf(molecule, basis, xcfun='hf'):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun
        scf_results = scf_drv.compute(molecule, basis)
        return scf_drv, scf_results

    @staticmethod
    def run_unrestricted_scf(molecule, basis, xcfun='hf'):

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun
        scf_results = scf_drv.compute(molecule, basis)
        return scf_drv, scf_results

    def test_compute_uses_internal_scf_results_for_numerical_gradient(
            self, monkeypatch):

        molecule, basis = self.get_h2_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_scf(molecule, basis)
        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.numerical = True

        seen = {}

        def fake_compute_numerical(self, mol, ao_basis, scf_results):
            seen['molecule'] = mol
            seen['basis'] = ao_basis
            seen['scf_results'] = scf_results
            self.gradient = np.zeros((mol.number_of_atoms(), 3))

        monkeypatch.setattr(ScfGradientDriver, 'compute_numerical',
                            fake_compute_numerical)

        grad_drv.compute(molecule, basis)

        assert seen['molecule'] is molecule
        assert seen['basis'] is basis
        assert seen['scf_results'] is scf_drv.scf_results
        assert np.allclose(grad_drv.get_gradient(), 0.0)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_restricted_openshell_marker(self):

        molecule, basis = self.get_h2_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_scf(molecule, basis)
        grad_drv = ScfGradientDriver(scf_drv)

        invalid_results = dict(scf_results)
        invalid_results['scf_type'] = 'restricted_openshell'

        with pytest.raises(AssertionError,
                           match='Not implemented for restricted open-shell'):
            grad_drv.compute(molecule, basis, invalid_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_restricted_ri_gradient_rejects_non_def2_basis(self):

        molecule, basis = self.get_h2o_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_scf(molecule, basis)
        scf_drv.ri_coulomb = True

        grad_drv = ScfGradientDriver(scf_drv)

        with pytest.raises(AssertionError, match='Invalid basis set for RI-J'):
            grad_drv.compute(molecule, basis, scf_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_unrestricted_ri_gradient_rejects_non_def2_basis(self):

        molecule, basis = self.get_ch3_molecule_and_basis()
        scf_drv, scf_results = self.run_unrestricted_scf(molecule, basis)
        scf_drv.ri_coulomb = True

        grad_drv = ScfGradientDriver(scf_drv)

        with pytest.raises(AssertionError, match='Invalid basis set for RI-J'):
            grad_drv.compute(molecule, basis, scf_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_restricted_ri_jk_gradient_rejects(self):

        molecule, basis = self.get_h2o_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_scf(molecule, basis)
        scf_drv.ri_jk = True

        grad_drv = ScfGradientDriver(scf_drv)

        with pytest.raises(AssertionError, match='RI-JK is not yet supported'):
            grad_drv.compute(molecule, basis, scf_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_unrestricted_ri_jk_gradient_rejects(self):

        molecule, basis = self.get_ch3_molecule_and_basis()
        scf_drv, scf_results = self.run_unrestricted_scf(molecule, basis)
        scf_drv.ri_jk = True

        grad_drv = ScfGradientDriver(scf_drv)

        with pytest.raises(AssertionError, match='RI-JK is not yet supported'):
            grad_drv.compute(molecule, basis, scf_results)

    def test_restricted_gradient_timing_output_path(self):

        molecule, basis = self.get_h2_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_scf(molecule, basis)
        scf_drv.timing = True

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis, scf_results)

        assert np.max(np.abs(grad_drv.get_gradient())) > 0.0

    def test_compute_energy_updates_results_and_restart_mode(self):

        molecule, basis = self.get_h2_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_scf(molecule, basis)
        grad_drv = ScfGradientDriver(scf_drv)

        original_compute = scf_drv.compute
        restart_states = []

        def wrapped_compute(*args, **kwargs):
            restart_states.append(scf_drv.restart)
            return original_compute(*args, **kwargs)

        scf_drv.compute = wrapped_compute

        updated_results = {'sentinel': 'keep'}
        energy = grad_drv.compute_energy(molecule, basis, updated_results)

        assert restart_states == [True]
        assert np.isclose(energy, scf_drv.get_scf_energy())
        if grad_drv.rank == mpi_master():
            assert 'scf_energy' in updated_results
            assert updated_results['sentinel'] == 'keep'

        grad_drv.numerical = True
        grad_drv.compute_energy(molecule, basis)

        assert restart_states == [True, False]

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_energy_raises_on_unconverged_scf(self, monkeypatch):

        molecule, basis = self.get_h2_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_scf(molecule, basis)
        grad_drv = ScfGradientDriver(scf_drv)

        def fake_compute(*args, **kwargs):
            scf_drv._is_converged = False
            return {}

        monkeypatch.setattr(scf_drv, 'compute', fake_compute)

        with pytest.raises(AssertionError, match='SCF did not converge'):
            grad_drv.compute_energy(molecule, basis)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_unrestricted_cpcm_rejects_smd_in_gradient(self):

        molecule, basis = self.get_ch3_molecule_and_basis()
        scf_drv, scf_results = self.run_unrestricted_scf(
            molecule, basis, 'b3lyp')
        scf_drv._cpcm = True
        scf_drv._smd = True

        grad_drv = ScfGradientDriver(scf_drv)

        with pytest.raises(AssertionError,
                           match='Cannot use SMD in gradient calculation'):
            grad_drv.compute_analytical_unrestricted(molecule, basis,
                                                     scf_results)

    def test_unrestricted_dispersion_gradient_matches_added_d4_component(self):

        molecule, basis = self.get_ch3_molecule_and_basis()

        base_scf_drv, base_results = self.run_unrestricted_scf(
            molecule, basis, 'b3lyp')
        base_grad_drv = ScfGradientDriver(base_scf_drv)
        base_grad_drv.compute(molecule, basis, base_results)
        base_gradient = base_grad_drv.get_gradient()

        disp_scf_drv, disp_results = self.run_unrestricted_scf(
            molecule, basis, 'b3lyp')
        disp_scf_drv.dispersion = True
        disp_grad_drv = ScfGradientDriver(disp_scf_drv)
        disp_grad_drv.compute(molecule, basis, disp_results)
        disp_gradient = disp_grad_drv.get_gradient()

        d4_model = DispersionModel()
        d4_model.compute(molecule, 'b3lyp')

        assert np.max(
            np.abs((disp_gradient - base_gradient) -
                   d4_model.get_gradient())) < 1.0e-8
