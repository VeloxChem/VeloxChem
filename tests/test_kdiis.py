import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.errorhandler import VeloxChemError
from veloxchem.kdiis import KDiis
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.veloxchemlib import mpi_master


class TestKDiis:

    @staticmethod
    def get_water_and_basis(charge=0, multiplicity=1):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """

        molecule = Molecule.read_xyz_string(xyz_string)
        molecule.set_charge(charge)
        molecule.set_multiplicity(multiplicity)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

        return molecule, basis

    def test_restricted_first_order_update(self):

        fock = np.array([[-1.0, 0.0, 0.10], [0.0, -0.5, 0.0], [0.10, 0.0, 0.5]])
        coefficients = np.eye(3)
        updater = KDiis(max_vectors=3)

        new_coefficients, _ = updater.update_restricted(
            fock, coefficients, 2, np.eye(3), np.eye(3))
        new_coefficients = new_coefficients[0]
        new_fock_mo = np.linalg.multi_dot(
            [new_coefficients.T, fock, new_coefficients])

        assert np.allclose(new_coefficients.T @ new_coefficients, np.eye(3))
        assert abs(new_fock_mo[2, 0]) < abs(fock[2, 0])
        assert len(updater.error_vectors) == 1
        assert len(updater.fock_matrices) == 1

    def test_history_is_bounded(self):

        updater = KDiis(max_vectors=2)
        coefficients = np.eye(2)

        for coupling in [0.10, 0.08, 0.06]:
            fock = np.array([[-1.0, coupling], [coupling, 0.5]])
            updater.update_restricted(fock, coefficients, 1, np.eye(2),
                                      np.eye(2))

        assert len(updater.error_vectors) == 2
        assert len(updater.fock_matrices) == 2

    def test_error_vector_is_invariant_to_mo_gauge(self):

        fock = np.array([[-1.2, 0.1, 0.2, -0.1],
                         [0.1, -0.7, 0.05, 0.15],
                         [0.2, 0.05, 0.4, 0.03],
                         [-0.1, 0.15, 0.03, 0.9]])
        occ_angle = 0.37
        vir_angle = -0.61
        occ_rotation = np.array([[np.cos(occ_angle), -np.sin(occ_angle)],
                                 [np.sin(occ_angle), np.cos(occ_angle)]])
        vir_rotation = np.array([[np.cos(vir_angle), -np.sin(vir_angle)],
                                 [np.sin(vir_angle), np.cos(vir_angle)]])
        gauge_rotation = np.zeros((4, 4))
        gauge_rotation[:2, :2] = occ_rotation
        gauge_rotation[2:, 2:] = vir_rotation

        reference = KDiis()
        rotated = KDiis()
        reference.update_restricted(fock, np.eye(4), 2, np.eye(4), np.eye(4))
        rotated.update_restricted(fock, gauge_rotation, 2, np.eye(4),
                                  np.eye(4))

        assert np.allclose(reference.error_vectors[0],
                           rotated.error_vectors[0])

    @pytest.mark.solvers
    @pytest.mark.parametrize('xcfun, reference_energy', [
        ('hf', -74.9629282848),
        ('slda', -74.9283891389),
    ])
    def test_restricted_scf(self, xcfun, reference_energy):

        molecule, basis = self.get_water_and_basis()
        driver = ScfRestrictedDriver()
        driver.ostream.mute()
        driver.acc_type = 'kdiis'
        driver.max_iter = 80
        driver.xcfun = xcfun
        driver.grid_level = 1
        results = driver.compute(molecule, basis)

        assert driver.is_converged
        if driver.rank == mpi_master():
            assert abs(results['scf_energy'] - reference_energy) < 1.0e-8

    @pytest.mark.solvers
    @pytest.mark.parametrize('xcfun, reference_energy', [
        ('hf', -74.6557065987),
        ('slda', -74.5498369093),
    ])
    def test_unrestricted_scf(self, xcfun, reference_energy):

        molecule, basis = self.get_water_and_basis(charge=1, multiplicity=2)
        driver = ScfUnrestrictedDriver()
        driver.ostream.mute()
        driver.acc_type = 'kdiis'
        driver.max_iter = 80
        driver.xcfun = xcfun
        driver.grid_level = 1
        results = driver.compute(molecule, basis)

        assert driver.is_converged
        if driver.rank == mpi_master():
            assert abs(results['scf_energy'] - reference_energy) < 1.0e-8

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize('modifier', ['pfon', 'density_damping', 'mom'])
    def test_incompatible_scf_modifiers(self, modifier):

        molecule, basis = self.get_water_and_basis()
        driver = ScfRestrictedDriver()
        driver.ostream.mute()
        driver.acc_type = 'kdiis'
        if modifier == 'mom':
            driver._mom = (None, None)
        else:
            setattr(driver, modifier, True)

        with pytest.raises(VeloxChemError, match='KDIIS is incompatible'):
            driver.compute(molecule, basis)

    @pytest.mark.solvers
    def test_level_shift_is_applied_to_denominators(self):

        molecule, basis = self.get_water_and_basis()
        driver = ScfRestrictedDriver()
        driver.ostream.mute()
        driver.acc_type = 'kdiis'
        driver.level_shifting = 0.2
        driver.max_iter = 80
        results = driver.compute(molecule, basis)

        assert driver.is_converged
        assert driver.level_shifting == 0.0
        if driver.rank == mpi_master():
            assert abs(results['scf_energy'] - (-74.9629282848)) < 1.0e-8

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_restricted_open_shell_is_rejected(self):

        molecule, basis = self.get_water_and_basis(charge=1, multiplicity=2)
        driver = ScfRestrictedOpenDriver()
        driver.ostream.mute()
        driver.acc_type = 'kdiis'

        with pytest.raises(
                VeloxChemError,
                match='available only for restricted and unrestricted'):
            driver.compute(molecule, basis)
