import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tddftorbitalresponse import TddftOrbitalResponse
from veloxchem.veloxchemlib import mpi_master


@pytest.mark.solvers
class TestExcitedStateAbsorption:

    def _get_water_and_basis(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, '6-31g', ostream=None)

        return mol, bas

    def _get_formaldehyde_and_basis(self):

        xyz_string = """4
        xyz
        C    0.0000    0.0000    0.0000
        O    0.0000    0.0000    1.2050
        H    0.0000    0.9430   -0.5880
        H    0.0000   -0.9430   -0.5880
        """

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, '6-31g*', ostream=None)

        return mol, bas

    def _get_esa_transition_dipoles(self, mol, bas, nstates, esafrom=None):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = nstates
        lr_drv.esa = True
        if esafrom is not None:
            lr_drv.esa_from_state = esafrom
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            return np.array([
                item['transition_dipole'] for item in lr_results['esa_results']
            ])
        return None

    @staticmethod
    def _match_signs(dipoles, ref_dipoles):

        for i in range(dipoles.shape[0]):
            if np.vdot(dipoles[i], ref_dipoles[i]) < 0.0:
                dipoles[i] *= -1.0

        return dipoles

    def test_esa_transition_density_of_identical_states(self):

        # For identical states the ESA transition density must reduce to the unrelaxed density.

        mol, bas = self._get_water_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 2
        lr_results = lr_drv.compute(mol, bas, scf_results)

        nstates = len(lr_results['eigenvalues'])
        unrelaxed_densities = []
        for s in range(nstates):
            orb_drv_s = TddftOrbitalResponse()
            orb_drv_s.ostream.mute()
            orb_drv_s.state_deriv_index = [s + 1]
            orb_drv_s.compute(mol, bas, scf_drv.scf_tensors, lr_results)
            if lr_drv.rank == mpi_master():
                unrelaxed_densities.append(
                    orb_drv_s.cphf_results['unrelaxed_density_ao'][0])
            else:
                unrelaxed_densities.append(None)

        for s in range(nstates):
            # collective, must be called on all ranks
            eigvec = LinearResponseEigenSolver.get_full_solution_vector(
                lr_results['eigenvectors_distributed'][s])

            if lr_drv.rank == mpi_master():
                nocc = mol.number_of_alpha_occupied_orbitals(bas)
                mo_occ, mo_vir = lr_drv._get_mo_occ_and_mo_vir(
                    scf_results, nocc)
                z_mat, y_mat = lr_drv._get_z_mat_and_y_mat(eigvec, nocc)

                esa_trans_dens = lr_drv._get_esa_transition_density(
                    z_mat, y_mat, z_mat, y_mat, mo_occ, mo_vir)

                assert np.max(np.abs(esa_trans_dens -
                                     unrelaxed_densities[s])) < 1.0e-12

    def test_esa_transition_dipoles_water(self):

        # regression test for the sign of the de-excitation (Y) terms

        mol, bas = self._get_water_and_basis()

        dipoles = self._get_esa_transition_dipoles(mol, bas, 3, esafrom=1)

        if dipoles is not None:
            ref_dipoles = np.array([
                [-0.0794549688, -0.5971935937, 1.6021820000],
                [0.1512160878, -0.0251332634, -0.0018690476],
            ])

            dipoles = self._match_signs(dipoles, ref_dipoles)

            assert dipoles.shape[0] == 2
            assert np.max(np.abs(dipoles - ref_dipoles)) < 1.0e-8

    def test_esa_transition_dipoles_formaldehyde(self):

        # formaldehyde has larger de-excitation amplitudes than water

        mol, bas = self._get_formaldehyde_and_basis()

        dipoles = self._get_esa_transition_dipoles(mol, bas, 4, esafrom=1)

        if dipoles is not None:
            ref_dipoles = np.array([
                [0.0000000000, -0.0599013291, 0.0000000000],
                [0.0000000000, 0.0000000000, 0.0000000000],
                [-0.0477763616, 0.0000000000, 0.0000000000],
            ])

            dipoles = self._match_signs(dipoles, ref_dipoles)

            assert dipoles.shape[0] == 3
            assert np.max(np.abs(dipoles - ref_dipoles)) < 1.0e-6
