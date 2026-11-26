import numpy as np
import h5py
from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.distributedarray import DistributedArray
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tddftorbitalresponse import TddftOrbitalResponse
from veloxchem.errorhandler import VeloxChemError


@pytest.mark.solvers
class TestRPA:

    def run_rpa(self,
                xcfun_label,
                basis_label,
                ri,
                ref_exc_enes,
                ref_osc_str,
                tol,
                max_subspace_dim=None):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)

        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 5
        lr_drv.max_subspace_dim = max_subspace_dim
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_openshell_molecule(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(3)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        # empty scf results just for testing
        scf_results = {}

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()

        with pytest.raises(
                VeloxChemError,
                match="LinearResponseEigenSolver: not implemented for unrestricted case"):
            lr_results_not_used = lr_drv.compute(mol, bas, scf_results)

    @pytest.mark.parametrize('ri', [False, True])
    def test_hf_svp(self, ri):

        # vlxtag: RHF, Absorption, TDHF

        ref_exc_enes = {
            False: np.array(
                [0.33973039, 0.40464346, 0.43325535, 0.49805052, 0.55325390]),
            True: np.array(
                [0.33973064, 0.40463219, 0.43324153, 0.49804893, 0.55325612])
        }[ri]

        ref_osc_str = np.array(
            [0.023665, 0.000000, 0.097765, 0.086454, 0.291919])

        self.run_rpa('hf',
                     'def2-svp',
                     ri,
                     ref_exc_enes,
                     ref_osc_str,
                     1.0e-6,
                     max_subspace_dim=50)

    @pytest.mark.parametrize('ri', [False, True])
    def test_slda_svp(self, ri):

        # vlxtag: RKS, Absorption, TDDFT

        ref_exc_enes = {
            False: np.array(
                [0.27459103, 0.34805152, 0.35179468, 0.43068396, 0.51123472]),
            True: np.array(
                [0.27451655, 0.34787592, 0.35172595, 0.43054087, 0.51122210])
        }[ri]

        ref_osc_str = {
            False: np.array([0.018005, 0.000000, 0.075418, 0.059909, 0.258165]),
            True: np.array([0.017964, 0.000000, 0.075335, 0.059940, 0.258031])
        }[ri]

        self.run_rpa('slda', 'def2-svp', ri, ref_exc_enes, ref_osc_str, 1.0e-5)

    @pytest.mark.parametrize('ri', [False, True])
    def test_b3lyp_svp(self, ri):

        # vlxtag: RKS, Absorption, TDDFT

        ref_exc_enes = {
            False: np.array(
                [0.27940792, 0.35031697, 0.36276301, 0.43718376, 0.51430795]),
            True: np.array(
                [0.27940590, 0.35030507, 0.36275840, 0.43717334, 0.51430704])
        }[ri]

        ref_osc_str = {
            False: np.array([0.018243, 0.000000, 0.078326, 0.061329, 0.272933]),
            True: np.array([0.018244, 0.000000, 0.078323, 0.061347, 0.272881])
        }[ri]

        self.run_rpa('b3lyp', 'def2-svp', ri, ref_exc_enes, ref_osc_str, 1.0e-5)

    @pytest.mark.parametrize('ri', [False, True])
    def test_camb3lyp_svp(self, ri):

        # vlxtag: RKS, Absorption, TDDFT

        ref_exc_enes = {
            False: np.array(
                [0.28249530, 0.35533592, 0.36664012, 0.44312170, 0.51548914]),
            True: np.array(
                [0.28249339, 0.35532424, 0.36663576, 0.44311156, 0.51548814])
        }[ri]

        ref_osc_str = {
            False: np.array([0.018219, 0.000000, 0.077391, 0.059649, 0.273690]),
            True: np.array([0.018220, 0.000000, 0.077399, 0.059668, 0.273638])
        }[ri]

        self.run_rpa('cam-b3lyp', 'def2-svp', ri, ref_exc_enes, ref_osc_str,
                     1.0e-5)

    @pytest.mark.parametrize('ri', [False, True])
    def test_camb3lyp_tzvp(self, ri):

        # vlxtag: RKS, Absorption, TDDFT

        ref_exc_enes = {
            False: np.array(
                [0.27984870, 0.35148791, 0.36310452, 0.43770042, 0.50331834]),
            True: np.array(
                [0.27985380, 0.35148692, 0.36311413, 0.43773206, 0.50333810])
        }[ri]

        ref_osc_str = np.array(
            [0.033878, 0.000000, 0.099329, 0.053065, 0.229533])

        self.run_rpa('cam-b3lyp', 'def2-tzvp', ri, ref_exc_enes, ref_osc_str,
                     1.0e-5)

    @pytest.mark.parametrize('ri', [False, True])
    def test_tpssh_svp(self, ri):

        # vlxtag: RKS, Absorption, TDDFT

        ref_exc_enes = {
            False: np.array(
                [0.28898749, 0.36043436, 0.37287451, 0.44801126, 0.52385373]),
            True: np.array(
                [0.28898571, 0.36042311, 0.37287141, 0.44799951, 0.52385237])
        }[ri]

        ref_osc_str = {
            False: np.array([0.020158, 0.000000, 0.085796, 0.074117, 0.272050]),
            True: np.array([0.020160, 0.000000, 0.085799, 0.074135, 0.271989])
        }[ri]

        self.run_rpa('tpssh', 'def2-svp', ri, ref_exc_enes, ref_osc_str, 1.0e-5)

    @pytest.mark.parametrize('ri', [False, True])
    def test_detach_attach_charges_water_631g(self, ri):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, '6-31g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 3
        lr_drv.detach_attach = True
        lr_drv.detach_attach_charges = True
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            detachment = lr_results['detachment_charges']
            attachment = lr_results['attachment_charges']
            ref_detachment = {
                False: np.array([
                    [-1.0029851374, -0.0000948340, -0.0000948341],
                    [-1.0027814588, -0.0002540890, -0.0002540890],
                    [-0.9143468349, -0.0447074693, -0.0447074695],
                ]),
                True: np.array([[-1.0029853675, -0.0000948707, -0.0000948708],
                                [-1.0027818697, -0.0002541342, -0.0002541342],
                                [-0.9143439234, -0.0447087981, -0.0447087983]])
            }[ri]
            ref_attachment = {
                False: np.array([
                    [0.2897746805, 0.3567000207, 0.3567001044],
                    [0.3056100033, 0.3488398628, 0.3488397706],
                    [0.2784322949, 0.3626647099, 0.3626647689],
                ]),
                True: np.array([[0.2897635053, 0.3567057600, 0.3567058437],
                                [0.3055885369, 0.3488508467, 0.3488507544],
                                [0.2784174710, 0.3626719950, 0.3626720539]])
            }[ri]

            assert detachment.shape == (3, mol.number_of_atoms())
            assert attachment.shape == (3, mol.number_of_atoms())
            assert np.max(np.abs(detachment - ref_detachment)) < 1.0e-6
            assert np.max(np.abs(attachment - ref_attachment)) < 1.0e-6

    @pytest.mark.parametrize('ri', [False, True])
    def test_detach_attach_density_consistency(self, ri):

        # detach + attach must reproduce the unrelaxed difference density

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, '6-31g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 1
        lr_results = lr_drv.compute(mol, bas, scf_results)

        orb_drv = TddftOrbitalResponse()
        orb_drv.ostream.mute()
        orb_drv.state_deriv_index = [1]
        orb_drv.compute(mol, bas, scf_drv.scf_tensors, lr_results)

        eigvec = lr_drv.get_full_solution_vector(
            lr_results['eigenvectors_distributed'][0])

        if lr_drv.rank == mpi_master():
            nocc = mol.number_of_alpha_occupied_orbitals(bas)
            mo_occ, mo_vir = lr_drv._get_mo_occ_and_mo_vir(scf_results, nocc)
            z_mat, y_mat = lr_drv._get_z_mat_and_y_mat(eigvec, nocc)

            dens_D, dens_A = lr_drv.get_detach_attach_densities(
                z_mat, y_mat, mo_occ, mo_vir)
            unrel_dens = orb_drv.cphf_results['unrelaxed_density_ao'][0]

            assert np.max(np.abs(dens_D + dens_A - unrel_dens)) < 1.0e-12

            dens_D_mo = np.linalg.multi_dot([mo_occ.T, dens_D, mo_occ])
            dens_A_mo = np.linalg.multi_dot([mo_vir.T, dens_A, mo_vir])
            assert np.max(np.linalg.eigvalsh(dens_D_mo)) < 1.0e-12
            assert np.min(np.linalg.eigvalsh(dens_A_mo)) > -1.0e-12

    @pytest.mark.parametrize('ri', [False, True])
    def test_esa_water_631g(self, ri):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, '6-31g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 3
        lr_drv.esa = True
        lr_drv.esa_from_state = 1
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            esa_results = lr_results['esa_results']

            ref_excitation_energies = {
                False: np.array([0.0706309931, 0.0887848965]),
                True: np.array([0.0706208215, 0.0887779539]),
            }[ri]
            ref_oscillator_strengths = {
                False: np.array([0.1379630687, 0.0013910509]),
                True: np.array([0.1379461067, 0.0013909689]),
            }[ri]
            ref_transition_dipoles = {
                False: np.array([
                    [-0.0794549688, -0.5971935937, 1.6021820000],
                    [0.1512160878, -0.0251332634, -0.0018690476],
                ]),
                True: np.array([
                    [-0.0794558058, -0.5971998843, 1.6021988767],
                    [0.1512175387, -0.0251335045, -0.0018690655],
                ]),
            }[ri]

            src_states = [item['from_state'] for item in esa_results]
            dest_states = [item['to_state'] for item in esa_results]

            excitation_energies = np.array(
                [item['excitation_energy'] for item in esa_results])
            oscillator_strengths = np.array(
                [item['oscillator_strength'] for item in esa_results])
            transition_dipoles = np.array(
                [item['transition_dipole'] for item in esa_results])

            for i in range(transition_dipoles.shape[0]):
                if np.vdot(transition_dipoles[i],
                           ref_transition_dipoles[i]) < 0.0:
                    transition_dipoles[i] *= -1.0

            assert len(esa_results) == 2
            assert src_states == ['S1', 'S1']
            assert dest_states == ['S2', 'S3']
            assert np.max(np.abs(excitation_energies -
                                 ref_excitation_energies)) < 1.0e-8
            assert np.max(
                np.abs(oscillator_strengths -
                       ref_oscillator_strengths)) < 1.0e-5
            assert np.max(np.abs(transition_dipoles -
                                 ref_transition_dipoles)) < 1.0e-5

    def run_rpa_with_ecp(self, ref_exc_enes, ref_osc_str, tol, ri):

        xyz_string = """2
        xyz
        Au 0 0 0
        Au 0 0 2.88
        """
        mol = Molecule.read_xyz_string(xyz_string)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 3
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    @pytest.mark.parametrize('ri', [
        pytest.param(False, marks=pytest.mark.timeconsuming),
        True,
    ])
    def test_hf_with_ecp(self, ri):

        # vlxtag: RHF, Absorption, TDHF

        ref_exc_enes = {
            False: np.array([0.10859702, 0.15291543, 0.16111886]),
            True: np.array([0.10859685, 0.15290990, 0.16111267])
        }[ri]
        ref_osc_str = np.array([0.4439, 0.0000, 0.0025])

        self.run_rpa_with_ecp(ref_exc_enes, ref_osc_str, 1.0e-6, ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_checkpoint_and_restart(self, tmp_path, ri):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()

        filename = str(tmp_path / "water_restart")
        # To avoid inconsistency across MPI ranks
        filename = scf_drv.comm.bcast(filename, root=mpi_master())

        scf_drv.filename = filename
        scf_drv.compute(mol, bas)
        scf_drv.restart = True
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        lr_drv = LinearResponseEigenSolver()
        lr_drv.filename = filename
        lr_drv.ostream.mute()

        lr_drv.nstates = 5
        lr_results_not_used = lr_drv.compute(mol, bas, scf_results)

        lr_drv.restart = True
        lr_drv.nstates = 10
        lr_results_first = lr_drv.compute(mol, bas, scf_results)
        assert lr_drv.restart is True

        lr_drv.restart = False
        lr_drv.nstates = 10
        lr_results_second = lr_drv.compute(mol, bas, scf_results)
        assert lr_drv.restart is False

        if scf_drv.rank == mpi_master():
            keys = ['eigenvalues', 'oscillator_strengths']
            tols = [1e-8, 1e-5]
            for key, tol in zip(keys, tols):
                assert np.max(
                    np.abs(lr_results_first[key] -
                           lr_results_second[key])) < tol

    @pytest.mark.parametrize('ri', [False, True])
    def test_restart_allows_missing_nstates(self, tmp_path, ri):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()

        filename = str(tmp_path / 'restricted_missing_nstates')
        filename = scf_drv.comm.bcast(filename, root=mpi_master())

        scf_drv.filename = filename
        scf_drv.compute(mol, bas)
        scf_drv.restart = True
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']
        scf_results['filename'] = filename

        reference_drv = LinearResponseEigenSolver()
        reference_drv.ostream.mute()
        reference_drv.filename = filename
        reference_drv.conv_thresh = 1.0e-5
        reference_drv.max_iter = 60
        reference_drv.nstates = 2
        reference_drv.initial_guess_multiplier = 2
        reference_drv.max_subspace_dim = 8
        reference_drv.ri_auxiliary_basis = 'def2-universal-jfit'
        reference_drv.compute(mol, bas, scf_results)

        if reference_drv.rank == mpi_master():
            with h5py.File(reference_drv.checkpoint_file, 'a') as h5f:
                del h5f['nstates']

        reference_drv.comm.barrier()

        restarted_drv = LinearResponseEigenSolver()
        restarted_drv.ostream.mute()
        restarted_drv.filename = filename
        restarted_drv.checkpoint_file = reference_drv.checkpoint_file
        restarted_drv.conv_thresh = 1.0e-5
        restarted_drv.max_iter = 60
        restarted_drv.nstates = 2
        restarted_drv.initial_guess_multiplier = 2
        restarted_drv.max_subspace_dim = 8
        restarted_drv.ri_auxiliary_basis = 'def2-universal-jfit'
        restarted_drv.restart = True

        restarted_results = restarted_drv.compute(mol, bas, scf_results)

        assert restarted_drv.restart is True

        if restarted_drv.rank == mpi_master():
            assert restarted_results[
                'number_of_states'] == restarted_drv.nstates

    @pytest.mark.parametrize('ri', [False, True])
    def test_get_e2_h2o(self, ri):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        e2_matrix = lr_drv.get_e2(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert e2_matrix.shape == (20, 20)
            assert np.max(np.abs(e2_matrix - e2_matrix.T)) < 1.0e-10
            assert np.all(np.diag(e2_matrix) > 0.0)

    @pytest.mark.parametrize('ri', [False, True])
    def test_get_e2_cpcm_matches_rpa(self, ri):
        """get_e2 must build the same E[2] as the RPA solver when CPCM is on.

        Regression test: get_e2 used to skip the solvation sanity checks, so
        _init_cpcm was a no-op and the CPCM response-Fock contribution was
        missing from E[2] (only the solvated MOs entered).  Diagonalizing E[2]
        then disagreed with the regular RPA by ~1e-2 a.u. in the presence of
        CPCM.  This test compares the full eigenvalue spectrum.
        """

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)
        norb = bas.get_dimensions_of_basis()
        nocc = mol.number_of_alpha_occupied_orbitals(bas)
        n_exc = nocc * (norb - nocc)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_epsilon = 10.0
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        # regular RPA over the full excitation space (reference)
        rpa_drv = LinearResponseEigenSolver()
        rpa_drv.ostream.mute()
        rpa_drv.nstates = n_exc
        rpa_drv.max_iter = 200
        rpa_drv.initial_guess_multiplier = 4
        rpa_drv.max_subspace_dim = 2 * n_exc + 10
        rpa_drv.cpcm_optical_epsilon = 5.0
        rpa_drv.conv_thresh = 1.0e-7
        rpa_results = rpa_drv.compute(mol, bas, scf_results)

        # E[2] via get_e2 + diagonalization of E2 v = w S v, S = diag(I,-I)
        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.cpcm_optical_epsilon = 5.0
        e2_matrix = lr_drv.get_e2(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            S = np.zeros_like(e2_matrix)
            S[:n_exc, :n_exc] = np.eye(n_exc)
            S[n_exc:, n_exc:] = -np.eye(n_exc)
            w_all = np.linalg.eigvals(S @ e2_matrix)
            assert np.max(np.abs(np.imag(w_all))) < 1.0e-5
            w_pos = np.sort(np.real(w_all[w_all > 1.0e-6]))
            w_rpa = np.sort(rpa_results['eigenvalues'])
            assert w_pos.shape == w_rpa.shape == (n_exc,)
            assert np.max(np.abs(w_pos - w_rpa)) < 1.0e-6

    @pytest.mark.parametrize('ri', [False, True])
    def test_get_e2_gostshyp_matches_rpa(self, ri):
        """get_e2 must build the same E[2] as the RPA solver when GOSTSHYP is on.

        Mirrors test_get_e2_cpcm_matches_rpa for hydrostatic pressure; before
        the get_e2 sanity-check fix the GOSTSHYP response-Fock contribution was
        missing and the eigenvalues disagreed with RPA by ~5e-4.
        """

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)
        norb = bas.get_dimensions_of_basis()
        nocc = mol.number_of_alpha_occupied_orbitals(bas)
        n_exc = nocc * (norb - nocc)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.pressure = 20000.0
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        rpa_drv = LinearResponseEigenSolver()
        rpa_drv.ostream.mute()
        rpa_drv.nstates = n_exc
        rpa_drv.max_iter = 200
        rpa_drv.initial_guess_multiplier = 4
        rpa_drv.max_subspace_dim = 2 * n_exc + 10
        rpa_drv.conv_thresh = 1.0e-7
        rpa_results = rpa_drv.compute(mol, bas, scf_results)

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        e2_matrix = lr_drv.get_e2(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            S = np.zeros_like(e2_matrix)
            S[:n_exc, :n_exc] = np.eye(n_exc)
            S[n_exc:, n_exc:] = -np.eye(n_exc)
            w_all = np.linalg.eigvals(S @ e2_matrix)
            assert np.max(np.abs(np.imag(w_all))) < 1.0e-5
            w_pos = np.sort(np.real(w_all[w_all > 1.0e-6]))
            w_rpa = np.sort(rpa_results['eigenvalues'])
            assert w_pos.shape == w_rpa.shape == (n_exc,)
            assert np.max(np.abs(w_pos - w_rpa)) < 1.0e-6

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize('ri', [False, True])
    def test_compute_rejects_gostshyp_with_subcomms(self, ri):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.pressure = 20000.0
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 1
        lr_drv.use_subcomms = True

        with pytest.raises(
                VeloxChemError,
                match='LinearSolver: Cannot use subcomms with GOSTSHYP'):
            lr_results_not_used = lr_drv.compute(mol, bas, scf_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_solvation_model_with_pressure(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        # empty scf results just for testing
        scf_results = {}

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.solvation_model = 'cpcm'
        lr_drv.pressure = 20000.0

        with pytest.raises(
                VeloxChemError,
                match="LinearResponseEigenSolver: The 'solvation_model' option "
                      "is incompatible with GOSTSHYP"):
            lr_results_not_used = lr_drv.compute(mol, bas, scf_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_get_e2_rejects_solvation_model_with_pressure(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        # empty scf results just for testing
        scf_results = {}

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.solvation_model = 'smd'
        lr_drv.pressure = 20000.0

        with pytest.raises(
                VeloxChemError,
                match="LinearResponseEigenSolver: The 'solvation_model' option "
                      "is incompatible with GOSTSHYP"):
            e2_matrix_not_used = lr_drv.get_e2(mol, bas, scf_results)

    def test_guess_and_preconditioner_helpers(self):

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 2
        lr_drv.initial_guess_multiplier = 2
        lr_drv.guess_scaling_threshold = 5

        full_guesses = lr_drv._initial_excitations(
            2,
            np.array([-0.9, -0.4, 0.2]),
            1,
            3,
        )
        assert len(full_guesses) == 2
        full_guess_matrix = full_guesses[0].get_full_matrix()
        if lr_drv.rank == mpi_master():
            assert full_guess_matrix.shape == (2, 2)

        lr_drv.core_excitation = True
        lr_drv.num_core_orbitals = 1
        core_guesses = lr_drv._initial_excitations(
            1,
            np.array([-0.9, -0.4, 0.2]),
            2,
            3,
        )
        assert len(core_guesses) == 1

        lr_drv.core_excitation = False
        lr_drv.restricted_subspace = True
        lr_drv.num_core_orbitals = 1
        lr_drv.num_valence_orbitals = 1
        lr_drv.num_virtual_orbitals = 1
        restricted_guesses = lr_drv._initial_excitations(
            1, np.array([-0.9, -0.4, 0.2, 0.6]), 2, 4)
        assert len(restricted_guesses) == 2

        precond = lr_drv._get_precond(np.array([-0.9, -0.4, 0.2, 0.6]), 2, 4,
                                      0.1)
        full_precond = precond.get_full_matrix()
        if lr_drv.rank == mpi_master():
            assert full_precond.shape == (2, 2)
            assert np.all(np.isfinite(full_precond))

        nrows = 16
        pa = np.arange(2.0, 2.0 + nrows)
        pb = np.ones(nrows)
        dist_precond = DistributedArray(
            np.column_stack((pa, pb)),
            lr_drv.comm,
        )
        v_rg = np.arange(1.0, nrows + 1.0)
        v_ru = np.arange(31.0, 31.0 + nrows)
        v_in = DistributedArray(
            np.column_stack((v_rg, v_ru)),
            lr_drv.comm,
        )

        v_out = lr_drv._preconditioning(dist_precond, v_in)
        full_v_out = v_out.get_full_matrix()

        if lr_drv.rank == mpi_master():
            expected_v_out = np.column_stack(
                (pa * v_rg + pb * v_ru, pb * v_rg + pa * v_ru))
            assert np.allclose(full_v_out, expected_v_out)

        lr_drv.norm_thresh = 0.5
        trial_vectors = {
            0: DistributedArray(
                np.column_stack((np.arange(1.0, nrows + 1.0), np.zeros(nrows))),
                lr_drv.comm,
            ),
            1: DistributedArray(
                np.column_stack(
                    (np.zeros(nrows), np.arange(101.0, 101.0 + nrows))),
                lr_drv.comm,
            ),
            2: DistributedArray(
                np.zeros((nrows, 2)),
                lr_drv.comm,
            ),
        }
        diag_precond = {
            k: DistributedArray(
                np.column_stack((np.ones(nrows), np.zeros(nrows))),
                lr_drv.comm,
            ) for k in (0, 1, 2)
        }

        new_ger, new_ung = lr_drv._precond_trials(trial_vectors, diag_precond)
        full_new_ger = new_ger.get_full_matrix()
        full_new_ung = new_ung.get_full_matrix()

        if lr_drv.rank == mpi_master():
            assert np.allclose(full_new_ger,
                               np.arange(1.0, nrows + 1.0).reshape(-1, 1))
            assert np.allclose(full_new_ung,
                               np.arange(101.0, 101.0 + nrows).reshape(-1, 1))

    @pytest.mark.parametrize('ri', [False, True])
    def test_unstable_reference_is_flagged(self, ri):

        # ethylene twisted by 89.9 degrees: the RHF solution is unstable towards
        # complex orbitals; the planar molecule is stable
        # NOTE: 89.9 rather than 90 degrees. At 90 degrees the highest occupied
        # orbital is doubly degenerate, which of the pair the SCF occupies
        # depends on rounding, and plain DIIS (which RI uses) then converges
        # to either of two RHF solutions, 31 mEh apart, depending on the number
        # of threads. At 89.9 degrees the lowest imaginary eigenvalue is
        # -4.9e-3, with RI and four-center alike.

        found = {}

        for twist, h_xy in [(0, ('0.920851', '0.000000')),
                            (89.9, ('0.001607', '0.920850'))]:
            xyz_string = f"""6
            xyz
            C   0.000000   0.000000  -0.665000
            C   0.000000   0.000000   0.665000
            H   0.920851   0.000000  -1.229298
            H  -0.920851   0.000000  -1.229298
            H   {h_xy[0]}   {h_xy[1]}   1.229298
            H  -{h_xy[0]}  -{h_xy[1]}   1.229298
            """
            mol = Molecule.read_xyz_string(xyz_string)
            bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

            scf_drv = ScfRestrictedDriver()
            scf_drv.ri = ri
            scf_drv.ostream.mute()
            scf_results = scf_drv.compute(mol, bas)

            if scf_drv.rank == mpi_master():
                if ri:
                    assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                    assert not scf_results['ri_fallback']

            lr_drv = LinearResponseEigenSolver()
            lr_drv.ostream.mute()
            lr_drv.nstates = 1
            lr_drv.compute(mol, bas, scf_results)

            found[twist] = lr_drv._reference_instability

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert found[0] is None
            assert len(found[89.9]) == 1
            rotation, lowest = found[89.9][0]
            assert rotation == 'imaginary'
            assert lowest < -0.003
