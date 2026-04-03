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


@pytest.mark.solvers
class TestRPA:

    def run_rpa(self,
                xcfun_label,
                basis_label,
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
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

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
                AssertionError,
                match="Molecule: Invalid multiplicity for restricted"):
            lr_results_not_used = lr_drv.compute(mol, bas, scf_results)

    def test_hf_svp(self):

        # vlxtag: RHF, Absorption, TDHF

        ref_exc_enes = np.array(
            [0.33973039, 0.40464346, 0.43325535, 0.49805052, 0.55325390])

        ref_osc_str = np.array(
            [0.023665, 0.000000, 0.097765, 0.086454, 0.291919])

        self.run_rpa('hf',
                     'def2-svp',
                     ref_exc_enes,
                     ref_osc_str,
                     1.0e-6,
                     max_subspace_dim=50)

    def test_slda_svp(self):

        # vlxtag: RKS, Absorption, TDDFT

        ref_exc_enes = np.array(
            [0.27459103, 0.34805152, 0.35179468, 0.43068396, 0.51123472])

        ref_osc_str = np.array(
            [0.018005, 0.000000, 0.075418, 0.059909, 0.258165])

        self.run_rpa('slda', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5)

    def test_b3lyp_svp(self):

        # vlxtag: RKS, Absorption, TDDFT

        ref_exc_enes = np.array(
            [0.27940792, 0.35031697, 0.36276301, 0.43718376, 0.51430795])

        ref_osc_str = np.array(
            [0.018243, 0.000000, 0.078326, 0.061329, 0.272933])

        self.run_rpa('b3lyp', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5)

    def test_camb3lyp_svp(self):

        # vlxtag: RKS, Absorption, TDDFT

        ref_exc_enes = np.array(
            [0.28249530, 0.35533592, 0.36664012, 0.44312170, 0.51548914])

        ref_osc_str = np.array(
            [0.018219, 0.000000, 0.077391, 0.059649, 0.273690])

        self.run_rpa('cam-b3lyp', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5)

    def test_camb3lyp_tzvp(self):

        # vlxtag: RKS, Absorption, TDDFT

        ref_exc_enes = np.array(
            [0.27984870, 0.35148791, 0.36310452, 0.43770042, 0.50331834])

        ref_osc_str = np.array(
            [0.033878, 0.000000, 0.099329, 0.053065, 0.229533])

        self.run_rpa('cam-b3lyp', 'def2-tzvp', ref_exc_enes, ref_osc_str,
                     1.0e-5)

    def test_tpssh_svp(self):

        # vlxtag: RKS, Absorption, TDDFT

        ref_exc_enes = np.array(
            [0.28898749, 0.36043436, 0.37287451, 0.44801126, 0.52385373])

        ref_osc_str = np.array(
            [0.020158, 0.000000, 0.085796, 0.074117, 0.272050])

        self.run_rpa('tpssh', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5)

    def test_detach_attach_charges_water_631g(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, '6-31g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 3
        lr_drv.detach_attach = True
        lr_drv.detach_attach_charges = True
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            detachment = lr_results['detachment_charges']
            attachment = lr_results['attachment_charges']
            ref_detachment = np.array([
                [-1.0001863412, 9.31706e-05, 9.31706e-05],
                [-1.0004727906, 2.363953e-04, 2.363953e-04],
                [-0.9115710485, -4.42144757e-02, -4.42144758e-02],
            ])
            ref_attachment = np.array([
                [0.2873499273, 0.3563249945, 0.3563250781],
                [0.3025686881, 0.3487157020, 0.3487156098],
                [0.2760202908, 0.3619898252, 0.3619898840],
            ])

            assert detachment.shape == (3, mol.number_of_atoms())
            assert attachment.shape == (3, mol.number_of_atoms())
            assert np.max(np.abs(detachment - ref_detachment)) < 1.0e-6
            assert np.max(np.abs(attachment - ref_attachment)) < 1.0e-6

    def test_esa_water_631g(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, '6-31g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 3
        lr_drv.esa = True
        lr_drv.esa_from_state = 1
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            esa_results = lr_results['esa_results']

            ref_excitation_energies = np.array([0.0706309931, 0.0887848965])
            ref_oscillator_strengths = np.array([0.1379665336, 0.0013975526])
            ref_transition_dipoles = np.array([
                [-0.0794559666, -0.5972010928, 1.6022021187],
                [0.1515690613, -0.0251919303, -0.0018734103],
            ])

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

    def run_rpa_with_ecp(self, ref_exc_enes, ref_osc_str, tol):

        xyz_string = """2
        xyz
        Au 0 0 0
        Au 0 0 2.88
        """
        mol = Molecule.read_xyz_string(xyz_string)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

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

    def test_hf_with_ecp(self):

        # vlxtag: RHF, Absorption, TDHF

        ref_exc_enes = np.array([0.10859702, 0.15291543, 0.16111886])
        ref_osc_str = np.array([0.4439, 0.0000, 0.0025])

        self.run_rpa_with_ecp(ref_exc_enes, ref_osc_str, 1.0e-6)

    def test_checkpoint_and_restart(self, tmp_path):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()

        filename = str(tmp_path / "water_restart")
        # To avoid inconsistency across MPI ranks
        filename = scf_drv.comm.bcast(filename, root=mpi_master())

        scf_drv.filename = filename
        scf_drv.compute(mol, bas)
        scf_drv.restart = True
        scf_results = scf_drv.compute(mol, bas)

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

    def test_restart_allows_missing_nstates(self, tmp_path):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()

        filename = str(tmp_path / 'restricted_missing_nstates')
        filename = scf_drv.comm.bcast(filename, root=mpi_master())

        scf_drv.filename = filename
        scf_drv.compute(mol, bas)
        scf_drv.restart = True
        scf_results = scf_drv.compute(mol, bas)
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

    def test_get_e2_h2o(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        e2_matrix = lr_drv.get_e2(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert e2_matrix.shape == (20, 20)
            assert np.max(np.abs(e2_matrix - e2_matrix.T)) < 1.0e-10
            assert np.all(np.diag(e2_matrix) > 0.0)

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
