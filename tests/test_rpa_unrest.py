import numpy as np
import h5py
from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.distributedarray import DistributedArray
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.lreigensolverunrest import LinearResponseUnrestrictedEigenSolver


@pytest.mark.solvers
class TestUnrestrictedRPA:

    def run_rpa(self,
                xcfun_label,
                basis_label,
                ref_exc_enes,
                ref_osc_str,
                tol,
                max_subspace_dim=None):

        xyz_string = """3
        xyz
        O      0.000000   0.000000   0.117790
        H      0.000000   0.755453  -0.471161
        H      0.000000  -0.755453  -0.471161
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(3)

        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 10
        lr_drv.max_subspace_dim = max_subspace_dim
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    def test_hf_svp(self):

        # vlxtag: UHF, Absorption, TDHF

        ref_exc_enes = np.array([
            0.07000019, 0.09444351, 0.23795183, 0.53171161, 0.54518734,
            0.55496039, 0.57527661, 0.58052563, 0.61929121, 0.63384529
        ])

        ref_osc_str = np.array([
            0.1803, 0.0013, 0.0000, 0.2184, 0.0358, 0.0002, 0.0187, 0.0000,
            0.0789, 0.0253
        ])

        self.run_rpa('hf', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-6)

    def test_pbe0_svp(self):

        # vlxtag: UKS, Absorption, TDDFT

        ref_exc_enes = np.array([
            0.08832630, 0.09472770, 0.22512802, 0.48348997, 0.48361504,
            0.51691424, 0.52828051, 0.54298284, 0.56229157, 0.60539714
        ])

        ref_osc_str = np.array([
            0.1596, 0.0012, 0.0000, 0.1985, 0.0338, 0.0072, 0.0000, 0.0308,
            0.0600, 0.1210
        ])

        self.run_rpa('pbe0',
                     'def2-svp',
                     ref_exc_enes,
                     ref_osc_str,
                     1.0e-5,
                     max_subspace_dim=120)

    def run_rpa_with_ecp(self, ref_exc_enes, ref_osc_str, tol):

        xyz_string = """2
        xyz
        Au 0 0 0
        H  0 0 1.55
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 5
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    def test_hf_with_ecp(self):

        # vlxtag: UHF, Absorption, TDHF

        ref_exc_enes = np.array(
            [0.08148954, 0.08148954, 0.09916297, 0.09916297, 0.12983930])
        ref_osc_str = np.array([0.0000, 0.0000, 0.0001, 0.0001, 0.0119])

        self.run_rpa_with_ecp(ref_exc_enes, ref_osc_str, 1.0e-6)

    def test_core_excitation_water_cation_sto3g(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 3
        lr_drv.core_excitation = True
        lr_drv.num_core_orbitals = 1
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            ref_exc_enes = np.array([19.686416638, 20.394042142, 20.435119894])
            ref_osc_str = np.array([0.051629296, 0.002790081, 0.013572433])

            assert np.max(np.abs(lr_results['eigenvalues'] -
                                 ref_exc_enes)) < 1.0e-7
            assert np.max(
                np.abs(lr_results['oscillator_strengths'] -
                       ref_osc_str)) < 1.0e-5

    def test_nto_and_detach_attach_water_cation_sto3g(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 2
        lr_drv.nto = True
        lr_drv.detach_attach = True
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            ref_nto_lambdas_a = np.array([
                [7.92386e-05, 0.0],
                [4.3593022e-03, 0.0],
            ])
            ref_nto_lambdas_b = np.array([
                [1.0002985115, 0.0, 0.0],
                [1.0029741042, 0.0, 0.0],
            ])

            assert 'nto_lambdas_a' in lr_results
            assert 'nto_lambdas_b' in lr_results
            assert len(lr_results['nto_lambdas_a']) == 2
            assert len(lr_results['nto_lambdas_b']) == 2
            assert np.max(
                np.abs(
                    np.array(lr_results['nto_lambdas_a']) -
                    ref_nto_lambdas_a)) < 1.0e-6
            assert np.max(
                np.abs(
                    np.array(lr_results['nto_lambdas_b']) -
                    ref_nto_lambdas_b)) < 1.0e-6
            assert lr_results['number_of_states'] == 2
            assert len(lr_results['eigenvalues']) == 2
            assert np.all(np.isfinite(lr_results['eigenvalues']))

    def test_checkpoint_and_restart(self, tmp_path):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1.0)
        mol.set_multiplicity(2)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()

        filename = str(tmp_path / "water_restart")
        # To avoid inconsistency across MPI ranks
        filename = scf_drv.comm.bcast(filename, root=mpi_master())

        scf_drv.filename = filename
        scf_drv.compute(mol, bas)
        scf_drv.restart = True
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseUnrestrictedEigenSolver()
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
        mol.set_charge(1)
        mol.set_multiplicity(2)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()

        filename = str(tmp_path / 'unrestricted_missing_nstates')
        filename = scf_drv.comm.bcast(filename, root=mpi_master())

        scf_drv.filename = filename
        scf_drv.compute(mol, bas)
        scf_drv.restart = True
        scf_results = scf_drv.compute(mol, bas)
        scf_results['filename'] = filename

        reference_drv = LinearResponseUnrestrictedEigenSolver()
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

        restarted_drv = LinearResponseUnrestrictedEigenSolver()
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

    def test_guess_and_preconditioner_helpers(self):

        lr_drv = LinearResponseUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 1
        lr_drv.initial_guess_multiplier = 2

        guesses = lr_drv._initial_excitations(
            1, (np.array([-0.8, 0.1]), np.array([-0.7, 0.2])), (1, 1), 2)
        assert len(guesses) == 2
        full_guess_0 = guesses[0].get_full_matrix()
        full_guess_1 = guesses[1].get_full_matrix()
        if lr_drv.rank == mpi_master():
            assert full_guess_0.shape == (2, 2)
            assert full_guess_1.shape == (2, 2)

        lr_drv.core_excitation = True
        lr_drv.num_core_orbitals = 1
        core_guesses = lr_drv._initial_excitations(
            1, (np.array([-0.9, -0.4, 0.2]), np.array([-0.8, -0.3, 0.3])),
            (2, 2), 3)
        assert len(core_guesses) == 2

        precond = lr_drv._get_precond(
            (np.array([-0.9, -0.4, 0.2]), np.array([-0.8, -0.3, 0.3])), (2, 2),
            3, 0.1)
        full_precond = precond.get_full_matrix()
        if lr_drv.rank == mpi_master():
            assert full_precond.shape == (2, 2)
            assert np.all(np.isfinite(full_precond))

        nrows = 20
        pa = np.arange(2.0, 2.0 + nrows)
        pb = np.ones(nrows)
        dist_precond = DistributedArray(
            np.column_stack((pa, pb)),
            lr_drv.comm,
        )
        v_rg = np.arange(1.0, nrows + 1.0)
        v_ru = np.arange(41.0, 41.0 + nrows)
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
                    (np.zeros(nrows), np.arange(201.0, 201.0 + nrows))),
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
                               np.arange(201.0, 201.0 + nrows).reshape(-1, 1))

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_rejects_unsupported_modes(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        nonlinear_drv = LinearResponseUnrestrictedEigenSolver()
        nonlinear_drv.ostream.mute()
        nonlinear_drv.nstates = 1
        nonlinear_drv.nonlinear = True

        with pytest.raises(AssertionError,
                           match='not implemented for nonlinear'):
            nonlinear_drv.compute(mol, bas, scf_results)

        restricted_drv = LinearResponseUnrestrictedEigenSolver()
        restricted_drv.ostream.mute()
        restricted_drv.nstates = 1
        restricted_drv.restricted_subspace = True

        with pytest.raises(AssertionError,
                           match='restricted_subspace not implemented'):
            restricted_drv.compute(mol, bas, scf_results)
