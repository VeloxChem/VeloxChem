import numpy as np
import h5py
from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
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
            assert restarted_results['number_of_states'] == restarted_drv.nstates

    def test_guess_and_preconditioner_helpers(self):

        lr_drv = LinearResponseUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 1
        lr_drv.initial_guess_multiplier = 2

        guesses = lr_drv._initial_excitations(
            1, (np.array([-0.8, 0.1]), np.array([-0.7, 0.2])), (1, 1), 2)
        assert len(guesses) == 2
        assert guesses[0].data.shape == (2, 2)
        assert guesses[1].data.shape == (2, 2)

        lr_drv.core_excitation = True
        lr_drv.num_core_orbitals = 1
        core_guesses = lr_drv._initial_excitations(
            1, (np.array([-0.9, -0.4, 0.2]), np.array([-0.8, -0.3, 0.3])), (2, 2),
            3)
        assert len(core_guesses) == 2

        precond = lr_drv._get_precond(
            (np.array([-0.9, -0.4, 0.2]), np.array([-0.8, -0.3, 0.3])), (2, 2), 3,
            0.1)
        assert precond.data.shape == (2, 2)
        assert np.all(np.isfinite(precond.data))

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

        with pytest.raises(AssertionError, match='not implemented for nonlinear'):
            nonlinear_drv.compute(mol, bas, scf_results)

        restricted_drv = LinearResponseUnrestrictedEigenSolver()
        restricted_drv.ostream.mute()
        restricted_drv.nstates = 1
        restricted_drv.restricted_subspace = True

        with pytest.raises(AssertionError,
                           match='restricted_subspace not implemented'):
            restricted_drv.compute(mol, bas, scf_results)
