import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver


@pytest.mark.solvers
class TestRPA:

    def run_rpa(self, xcfun_label, basis_label, ref_exc_enes, ref_osc_str, tol):

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
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    def test_hf_svp(self):

        # vlxtag: RHF, Absorption, TDHF

        ref_exc_enes = np.array(
            [0.33973039, 0.40464346, 0.43325535, 0.49805052, 0.55325390])

        ref_osc_str = np.array(
            [0.023665, 0.000000, 0.097765, 0.086454, 0.291919])

        self.run_rpa('hf', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-6)

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
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseEigenSolver()
        lr_drv.filename = filename
        lr_drv.ostream.mute()

        lr_drv.nstates = 5
        lr_results_not_used = lr_drv.compute(mol, bas, scf_results)

        lr_drv.restart = True
        lr_drv.nstates = 10
        lr_results_first = lr_drv.compute(mol, bas, scf_results)

        lr_drv.restart = False
        lr_drv.nstates = 10
        lr_results_second = lr_drv.compute(mol, bas, scf_results)

        if scf_drv.rank == mpi_master():
            for key in [
                    'eigenvalues',
                    'oscillator_strengths',
                    'rotatory_strengths',
            ]:
                assert np.max(
                    np.abs(lr_results_first[key] -
                           lr_results_second[key])) < 1e-10
