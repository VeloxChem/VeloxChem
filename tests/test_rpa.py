import numpy as np
import pytest

from veloxchem import mpi_master
from veloxchem import Molecule, MolecularBasis
from veloxchem import ScfRestrictedDriver
from veloxchem import LinearResponseEigenSolver


@pytest.mark.solvers
class TestRPA:

    def run_rpa(self, xcfun_label, ref_exc_enes, ref_osc_str, tol):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.check_multiplicity()

        basis_label = 'def2-svp'
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

    def test_hf(self):

        ref_exc_enes = np.array(
            [0.33973039, 0.40464346, 0.43325535, 0.49805052, 0.55325390])

        ref_osc_str = np.array(
            [0.023665, 0.000000, 0.097765, 0.086454, 0.291919])

        self.run_rpa('hf', ref_exc_enes, ref_osc_str, 1.0e-6)

    def test_slda(self):

        ref_exc_enes = np.array(
            [0.27459103, 0.34805152, 0.35179468, 0.43068396, 0.51123472])

        ref_osc_str = np.array(
            [0.018005, 0.000000, 0.075418, 0.059909, 0.258165])

        self.run_rpa('slda', ref_exc_enes, ref_osc_str, 1.0e-5)

    def test_b3lyp(self):

        ref_exc_enes = np.array(
            [0.27940792, 0.35031697, 0.36276301, 0.43718376, 0.51430795])

        ref_osc_str = np.array(
            [0.018243, 0.000000, 0.078326, 0.061329, 0.272933])

        self.run_rpa('b3lyp', ref_exc_enes, ref_osc_str, 1.0e-5)

    def test_tpssh(self):

        ref_exc_enes = np.array(
            [0.28898749, 0.36043436, 0.37287451, 0.44801126, 0.52385373])

        ref_osc_str = np.array(
            [0.020158, 0.000000, 0.085796, 0.074117, 0.272050])

        self.run_rpa('tpssh', ref_exc_enes, ref_osc_str, 1.0e-5)
