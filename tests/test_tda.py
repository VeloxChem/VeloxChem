import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver


@pytest.mark.solvers
class TestTDA:

    def run_tda(self, xcfun_label, basis_label, ref_exc_enes, ref_osc_str, tol):

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

        lr_drv = TdaEigenSolver()
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

        # vlxtag: RHF, Absorption, CIS

        ref_exc_enes = np.array(
            [0.34189725, 0.40717708, 0.43577895, 0.50150019, 0.55485041])

        ref_osc_str = np.array([0.0229, 0.0000, 0.1039, 0.0975, 0.3060])

        self.run_tda('hf', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-6)

    def test_b3lyp_svp(self):

        # vlxtag: RKS, Absorption, TDA

        ref_exc_enes = np.array(
            [0.28045483, 0.35053288, 0.36504833, 0.43904962, 0.51570229])

        ref_osc_str = np.array([0.0180, 0.0000, 0.0854, 0.0695, 0.3006])

        self.run_tda('b3lyp', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5)

    def test_camb3lyp_svp(self):

        # vlxtag: RKS, Absorption, TDA

        ref_exc_enes = np.array(
            [0.28370223, 0.35557650, 0.36883370, 0.44496286, 0.51683587])

        ref_osc_str = np.array(
            [0.017966, 0.000000, 0.084011, 0.068027, 0.300331])

        self.run_tda('cam-b3lyp', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5)

    def test_camb3lyp_tzvp(self):

        # vlxtag: RKS, Absorption, TDA

        ref_exc_enes = np.array(
            [0.28085567, 0.35172399, 0.36474654, 0.43912355, 0.50442105])

        ref_osc_str = np.array(
            [0.034572, 0.000000, 0.107302, 0.059653, 0.249276])

        self.run_tda('cam-b3lyp', 'def2-tzvp', ref_exc_enes, ref_osc_str,
                     1.0e-5)
