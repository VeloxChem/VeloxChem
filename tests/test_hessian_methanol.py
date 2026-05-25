import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.vibrationalanalysis import VibrationalAnalysis


class TestHessianMethanol:

    def run_hessian_methanol(self, numerical_hessian):

        mol = Molecule.read_xyz_string("""6
        methanol
        H           1.233468597699        0.036555995513        0.834256310250
        C           0.688882904364        0.008190488215       -0.121375210387
        H           1.004038901495        0.896197692302       -0.704487055168
        H           1.031928585196       -0.893040561619       -0.667271522601
        O          -0.677956408466       -0.007094311055        0.168195416310
        H          -1.160565041296       -0.031866350457       -0.663142361588
        """)
        bas = MolecularBasis.read(mol, 'def2-svp')

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = 'pbe0'
        scf_drv.ostream.mute()
        scf_drv.compute(mol, bas)

        vibanalysis_drv = VibrationalAnalysis(scf_drv)
        vibanalysis_drv.numerical_hessian = numerical_hessian
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            ref_vib_freqs = np.array([
                340.605, 1076.485, 1145.587, 1177.917, 1369.167, 1470.522,
                1475.389, 1496.336, 2992.628, 3054.415, 3142.981, 3897.690
            ])

            max_rel_diff_freqs = np.max(
                np.abs(vibanalysis_drv.vib_frequencies / ref_vib_freqs - 1.0))

            assert max_rel_diff_freqs < 0.02

            ref_ir_intens = np.array([
                128.9167, 66.1455, 65.0787, 0.5474, 28.7916, 2.5639, 6.1060,
                6.6203, 63.3375, 72.5923, 25.9392, 36.6984
            ])

            max_rel_diff_ir_intens = np.max(
                np.abs(vibanalysis_drv.ir_intensities / ref_ir_intens - 1.0))

            assert max_rel_diff_ir_intens < 0.02

    @pytest.mark.timeconsuming
    def test_numerical_hessian_methanol(self):

        self.run_hessian_methanol(numerical_hessian=True)

    @pytest.mark.solvers
    def test_analytical_hessian_methanol(self):

        self.run_hessian_methanol(numerical_hessian=False)
