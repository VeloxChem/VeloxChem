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
        vibanalysis_drv.update_settings(method_dict={'xcfun': 'pbe0'})
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            ref_vib_freqs = np.array([
                340.790, 1076.524, 1145.529, 1177.911, 1369.222, 1470.513,
                1475.423, 1496.391, 2992.404, 3054.107, 3142.856, 3897.554
            ])

            max_rel_diff_freqs = np.max(
                np.abs(vibanalysis_drv.vib_frequencies / ref_vib_freqs - 1.0))

            assert max_rel_diff_freqs < 0.02

    @pytest.mark.timeconsuming
    def test_numerical_hessian_methanol(self):

        self.run_hessian_methanol(numerical_hessian=True)

    @pytest.mark.solvers
    def test_analytical_hessian_methanol(self):

        self.run_hessian_methanol(numerical_hessian=False)
