import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.dispersionmodel import DispersionModel
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.vibrationalanalysis import VibrationalAnalysis


class TestHessianFrequencies:

    def run_hessian_freq(self, dispersion, ref_vib_freqs):

        mol = Molecule.read_xyz_string("""14
        xyz
        O     -1.187595563182     1.411836248568    -2.162374839461
        C     -0.439137157819     0.932428966353    -1.332071820972
        N     -0.614791063868     1.150434967949     0.070819254318
        C      0.219478691040     0.619068575821     1.041838969448
        N      0.136993213844     0.762470462293     2.321996830270
        C      1.293188731981    -0.030588456682     2.809691691040
        O      1.617527294188    -0.185955707074     3.962814461787
        N      2.048135131575    -0.646588149798     1.641461874103
        C      1.382415205242    -0.238515458278     0.625198103541
        C      1.637777746106    -0.517149891604    -0.834247645230
        O      2.543913428485    -1.193727604255    -1.285579168444
        N      0.674840114409     0.114341639424    -1.679807817737
        H     -1.403361956671     1.731970220005     0.353632415076
        H      0.790123185526    -0.034140813214    -2.682359312499
        """)
        bas = MolecularBasis.read(mol, 'sto-3g')

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.dispersion = dispersion
        scf_drv.compute(mol, bas)

        vibanalysis_drv = VibrationalAnalysis(scf_drv)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            n_freqs = len(ref_vib_freqs)
            calc_vib_freqs = vibanalysis_drv.vib_frequencies
            assert np.max(np.abs(calc_vib_freqs[:n_freqs] -
                                 ref_vib_freqs)) < 0.1

    @pytest.mark.solvers
    @pytest.mark.skipif(not DispersionModel.is_available(),
                        reason='dftd4-python not available')
    def test_hessian_freq_d4(self):

        dispersion = True
        ref_vib_freqs = np.array([52.90, 57.96, 154.53, 199.51, 275.38])
        self.run_hessian_freq(dispersion, ref_vib_freqs)

    @pytest.mark.timeconsuming
    def test_hessian_freq(self):

        dispersion = False
        ref_vib_freqs = np.array([37.77, 42.84, 149.94, 195.22, 274.18])
        self.run_hessian_freq(dispersion, ref_vib_freqs)
