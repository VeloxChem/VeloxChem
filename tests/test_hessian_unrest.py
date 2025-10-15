import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.vibrationalanalysis import VibrationalAnalysis


class TestScfUnrestrictedHessian:

    def run_hessian_unrest(self, xcfun_label, ref_vib_freqs):

        mol = Molecule.read_xyz_string("""4
        xyz
        C       -0.8246083505   -2.6364102075   -0.0255277471
        H       -0.4137775504   -3.6364807662   -0.1702755824
        H       -1.8972335969   -2.4672317696   -0.1290827158
        H       -0.1628531453   -1.8055908756    0.2228478391
        """)
        mol.set_multiplicity(2)

        bas = MolecularBasis.read(mol, 'def2-svp')

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.compute(mol, bas)

        vibanalysis_drv = VibrationalAnalysis(scf_drv)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            calc_vib_freqs = vibanalysis_drv.vib_frequencies
            assert np.max(np.abs(calc_vib_freqs - ref_vib_freqs)) < 0.1

    @pytest.mark.solvers
    def test_hessian_unrest(self):

        xcfun_label = 'b3lyp'

        ref_vib_freqs = np.array(
            [469.94, 1373.71, 1373.75, 3096.45, 3296.31, 3296.34])

        self.run_hessian_unrest(xcfun_label, ref_vib_freqs)
