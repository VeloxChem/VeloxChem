import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.tdaeigensolverunrest import TdaUnrestrictedEigenSolver


@pytest.mark.solvers
class TestUnrestrictedTDA:

    def run_tda(self, xcfun_label, basis_label, ref_exc_enes, ref_osc_str, tol):

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

        lr_drv = TdaUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 10
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    def test_camb3lyp_svp(self):

        # vlxtag: UKS, Absorption, TDA

        ref_exc_enes = np.array([
            0.09065604, 0.09500459, 0.21881146, 0.47239371, 0.48842790,
            0.51728556, 0.52413401, 0.55026290, 0.55530527, 0.59463209
        ])

        ref_osc_str = np.array([
            0.0013, 0.2065, 0.0000, 0.0325, 0.1953, 0.0015, 0.0000, 0.0256,
            0.0701, 0.1727
        ])

        self.run_tda('cam-b3lyp', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5)
