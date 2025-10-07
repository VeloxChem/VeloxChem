import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.lreigensolverunrest import LinearResponseUnrestrictedEigenSolver


@pytest.mark.solvers
class TestUnrestrictedRPA:

    def run_rpa(self, xcfun_label, basis_label, ref_exc_enes, ref_osc_str, tol):

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

        self.run_rpa('pbe0', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5)
