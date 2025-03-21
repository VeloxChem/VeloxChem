import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver


@pytest.mark.solvers
class TestScfGradientDriverWithRI:

    def run_scf_grad(self, mol_xyz, xcfun_label, ref_scf_energy, ref_scf_grad,
                     tol_energy, tol_grad):

        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(mol_xyz)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.ri_coulomb = True
        scf_results = scf_drv.compute(mol, bas)

        scf_grad_drv = ScfGradientDriver(scf_drv)
        scf_grad_drv.compute(mol, bas, scf_results)
        grad = scf_grad_drv.get_gradient()

        if scf_drv.rank == mpi_master():
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol_energy
            assert np.max(np.abs(ref_scf_grad - grad)) < tol_grad

    def test_ri_blyp(self):

        mol_xyz = """24
            xyz
            O    -0.537797810763    2.548484593229    0.000000000000
            O     3.060053657153   -0.342736262655   -0.000000000000
            N     0.981027030126   -1.286353215374   -0.000000000000
            N    -2.292611960807    0.008027029000    0.000000000000
            N     1.262169862164    1.085372972620    0.000000000000
            N    -1.346854892845   -2.034482811562   -0.000000000000
            C    -0.915797659810    0.195858564021    0.000000000000
            C    -0.379747568515   -1.078485797996   -0.000000000000
            C    -0.121435139475    1.390724912803    0.000000000000
            C     1.846035071709   -0.196767477313   -0.000000000000
            C    -2.478662573702   -1.335323738891   -0.000000000000
            C     1.562698463388   -2.628418219958   -0.000000000000
            C    -3.317085271030    1.051202855269    0.000000000000
            C     2.214919159110    2.200764309972    0.000000000000
            H    -3.468470456205   -1.773774496739   -0.000000000000
            H     0.742628036719   -3.346033930132   -0.000000000000
            H     2.186800646177   -2.766305911953    0.886988671775
            H     2.186800646177   -2.766305911953   -0.886988671775
            H    -2.808345859284    2.015050628143    0.000000000000
            H    -3.944028242129    0.965315311095    0.892919708821
            H    -3.944028242129    0.965315311095   -0.892919708821
            H     1.635263989242    3.121966936999    0.000000000000
            H     2.852229840507    2.147373554317   -0.886283579224
            H     2.852229840507    2.147373554317    0.886283579224
            """

        ref_scf_energy = -679.6640733004708
        tol_energy = 1.0e-3

        ref_scf_grad = np.array([[
            [0.00524178, -0.01338491, 0.00000000],
            [-0.01356798, 0.00198281, -0.00000000],
            [0.00094571, 0.005083, 0.00000000],
            [0.00273102, -0.00834755, 0.00000000],
            [-0.00255591, -0.00333202, -0.00000000],
            [-0.00297843, 0.00924393, -0.00000000],
            [0.00303613, -0.0076707, 0.00000000],
            [-0.00472998, 0.00938578, 0.00000000],
            [0.0009269, 0.00158823, -0.00000000],
            [-0.00112704, -0.00106581, 0.00000000],
            [0.00363025, 0.00599441, -0.00000000],
            [0.00292737, -0.00351007, -0.00000000],
            [-0.00255359, 0.00378607, 0.00000000],
            [0.00362135, 0.0028869, -0.00000000],
            [0.00987238, 0.00387719, 0.00000000],
            [0.00646656, 0.00609804, 0.00000000],
            [-0.00576499, 0.00168603, -0.00678763],
            [-0.00576499, 0.00168603, 0.00678763],
            [-0.00374757, -0.00847489, 0.00000000],
            [0.00547476, 0.00015532, -0.00688561],
            [0.00547476, 0.00015532, 0.00688561],
            [0.00445278, -0.00805418, -0.00000000],
            [-0.00598365, 0.00013237, 0.00672354],
            [-0.00598365, 0.00013237, -0.00672354],
        ]])
        tol_grad = 1.0e-4

        self.run_scf_grad(mol_xyz, 'blyp', ref_scf_energy, ref_scf_grad,
                          tol_energy, tol_grad)

    def test_ri_blyp_dimer(self):

        mol_xyz = """2
            xyz
            Cu 0.1000 0.3000  0.4000
            Cu 0.4000 0.5100  2.1700
        """

        ref_scf_energy = -3280.675905579988
        tol_energy = 5.0e-4

        ref_scf_grad = np.array([[0.03204882, 0.02245542, 0.18911158],
                                 [-0.03204882, -0.02245542, -0.18911158]])
        tol_grad = 5.0e-5

        self.run_scf_grad(mol_xyz, 'blyp', ref_scf_energy, ref_scf_grad,
                          tol_energy, tol_grad)

        ref_scf_energy = -3280.5482717437712
        tol_energy = 5.0e-4

        ref_scf_grad = np.array([[0.03093886, 0.0216664, 0.1825468],
                                 [-0.03093886, -0.0216664, -0.1825468]])
        tol_grad = 5.0e-5

        self.run_scf_grad(mol_xyz, 'tpss', ref_scf_energy, ref_scf_grad,
                          tol_energy, tol_grad)
