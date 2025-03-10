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
        scf_drv.conv_thresh = 1.0e-5
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
        
        ref_scf_energy = -679.6642908122
        tol_energy = 1.0e-4

        ref_scf_grad = np.array([
            [0.005242277507, -0.013389592407, 0.000000000207],
            [-0.013571482652, 0.001982786219, 0.000000000082],
            [0.000943078535, 0.005079442325, 0.000000000263],
            [0.002733515659, -0.008339325258, -0.000000000625],
            [-0.002557913955, -0.003331008468, -0.000000000155],
            [-0.002973139885, 0.009249243520, -0.000000000134],
            [0.003034870060, -0.007666853630, -0.000000000551],
            [-0.004728363228, 0.009380264234, 0.000000000405],
            [0.000925234042, 0.001588007019, 0.000000000697],
            [-0.001123130832, -0.001065272766, -0.000000000294],
            [0.003638552597, 0.005998941873, 0.000000000278],
            [0.002930643512, -0.003498997249, -0.000000000837],
            [-0.002557785393, 0.003774512513, 0.000000000286],
            [0.003625755498, 0.002874887696, 0.000000000149],
            [0.009863853285, 0.003873312821, -0.000000000113],
            [0.006467040241, 0.006095849228, -0.000000000090],
            [-0.005768983575, 0.001686344798, -0.006795523826],
            [-0.005768984237, 0.001686344964, 0.006795524536],
            [-0.003746835573, -0.008472901572, 0.000000000122],
            [0.005479949832, 0.000156457072, -0.006894223131],
            [0.005479949426, 0.000156457139, 0.006894222634],
            [0.004451864845, -0.008051350930, 0.000000000077],
            [-0.005988005912, 0.000133050116, 0.006732158194],
            [-0.005988005959, 0.000133050092, -0.006732158173],
        ])
        tol_grad = 1.0e-5

        self.run_scf_grad(mol_xyz, 'blyp', ref_scf_energy, ref_scf_grad, tol_energy,
                          tol_grad)

    def test_ri_blyp_dimer(self):

        mol_xyz = """2
            xyz
            Cu 0.1000 0.3000  0.4000
            Cu 0.4000 0.5100  2.1700
        """
        
        ref_scf_energy = -3280.6759055800
        tol_energy = 1.0e-4

        # taken fron numerical gradient
        ref_scf_grad = np.array([
            [ 0.031992  ,  0.02250086,  0.18903518],
            [-0.031992  , -0.02250086, -0.18903518]
        ])
        # this need to be checked
        tol_grad = 2.0e-3

        self.run_scf_grad(mol_xyz, 'blyp', ref_scf_energy, ref_scf_grad, tol_energy,
                          tol_grad)

