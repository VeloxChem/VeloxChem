import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestScfRestrictedSemiNumDriver:

    def run_scf_rest(self, xcfun_label, ref_scf_energy, tol):

        #xyz_string = """6
        #xyz
        #O   -0.1858140  -1.1749469   0.7662596
        #H   -0.1285513  -0.8984365   1.6808606
        #H   -0.0582782  -0.3702550   0.2638279
        #O    0.1747051  11.1050002  -0.7244430
        #H   -0.5650842  11.3134964  -1.2949455
        #H    0.9282185  11.0652990  -1.3134026
        #"""
        xyz_string = """24
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
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)

        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.ri_coulomb = True
        scf_drv.semi_num_exchange = True
        scf_results = scf_drv.compute(mol, bas, min_bas)
        
        if scf_drv.rank == mpi_master():
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol

    def test_b3lyp(self):

        self.run_scf_rest('b3lyp', -152.7162995916, 1.0e-6)
