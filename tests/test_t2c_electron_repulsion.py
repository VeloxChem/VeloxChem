from pathlib import Path
import numpy as np

from veloxchem import TwoCenterElectronRepulsionDriver
from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import SubMatrix


class TestTwoCenterElectronRepulsionDriver:

    def get_data_svp(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-svp')

        return mol, bas

    def get_data_qzvp(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-qzvp')

        return mol, bas

    def test_electron_repulsion_h2o_svp(self):

        mol, bas = self.get_data_svp()

        # compute electron repulsion matrix
        eri_drv = TwoCenterElectronRepulsionDriver()
        eri_mat = eri_drv.compute(mol, bas)

        # load reference kinetic energy data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.svp.int2c2e.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(3)
        basdims = [0, 7, 19, 24]

        # check individual electron repulsion submatrices
        for i, j in zip(indexes[0], indexes[1]):
            sbra = basdims[i]
            ebra = basdims[i + 1]
            sket = basdims[j]
            eket = basdims[j + 1]
            cmat = eri_mat.submatrix((i, j))
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            assert cmat == rmat

        # check full electron repulsion matrix
        fmat = eri_mat.full_matrix()
        fref = SubMatrix([0, 0, 24, 24])
        fref.set_values(np.ascontiguousarray(ref_mat))
        assert fmat == fref

    def test_electron_repulsion_h2o_qzvp(self):

        mol, bas = self.get_data_qzvp()

        # compute electron repulsion matrix
        eri_drv = TwoCenterElectronRepulsionDriver()
        eri_mat = eri_drv.compute(mol, bas)

        # load reference kinetic energy data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.qzvp.int2c2e.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 15, 45, 80, 108, 117]

        # check individual electron repulsion submatrices
        for i, j in zip(indexes[0], indexes[1]):
            sbra = basdims[i]
            ebra = basdims[i + 1]
            sket = basdims[j]
            eket = basdims[j + 1]
            cmat = eri_mat.submatrix((i, j))
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            assert cmat == rmat

        # check full electron repulsion matrix
        fmat = eri_mat.full_matrix()
        fref = SubMatrix([0, 0, 117, 117])
        fref.set_values(np.ascontiguousarray(ref_mat))
        assert fmat == fref
