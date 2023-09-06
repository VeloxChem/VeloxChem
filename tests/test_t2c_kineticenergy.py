from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import KineticEnergyDriver
from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import Molecule
from veloxchem.submatrix import SubMatrix
from tester import Tester


class TestKineticEnergyDriver:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')

        here = Path(__file__).parent
        basis_path = str(here.parent / 'basis')
        bas = MolecularBasis.read(mol, 'def2-tzvpp', basis_path, ostream=None)

        return mol, bas

    def test_kinetic_energy_co_tzvpp(self):

        mol_co, bas_tzvpp = self.get_data()

        # compute kinetic energy matrix
        kin_drv = KineticEnergyDriver()
        kin_mat = kin_drv.compute(mol_co, bas_tzvpp)

        # load reference kinetic energy data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.tzvpp.kinetic.energy.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(4)
        basdims = [0, 10, 28, 48, 62]

        # check individual kinetic energy submatrices
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = kin_mat.get_submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
            # compare submatrices
            Tester.compare_submatrices(cmat, rmat)

        # check full kinetic energy matrix
        fmat = kin_mat.get_full_matrix()
        fref = SubMatrix([0, 0, 62, 62])
        fref.set_values(np.ascontiguousarray(ref_mat))
        Tester.compare_submatrices(fmat, fref)
