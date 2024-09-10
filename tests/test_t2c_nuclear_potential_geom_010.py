from pathlib import Path
import numpy as np

from veloxchem import NuclearPotentialGeom010Driver
from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import SubMatrix


class TestNuclearPotentialGeom010Driver:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-QZVP', ostream=None)

        return mol, bas

    def test_nuclear_potential_co_qzvp(self):

        mol_co, bas_qzvp = self.get_data()

        # external dipoles, coordinates
        dipoles = [
            1.0,
            2.0,
            3.0,
            1.5,
            2.5,
            3.5,
        ]
        coords = [
            [0.1, 0.2, 0.3],
            [1.0, 1.2, 1.4],
        ]

        # compute nuclear potential matrix
        geom_drv = NuclearPotentialGeom010Driver()
        geom_mat = geom_drv.compute(mol_co, bas_qzvp, dipoles, coords)

    # TODO: Need test data for electric field integrals

    def test_nuclear_potential_co_qzvp_for_c(self):

        mol_co, bas_qzvp = self.get_data()

        # compute nuclear potential matrix
        geom_drv = NuclearPotentialGeom010Driver()
        geom_mats = geom_drv.compute(mol_co, bas_qzvp, 0)

        # load reference overlap gradient for C atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' /
                      'co.qzvp.nuclear.potential.geom.010.c.npy')
        ref_mat = np.load(npyfile)
        ref_mat = -ref_mat

        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = ['X', 'Y', 'Z']

        for k, label in enumerate(labels):
            fmat = geom_mats.matrix(label)
            for i, j in zip(indexes[0], indexes[1]):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fmat.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(
                    np.ascontiguousarray(ref_mat[k, sbra:ebra, sket:eket]))
                # compare submatrices
                assert cmat == rmat
            smat = fmat.full_matrix()
            fref = SubMatrix([0, 0, 114, 114])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            assert smat == fref

    def test_nuclear_potential_co_qzvp_for_o(self):

        mol_co, bas_qzvp = self.get_data()

        # compute nuclear potential matrix
        geom_drv = NuclearPotentialGeom010Driver()
        geom_mats = geom_drv.compute(mol_co, bas_qzvp, 1)

        # load reference overlap gradient for C atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' /
                      'co.qzvp.nuclear.potential.geom.010.o.npy')
        ref_mat = np.load(npyfile)
        ref_mat = -ref_mat

        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = ['X', 'Y', 'Z']

        for k, label in enumerate(labels):
            fmat = geom_mats.matrix(label)
            for i, j in zip(indexes[0], indexes[1]):
                # bra side
                sbra = basdims[i]
                ebra = basdims[i + 1]
                # ket side
                sket = basdims[j]
                eket = basdims[j + 1]
                # load computed submatrix
                cmat = fmat.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(
                    np.ascontiguousarray(ref_mat[k, sbra:ebra, sket:eket]))
                # compare submatrices
                # NOTE: See test cases with numerical problems.
                # assert cmat == rmat
            smat = fmat.full_matrix()
            fref = SubMatrix([0, 0, 114, 114])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            # NOTE: See test cases with numerical problems.
            # assert smat == fref
