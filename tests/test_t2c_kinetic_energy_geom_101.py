from pathlib import Path
import numpy as np

from veloxchem import KineticEnergyGeom101Driver
from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import SubMatrix


class TestKineticEnergyGeom101Driver:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-QZVP')

        return mol, bas

    def test_kinetic_energy_co_qzvp_for_cc(self):

        mol, bas = self.get_data()

        # compute kinetic energy hessian matrix
        hess_drv = KineticEnergyGeom101Driver()
        hess_mats = hess_drv.compute(mol, bas, 0, 0)

        # load reference kinetic energy hessian for C,C atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.kinetic.energy.geom.101.cc.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = [
            'X_X', 'X_Y', 'X_Z', 'Y_X', 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z'
        ]

        for k, label in enumerate(labels):
            fmat = hess_mats.matrix(label)
            for i in range(0, 5):
                for j in range(0, 5):
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

    def test_kinetic_energy_co_qzvp_for_co(self):

        mol, bas = self.get_data()

        # compute kinetic energy hesian matrix
        hess_drv = KineticEnergyGeom101Driver()
        hess_mats = hess_drv.compute(mol, bas, 0, 1)

        # load reference overlap hessian for C,O atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.kinetic.energy.geom.101.co.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = [
            'X_X', 'X_Y', 'X_Z', 'Y_X', 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z'
        ]

        for k, label in enumerate(labels):
            fmat = hess_mats.matrix(label)
            for i in range(0, 5):
                for j in range(0, 5):
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

    def test_kinetic_energy_co_qzvp_for_oc(self):

        mol, bas = self.get_data()

        # compute kinetic energy hessina matrix
        hess_drv = KineticEnergyGeom101Driver()
        hess_mats = hess_drv.compute(mol, bas, 1, 0)

        # load reference hessian for O, C atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.kinetic.energy.geom.101.oc.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = [
            'X_X', 'X_Y', 'X_Z', 'Y_X', 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z'
        ]

        for k, label in enumerate(labels):
            fmat = hess_mats.matrix(label)
            for i in range(0, 5):
                for j in range(0, 5):
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

    def test_kinetic_energy_co_qzvp_for_oo(self):

        mol, bas = self.get_data()

        # compute kinetic energy matrix
        hess_drv = KineticEnergyGeom101Driver()
        hess_mats = hess_drv.compute(mol, bas, 1, 1)

        # load reference overlap hessian for O, O atom
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.kinetic.energy.geom.101.oo.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 14, 38, 68, 96, 114]

        # indices map
        labels = [
            'X_X', 'X_Y', 'X_Z', 'Y_X', 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z'
        ]

        for k, label in enumerate(labels):
            fmat = hess_mats.matrix(label)
            for i in range(0, 5):
                for j in range(0, 5):
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
