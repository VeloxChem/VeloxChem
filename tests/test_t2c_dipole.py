from pathlib import Path
import numpy as np
import h5py

from veloxchem.veloxchemlib import DipoleDriver
from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import Molecule
from veloxchem.veloxchemlib import Matrices
from veloxchem.submatrix import SubMatrix
from veloxchem.matrix import Matrix
from tester import Tester


class TestDipoleDriver:

    def get_data(self):

        of2str = """
            O   0.100  -0.400  -1.000
            F   0.300   1.400  -2.100
            F   0.200  -1.300  -2.000
        """
        mol = Molecule.read_str(of2str, 'au')

        here = Path(__file__).parent
        basis_path = str(here.parent / 'basis')
        bas = MolecularBasis.read(mol, 'def2-qzvpp', basis_path, ostream=None)

        return mol, bas

    def test_dipole_of2_qzvpp(self):

        mol_h2o, bas_qzvpp = self.get_data()

        dip_drv = DipoleDriver()
        dip_mats = dip_drv.compute(bas_qzvpp, mol_h2o, [0.0, 0.0, 0.0])

        dipx_mat = dip_mats.get_matrix('x')
        dipy_mat = dip_mats.get_matrix('y')
        dipz_mat = dip_mats.get_matrix('z')

        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'of2.qzvpp.electric.dipole.npy')
        ref_mats = np.load(npyfile)

        # (s|r|s) integrals
        ssmat_x = dipx_mat.get_submatrix((0, 0))
        ssmat_y = dipy_mat.get_submatrix((0, 0))
        ssmat_z = dipz_mat.get_submatrix((0, 0))
        ssref_x = SubMatrix([0, 0, 21, 21])
        ssref_x.set_values(np.ascontiguousarray(ref_mats[0, :21, :21]))
        ssref_y = SubMatrix([0, 0, 21, 21])
        ssref_y.set_values(np.ascontiguousarray(ref_mats[1, :21, :21]))
        ssref_z = SubMatrix([0, 0, 21, 21])
        ssref_z.set_values(np.ascontiguousarray(ref_mats[2, :21, :21]))
        Tester.compare_submatrices(ssmat_x, ssref_x)
        Tester.compare_submatrices(ssmat_y, ssref_y)
        Tester.compare_submatrices(ssmat_z, ssref_z)

        # (s|r|p) integrals
        spmat_x = dipx_mat.get_submatrix((0, 1))
        spmat_y = dipy_mat.get_submatrix((0, 1))
        spmat_z = dipz_mat.get_submatrix((0, 1))
        spref_x = SubMatrix([0, 21, 21, 36])
        spref_x.set_values(np.ascontiguousarray(ref_mats[0, :21, 21:57]))
        spref_y = SubMatrix([0, 21, 21, 36])
        spref_y.set_values(np.ascontiguousarray(ref_mats[1, :21, 21:57]))
        spref_z = SubMatrix([0, 21, 21, 36])
        spref_z.set_values(np.ascontiguousarray(ref_mats[2, :21, 21:57]))
        Tester.compare_submatrices(spmat_x, spref_x)
        Tester.compare_submatrices(spmat_y, spref_y)
        Tester.compare_submatrices(spmat_z, spref_z)

        # (s|r|d) integrals
        sdmat_x = dipx_mat.get_submatrix((0, 2))
        sdmat_y = dipy_mat.get_submatrix((0, 2))
        sdmat_z = dipz_mat.get_submatrix((0, 2))
        sdref_x = SubMatrix([0, 57, 21, 45])
        sdref_x.set_values(np.ascontiguousarray(ref_mats[0, :21, 57:102]))
        sdref_y = SubMatrix([0, 57, 21, 45])
        sdref_y.set_values(np.ascontiguousarray(ref_mats[1, :21, 57:102]))
        sdref_z = SubMatrix([0, 57, 21, 45])
        sdref_z.set_values(np.ascontiguousarray(ref_mats[2, :21, 57:102]))
        Tester.compare_submatrices(sdmat_x, sdref_x)
        Tester.compare_submatrices(sdmat_y, sdref_y)
        Tester.compare_submatrices(sdmat_z, sdref_z)

        # (s|r|f) integrals
        sfmat_x = dipx_mat.get_submatrix((0, 3))
        sfmat_y = dipy_mat.get_submatrix((0, 3))
        sfmat_z = dipz_mat.get_submatrix((0, 3))
        sfref_x = SubMatrix([0, 102, 21, 42])
        sfref_x.set_values(np.ascontiguousarray(ref_mats[0, :21, 102:144]))
        sfref_y = SubMatrix([0, 102, 21, 42])
        sfref_y.set_values(np.ascontiguousarray(ref_mats[1, :21, 102:144]))
        sfref_z = SubMatrix([0, 102, 21, 42])
        sfref_z.set_values(np.ascontiguousarray(ref_mats[2, :21, 102:144]))
        Tester.compare_submatrices(sfmat_x, sfref_x)
        Tester.compare_submatrices(sfmat_y, sfref_y)
        Tester.compare_submatrices(sfmat_z, sfref_z)

        # (s|r|g) integrals
        sgmat_x = dipx_mat.get_submatrix((0, 4))
        sgmat_y = dipy_mat.get_submatrix((0, 4))
        sgmat_z = dipz_mat.get_submatrix((0, 4))
        sgref_x = SubMatrix([0, 144, 21, 27])
        sgref_x.set_values(np.ascontiguousarray(ref_mats[0, :21, 144:171]))
        sgref_y = SubMatrix([0, 144, 21, 27])
        sgref_y.set_values(np.ascontiguousarray(ref_mats[1, :21, 144:171]))
        sgref_z = SubMatrix([0, 144, 21, 27])
        sgref_z.set_values(np.ascontiguousarray(ref_mats[2, :21, 144:171]))
        Tester.compare_submatrices(sgmat_x, sgref_x)
        Tester.compare_submatrices(sgmat_y, sgref_y)
        Tester.compare_submatrices(sgmat_z, sgref_z)

        # (p|r|p) integrals
        ppmat_x = dipx_mat.get_submatrix((1, 1))
        ppmat_y = dipy_mat.get_submatrix((1, 1))
        ppmat_z = dipz_mat.get_submatrix((1, 1))
        ppref_x = SubMatrix([21, 21, 36, 36])
        ppref_x.set_values(np.ascontiguousarray(ref_mats[0, 21:57, 21:57]))
        ppref_y = SubMatrix([21, 21, 36, 36])
        ppref_y.set_values(np.ascontiguousarray(ref_mats[1, 21:57, 21:57]))
        ppref_z = SubMatrix([21, 21, 36, 36])
        ppref_z.set_values(np.ascontiguousarray(ref_mats[2, 21:57, 21:57]))
        Tester.compare_submatrices(ppmat_x, ppref_x)
        Tester.compare_submatrices(ppmat_y, ppref_y)
        Tester.compare_submatrices(ppmat_z, ppref_z)

        # (p|r|d) integrals
        pdmat_x = dipx_mat.get_submatrix((1, 2))
        pdmat_y = dipy_mat.get_submatrix((1, 2))
        pdmat_z = dipz_mat.get_submatrix((1, 2))
        pdref_x = SubMatrix([21, 57, 36, 45])
        pdref_x.set_values(np.ascontiguousarray(ref_mats[0, 21:57, 57:102]))
        pdref_y = SubMatrix([21, 57, 36, 45])
        pdref_y.set_values(np.ascontiguousarray(ref_mats[1, 21:57, 57:102]))
        pdref_z = SubMatrix([21, 57, 36, 45])
        pdref_z.set_values(np.ascontiguousarray(ref_mats[2, 21:57, 57:102]))
        Tester.compare_submatrices(pdmat_x, pdref_x)
        Tester.compare_submatrices(pdmat_y, pdref_y)
        Tester.compare_submatrices(pdmat_z, pdref_z)

        # (p|r|f) integrals
        pfmat_x = dipx_mat.get_submatrix((1, 3))
        pfmat_y = dipy_mat.get_submatrix((1, 3))
        pfmat_z = dipz_mat.get_submatrix((1, 3))
        pfref_x = SubMatrix([21, 102, 36, 42])
        pfref_x.set_values(np.ascontiguousarray(ref_mats[0, 21:57, 102:144]))
        pfref_y = SubMatrix([21, 102, 36, 42])
        pfref_y.set_values(np.ascontiguousarray(ref_mats[1, 21:57, 102:144]))
        pfref_z = SubMatrix([21, 102, 36, 42])
        pfref_z.set_values(np.ascontiguousarray(ref_mats[2, 21:57, 102:144]))
        Tester.compare_submatrices(pfmat_x, pfref_x)
        Tester.compare_submatrices(pfmat_y, pfref_y)
        Tester.compare_submatrices(pfmat_z, pfref_z)

        # (p|r|g) integrals
        pgmat_x = dipx_mat.get_submatrix((1, 4))
        pgmat_y = dipy_mat.get_submatrix((1, 4))
        pgmat_z = dipz_mat.get_submatrix((1, 4))
        pgref_x = SubMatrix([21, 144, 36, 27])
        pgref_x.set_values(np.ascontiguousarray(ref_mats[0, 21:57, 144:171]))
        pgref_y = SubMatrix([21, 144, 36, 27])
        pgref_y.set_values(np.ascontiguousarray(ref_mats[1, 21:57, 144:171]))
        pgref_z = SubMatrix([21, 144, 36, 27])
        pgref_z.set_values(np.ascontiguousarray(ref_mats[2, 21:57, 144:171]))
        Tester.compare_submatrices(pgmat_x, pgref_x)
        Tester.compare_submatrices(pgmat_y, pgref_y)
        Tester.compare_submatrices(pgmat_z, pgref_z)

        # (d|r|d) integrals
        ddmat_x = dipx_mat.get_submatrix((2, 2))
        ddmat_y = dipy_mat.get_submatrix((2, 2))
        ddmat_z = dipz_mat.get_submatrix((2, 2))
        ddref_x = SubMatrix([57, 57, 45, 45])
        ddref_x.set_values(np.ascontiguousarray(ref_mats[0, 57:102, 57:102]))
        ddref_y = SubMatrix([57, 57, 45, 45])
        ddref_y.set_values(np.ascontiguousarray(ref_mats[1, 57:102, 57:102]))
        ddref_z = SubMatrix([57, 57, 45, 45])
        ddref_z.set_values(np.ascontiguousarray(ref_mats[2, 57:102, 57:102]))
        Tester.compare_submatrices(ddmat_x, ddref_x)
        Tester.compare_submatrices(ddmat_y, ddref_y)
        Tester.compare_submatrices(ddmat_z, ddref_z)

        # (d|r|f) integrals
        dfmat_x = dipx_mat.get_submatrix((2, 3))
        dfmat_y = dipy_mat.get_submatrix((2, 3))
        dfmat_z = dipz_mat.get_submatrix((2, 3))
        dfref_x = SubMatrix([57, 102, 45, 42])
        dfref_x.set_values(np.ascontiguousarray(ref_mats[0, 57:102, 102:144]))
        dfref_y = SubMatrix([57, 102, 45, 42])
        dfref_y.set_values(np.ascontiguousarray(ref_mats[1, 57:102, 102:144]))
        dfref_z = SubMatrix([57, 102, 45, 42])
        dfref_z.set_values(np.ascontiguousarray(ref_mats[2, 57:102, 102:144]))
        Tester.compare_submatrices(dfmat_x, dfref_x)
        Tester.compare_submatrices(dfmat_y, dfref_y)
        Tester.compare_submatrices(dfmat_z, dfref_z)

        # (d|r|g) integrals
        dgmat_x = dipx_mat.get_submatrix((2, 4))
        dgmat_y = dipy_mat.get_submatrix((2, 4))
        dgmat_z = dipz_mat.get_submatrix((2, 4))
        dgref_x = SubMatrix([57, 144, 45, 27])
        dgref_x.set_values(np.ascontiguousarray(ref_mats[0, 57:102, 144:171]))
        dgref_y = SubMatrix([57, 144, 45, 27])
        dgref_y.set_values(np.ascontiguousarray(ref_mats[1, 57:102, 144:171]))
        dgref_z = SubMatrix([57, 144, 45, 27])
        dgref_z.set_values(np.ascontiguousarray(ref_mats[2, 57:102, 144:171]))
        Tester.compare_submatrices(dgmat_x, dgref_x)
        Tester.compare_submatrices(dgmat_y, dgref_y)
        Tester.compare_submatrices(dgmat_z, dgref_z)

        # (f|r|f) integrals
        ffmat_x = dipx_mat.get_submatrix((3, 3))
        ffmat_y = dipy_mat.get_submatrix((3, 3))
        ffmat_z = dipz_mat.get_submatrix((3, 3))
        ffref_x = SubMatrix([102, 102, 42, 42])
        ffref_x.set_values(np.ascontiguousarray(ref_mats[0, 102:144, 102:144]))
        ffref_y = SubMatrix([102, 102, 42, 42])
        ffref_y.set_values(np.ascontiguousarray(ref_mats[1, 102:144, 102:144]))
        ffref_z = SubMatrix([102, 102, 42, 42])
        ffref_z.set_values(np.ascontiguousarray(ref_mats[2, 102:144, 102:144]))
        Tester.compare_submatrices(ffmat_x, ffref_x)
        Tester.compare_submatrices(ffmat_y, ffref_y)
        Tester.compare_submatrices(ffmat_z, ffref_z)

        # (f|r|g) integrals
        fgmat_x = dipx_mat.get_submatrix((3, 4))
        fgmat_y = dipy_mat.get_submatrix((3, 4))
        fgmat_z = dipz_mat.get_submatrix((3, 4))
        fgref_x = SubMatrix([102, 144, 42, 27])
        fgref_x.set_values(np.ascontiguousarray(ref_mats[0, 102:144, 144:171]))
        fgref_y = SubMatrix([102, 144, 42, 27])
        fgref_y.set_values(np.ascontiguousarray(ref_mats[1, 102:144, 144:171]))
        fgref_z = SubMatrix([102, 144, 42, 27])
        fgref_z.set_values(np.ascontiguousarray(ref_mats[2, 102:144, 144:171]))
        Tester.compare_submatrices(fgmat_x, fgref_x)
        Tester.compare_submatrices(fgmat_y, fgref_y)
        Tester.compare_submatrices(fgmat_z, fgref_z)

        # (g|r|g) integrals
        ggmat_x = dipx_mat.get_submatrix((4, 4))
        ggmat_y = dipy_mat.get_submatrix((4, 4))
        ggmat_z = dipz_mat.get_submatrix((4, 4))
        ggref_x = SubMatrix([144, 144, 27, 27])
        ggref_x.set_values(np.ascontiguousarray(ref_mats[0, 144:171, 144:171]))
        ggref_y = SubMatrix([144, 144, 27, 27])
        ggref_y.set_values(np.ascontiguousarray(ref_mats[1, 144:171, 144:171]))
        ggref_z = SubMatrix([144, 144, 27, 27])
        ggref_z.set_values(np.ascontiguousarray(ref_mats[2, 144:171, 144:171]))
        Tester.compare_submatrices(ggmat_x, ggref_x)
        Tester.compare_submatrices(ggmat_y, ggref_y)
        Tester.compare_submatrices(ggmat_z, ggref_z)

        # full matrix
        fullmat_x = dipx_mat.get_full_matrix()
        fullmat_y = dipy_mat.get_full_matrix()
        fullmat_z = dipz_mat.get_full_matrix()
        fullref_x = SubMatrix([0, 0, 171, 171])
        fullref_x.set_values(np.ascontiguousarray(ref_mats[0, :171, :171]))
        fullref_y = SubMatrix([0, 0, 171, 171])
        fullref_y.set_values(np.ascontiguousarray(ref_mats[1, :171, :171]))
        fullref_z = SubMatrix([0, 0, 171, 171])
        fullref_z.set_values(np.ascontiguousarray(ref_mats[2, :171, :171]))
        Tester.compare_submatrices(fullmat_x, fullref_x)
        Tester.compare_submatrices(fullmat_y, fullref_y)
        Tester.compare_submatrices(fullmat_z, fullref_z)
