from pathlib import Path
import numpy as np
import h5py

from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import Molecule
from veloxchem.submatrix import SubMatrix
from tester import Tester


class TestOverlapDriver:

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

    def test_overlap_of2_qzvpp(self):

        mol_h2o, bas_qzvpp = self.get_data()

        ovl_drv = OverlapDriver()
        ovl_mat = ovl_drv.compute(bas_qzvpp, mol_h2o)

        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'of2.qzvpp.overlap.npy')
        ref_mat = np.load(npyfile)

        # (s|s) integrals
        ssmat = ovl_mat.get_submatrix((0, 0))
        ssref = SubMatrix([0, 0, 21, 21])
        ssref.set_values(np.ascontiguousarray(ref_mat[:21, :21]))
        Tester.compare_submatrices(ssmat, ssref)

        # (s|p) integrals
        spmat = ovl_mat.get_submatrix((0, 1))
        spref = SubMatrix([0, 21, 21, 36])
        spref.set_values(np.ascontiguousarray(ref_mat[:21, 21:57]))
        Tester.compare_submatrices(spmat, spref)

        # (s|d) integrals
        sdmat = ovl_mat.get_submatrix((0, 2))
        sdref = SubMatrix([0, 57, 21, 45])
        sdref.set_values(np.ascontiguousarray(ref_mat[:21, 57:102]))
        Tester.compare_submatrices(sdmat, sdref)

        # (s|f) integrals
        sfmat = ovl_mat.get_submatrix((0, 3))
        sfref = SubMatrix([0, 102, 21, 42])
        sfref.set_values(np.ascontiguousarray(ref_mat[:21, 102:144]))
        Tester.compare_submatrices(sfmat, sfref)

        # (s|g) integrals
        sgmat = ovl_mat.get_submatrix((0, 4))
        sgref = SubMatrix([0, 144, 21, 27])
        sgref.set_values(np.ascontiguousarray(ref_mat[:21, 144:171]))
        Tester.compare_submatrices(sgmat, sgref)

        # (p|p) integrals
        ppmat = ovl_mat.get_submatrix((1, 1))
        ppref = SubMatrix([21, 21, 36, 36])
        ppref.set_values(np.ascontiguousarray(ref_mat[21:57, 21:57]))
        Tester.compare_submatrices(ppmat, ppref)

        # (p|d) integrals
        pdmat = ovl_mat.get_submatrix((1, 2))
        pdref = SubMatrix([21, 57, 36, 45])
        pdref.set_values(np.ascontiguousarray(ref_mat[21:57, 57:102]))
        Tester.compare_submatrices(pdmat, pdref)

        # (p|f) integrals
        pfmat = ovl_mat.get_submatrix((1, 3))
        pfref = SubMatrix([21, 102, 36, 42])
        pfref.set_values(np.ascontiguousarray(ref_mat[21:57, 102:144]))
        Tester.compare_submatrices(pfmat, pfref)

        # (p|g) integrals
        pgmat = ovl_mat.get_submatrix((1, 4))
        pgref = SubMatrix([21, 144, 36, 27])
        pgref.set_values(np.ascontiguousarray(ref_mat[21:57, 144:171]))
        Tester.compare_submatrices(pgmat, pgref)

        # (d|d) integrals
        ddmat = ovl_mat.get_submatrix((2, 2))
        ddref = SubMatrix([57, 57, 45, 45])
        ddref.set_values(np.ascontiguousarray(ref_mat[57:102, 57:102]))
        Tester.compare_submatrices(ddmat, ddref)

        # (d|f) integrals
        dfmat = ovl_mat.get_submatrix((2, 3))
        dfref = SubMatrix([57, 102, 45, 42])
        dfref.set_values(np.ascontiguousarray(ref_mat[57:102, 102:144]))
        Tester.compare_submatrices(dfmat, dfref)

        # (d|g) integrals
        dgmat = ovl_mat.get_submatrix((2, 4))
        dgref = SubMatrix([57, 144, 45, 27])
        dgref.set_values(np.ascontiguousarray(ref_mat[57:102, 144:171]))
        Tester.compare_submatrices(dgmat, dgref)

        # (f|f) integrals
        ffmat = ovl_mat.get_submatrix((3, 3))
        ffref = SubMatrix([102, 102, 42, 42])
        ffref.set_values(np.ascontiguousarray(ref_mat[102:144, 102:144]))
        Tester.compare_submatrices(ffmat, ffref)

        # (f|g) integrals
        fgmat = ovl_mat.get_submatrix((3, 4))
        fgref = SubMatrix([102, 144, 42, 27])
        fgref.set_values(np.ascontiguousarray(ref_mat[102:144, 144:171]))
        Tester.compare_submatrices(fgmat, fgref)

        # (g|g) integrals
        ggmat = ovl_mat.get_submatrix((4, 4))
        ggref = SubMatrix([144, 144, 27, 27])
        ggref.set_values(np.ascontiguousarray(ref_mat[144:171, 144:171]))
        Tester.compare_submatrices(ggmat, ggref)

        # full matrix
        full_mat = ovl_mat.get_full_matrix()
        full_ref = SubMatrix([0, 0, 171, 171])
        full_ref.set_values(np.ascontiguousarray(ref_mat))
        Tester.compare_submatrices(full_mat, full_ref)
