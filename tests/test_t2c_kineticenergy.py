from pathlib import Path
import numpy as np
import h5py

from veloxchem.veloxchemlib import KineticEnergyDriver
from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import Molecule
from veloxchem.submatrix import SubMatrix
from tester import Tester


class TestKineticEnergyDriver:

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

    def test_kinetic_energy_of2_qzvpp(self):

        mol_h2o, bas_qzvpp = self.get_data()

        kin_drv = KineticEnergyDriver()
        kin_mat = kin_drv.compute(bas_qzvpp, mol_h2o)

        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'of2.qzvpp.kinetic.energy.npy')
        ref_mat = np.load(npyfile)

        # (s|T|s) integrals
        ssmat = kin_mat.get_submatrix((0, 0))
        ssref = SubMatrix([0, 0, 21, 21])
        ssref.set_values(np.ascontiguousarray(ref_mat[:21, :21]))
        Tester.compare_submatrices(ssmat, ssref)

        # (s|T|p) integrals
        spmat = kin_mat.get_submatrix((0, 1))
        spref = SubMatrix([0, 21, 21, 36])
        spref.set_values(np.ascontiguousarray(ref_mat[:21, 21:57]))
        Tester.compare_submatrices(spmat, spref)

        # (s|T|d) integrals
        sdmat = kin_mat.get_submatrix((0, 2))
        sdref = SubMatrix([0, 57, 21, 45])
        sdref.set_values(np.ascontiguousarray(ref_mat[:21, 57:102]))
        Tester.compare_submatrices(sdmat, sdref)

        # (s|T|f) integrals
        sfmat = kin_mat.get_submatrix((0, 3))
        sfref = SubMatrix([0, 102, 21, 42])
        sfref.set_values(np.ascontiguousarray(ref_mat[:21, 102:144]))
        Tester.compare_submatrices(sfmat, sfref)

        # (s|T|g) integrals
        sgmat = kin_mat.get_submatrix((0, 4))
        sgref = SubMatrix([0, 144, 21, 27])
        sgref.set_values(np.ascontiguousarray(ref_mat[:21, 144:171]))
        Tester.compare_submatrices(sgmat, sgref)

        # (p|T|p) integrals
        ppmat = kin_mat.get_submatrix((1, 1))
        ppref = SubMatrix([21, 21, 36, 36])
        ppref.set_values(np.ascontiguousarray(ref_mat[21:57, 21:57]))
        Tester.compare_submatrices(ppmat, ppref)

        # (p|T|d) integrals
        pdmat = kin_mat.get_submatrix((1, 2))
        pdref = SubMatrix([21, 57, 36, 45])
        pdref.set_values(np.ascontiguousarray(ref_mat[21:57, 57:102]))
        Tester.compare_submatrices(pdmat, pdref)

        # (p|T|f) integrals
        pfmat = kin_mat.get_submatrix((1, 3))
        pfref = SubMatrix([21, 102, 36, 42])
        pfref.set_values(np.ascontiguousarray(ref_mat[21:57, 102:144]))
        Tester.compare_submatrices(pfmat, pfref)

        # (p|T|g) integrals
        pgmat = kin_mat.get_submatrix((1, 4))
        pgref = SubMatrix([21, 144, 36, 27])
        pgref.set_values(np.ascontiguousarray(ref_mat[21:57, 144:171]))
        Tester.compare_submatrices(pgmat, pgref)

        # (d|T|d) integrals
        ddmat = kin_mat.get_submatrix((2, 2))
        ddref = SubMatrix([57, 57, 45, 45])
        ddref.set_values(np.ascontiguousarray(ref_mat[57:102, 57:102]))
        Tester.compare_submatrices(ddmat, ddref)

        # (d|T|f) integrals
        dfmat = kin_mat.get_submatrix((2, 3))
        dfref = SubMatrix([57, 102, 45, 42])
        dfref.set_values(np.ascontiguousarray(ref_mat[57:102, 102:144]))
        Tester.compare_submatrices(dfmat, dfref)

        # (d|T|g) integrals
        dgmat = kin_mat.get_submatrix((2, 4))
        dgref = SubMatrix([57, 144, 45, 27])
        dgref.set_values(np.ascontiguousarray(ref_mat[57:102, 144:171]))
        Tester.compare_submatrices(dgmat, dgref)

        # (f|T|f) integrals
        ffmat = kin_mat.get_submatrix((3, 3))
        ffref = SubMatrix([102, 102, 42, 42])
        ffref.set_values(np.ascontiguousarray(ref_mat[102:144, 102:144]))
        Tester.compare_submatrices(ffmat, ffref)

        # (f|T|g) integrals
        fgmat = kin_mat.get_submatrix((3, 4))
        fgref = SubMatrix([102, 144, 42, 27])
        fgref.set_values(np.ascontiguousarray(ref_mat[102:144, 144:171]))
        Tester.compare_submatrices(fgmat, fgref)

        # (g|T|g) integrals
        ggmat = kin_mat.get_submatrix((4, 4))
        ggref = SubMatrix([144, 144, 27, 27])
        ggref.set_values(np.ascontiguousarray(ref_mat[144:171, 144:171]))
        Tester.compare_submatrices(ggmat, ggref)

        # full matrix
        full_mat = kin_mat.get_full_matrix()
        full_ref = SubMatrix([0, 0, 171, 171])
        full_ref.set_values(np.ascontiguousarray(ref_mat))
        Tester.compare_submatrices(full_mat, full_ref)
