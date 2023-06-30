from pathlib import Path
import numpy as np
import h5py

from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import Molecule


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
        bas = MolecularBasis.read(mol, 'sto-3g', basis_path, ostream=None)

        return mol, bas

    def test_overlap_of2_svp(self):

        mol_h2o, bas_svp = self.get_data()

        ovl_drv = OverlapDriver()
        ovl_mat = ovl_drv.compute(bas_svp, mol_h2o)

        fullmat = ovl_mat.get_full_matrix().to_numpy()

        ssmat = ovl_mat.get_submatrix((0, 0)).to_numpy()
        spmat = ovl_mat.get_submatrix((0, 1)).to_numpy()
        ppmat = ovl_mat.get_submatrix((1, 1)).to_numpy()

        nao_s = ssmat.shape[0]
        nao_p = ppmat.shape[0]

        assert ssmat.shape == (nao_s, nao_s)
        assert spmat.shape == (nao_s, nao_p)
        assert ppmat.shape == (nao_p, nao_p)

        assembled_mat = np.zeros((nao_s + nao_p, nao_s + nao_p))
        assembled_mat[:nao_s, :nao_s] = ssmat[:, :]
        assembled_mat[:nao_s, nao_s:] = spmat[:, :]
        assembled_mat[nao_s:, :nao_s] = spmat.T[:, :]
        assembled_mat[nao_s:, nao_s:] = ppmat[:, :]

        here = Path(__file__).parent
        h5file = str(here / 'data' / 'of2_overlap_sto3g.h5')
        hf = h5py.File(h5file, 'r')
        ref_mat = np.array(hf.get('of2_overlap_sto3g'))
        hf.close()

        assert np.max(np.abs(ref_mat - fullmat)) < 1.0e-12
        assert np.max(np.abs(ref_mat - assembled_mat)) < 1.0e-12
