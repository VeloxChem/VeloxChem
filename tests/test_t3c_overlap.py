from pathlib import Path
import numpy as np

from veloxchem import ThreeCenterOverlapDriver
from veloxchem import MolecularBasis
from veloxchem import AtomBasis
from veloxchem import BasisFunction
from veloxchem import Molecule
from veloxchem import SubMatrix


class TestThreeCenterOverlapDriver:

    def get_data(self):

        h2str = """
         H   0.0   0.0   0.0
         H   0.3   0.2   1.0   
        """
        mol = Molecule.read_str(h2str, 'au')
        
        h1_1s = BasisFunction([2.0000, ], [1.0000, ], 0)
        h1_1s.normalize()
        h2_1s = BasisFunction([3.2000, ], [1.0000, ], 0)
        h2_1s.normalize()
        
        h1_bas = AtomBasis([h1_1s, ], 'MINI-A', '', 1)
        h2_bas = AtomBasis([h2_1s, ], 'MINI-B', '', 1)
        bas = MolecularBasis([h1_bas, h2_bas], [0, 1])
        
        return mol, bas

    def test_overlap_h2_minimal(self):

        mol_h2, bas_min = self.get_data()

        # compute overlap matrix
        ovl_drv = ThreeCenterOverlapDriver()
        ovl_mat = ovl_drv.compute(mol_h2, bas_min, [0.8, ], [1.0, ], [[1.1, 1.2, 1.3], ])
        
        # check full overlap matrix
        fmat = ovl_mat.full_matrix()
        fref = SubMatrix([0, 0, 2, 2])
        ref_mat = np.array([[0.042137564466, 0.034857221722],
                            [0.034857221722, 0.244902091288]])
        fref.set_values(np.ascontiguousarray(ref_mat))
        assert fmat == fref
