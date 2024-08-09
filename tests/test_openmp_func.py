from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import set_number_of_threads
from veloxchem import get_number_of_threads
from veloxchem import make_gto_blocks
from veloxchem import make_work_tasks


class TestOpenMPFunc:
    """
    Implements tests for src/general/OpenMPFunc.hpp
    """

    def get_data(self):

        h2ostr = """O   0.000   0.000  -1.000
                    H   0.000   1.400  -2.100
                    H   0.000  -1.400  -2.100"""

        mol = Molecule.read_str(h2ostr, 'au')
        bas = MolecularBasis.read(mol, 'DEF2-SVP')

        return (mol, bas)

    def test_number_of_threads(self):

        nthreads = get_number_of_threads()

        set_number_of_threads(128)
        assert get_number_of_threads() == 128

        set_number_of_threads(nthreads)

    def test_make_work_tasks_for_gto_blocks(self):

        nthreads = get_number_of_threads()
        mol_h2o, bas_svp = self.get_data()
        set_number_of_threads(3)
        gblocks = make_gto_blocks(bas_svp, mol_h2o)
        wtasks = make_work_tasks(gblocks)
        rtasks = [[0, 0, 0, 4, 0, 4], [0, 1, 0, 4, 0, 2], [0, 2, 0, 4, 0, 1],
                  [0, 3, 0, 4, 0, 3], [0, 4, 0, 4, 0, 1], [0, 5, 0, 4, 0, 1],
                  [1, 1, 0, 2, 0, 2], [1, 2, 0, 2, 0, 1], [1, 3, 0, 2, 0, 3],
                  [1, 4, 0, 2, 0, 1], [1, 5, 0, 2, 0, 1], [2, 2, 0, 1, 0, 1],
                  [2, 3, 0, 1, 0, 3], [2, 4, 0, 1, 0, 1], [2, 5, 0, 1, 0, 1],
                  [3, 3, 0, 3, 0, 3], [3, 4, 0, 3, 0, 1], [3, 5, 0, 3, 0, 1],
                  [4, 4, 0, 1, 0, 1], [4, 5, 0, 1, 0, 1], [5, 5, 0, 1, 0, 1]]
        assert wtasks == rtasks

        wtasks = make_work_tasks(gblocks, gblocks)
        print(wtasks)
        rtasks = [[0, 0, 0, 4, 0, 4], [0, 1, 0, 4, 0, 2], [0, 2, 0, 4, 0, 1],
                  [0, 3, 0, 4, 0, 3], [0, 4, 0, 4, 0, 1], [0, 5, 0, 4, 0, 1],
                  [1, 0, 0, 2, 0, 4], [1, 1, 0, 2, 0, 2], [1, 2, 0, 2, 0, 1],
                  [1, 3, 0, 2, 0, 3], [1, 4, 0, 2, 0, 1], [1, 5, 0, 2, 0, 1],
                  [2, 0, 0, 1, 0, 4], [2, 1, 0, 1, 0, 2], [2, 2, 0, 1, 0, 1],
                  [2, 3, 0, 1, 0, 3], [2, 4, 0, 1, 0, 1], [2, 5, 0, 1, 0, 1],
                  [3, 0, 0, 3, 0, 4], [3, 1, 0, 3, 0, 2], [3, 2, 0, 3, 0, 1],
                  [3, 3, 0, 3, 0, 3], [3, 4, 0, 3, 0, 1], [3, 5, 0, 3, 0, 1],
                  [4, 0, 0, 1, 0, 4], [4, 1, 0, 1, 0, 2], [4, 2, 0, 1, 0, 1],
                  [4, 3, 0, 1, 0, 3], [4, 4, 0, 1, 0, 1], [4, 5, 0, 1, 0, 1],
                  [5, 0, 0, 1, 0, 4], [5, 1, 0, 1, 0, 2], [5, 2, 0, 1, 0, 1],
                  [5, 3, 0, 1, 0, 3], [5, 4, 0, 1, 0, 1], [5, 5, 0, 1, 0, 1]]
        assert wtasks == rtasks

        set_number_of_threads(nthreads)
