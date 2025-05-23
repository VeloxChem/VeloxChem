from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.veloxchemlib import set_number_of_threads
from veloxchem.veloxchemlib import get_number_of_threads
from veloxchem.veloxchemlib import make_gto_blocks
from veloxchem.veloxchemlib import make_gto_pair_blocks
from veloxchem.veloxchemlib import make_work_tasks
from veloxchem.veloxchemlib import make_diag_work_tasks
from veloxchem.veloxchemlib import make_work_group
from veloxchem.veloxchemlib import T4CScreener


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
        #print(wtasks)
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

    def test_make_diag_work_tasks_for_gto_pair_blocks(self):

        nthreads = get_number_of_threads()
        mol_h2o, bas_svp = self.get_data()
        set_number_of_threads(3)

        gpairs = make_gto_pair_blocks(bas_svp, mol_h2o)
        wtasks = make_diag_work_tasks(gpairs)
        rtasks = [[0, 0, 10], [1, 0, 8], [2, 0, 4], [3, 0, 12], [4, 0, 4],
                  [5, 0, 4], [6, 0, 3], [7, 0, 2], [8, 0, 6], [9, 0, 2],
                  [10, 0, 2], [11, 0, 1], [12, 0, 3], [13, 0, 1], [14, 0, 1],
                  [15, 0, 6], [16, 0, 3], [17, 0, 3], [18, 0, 1], [19, 0, 1],
                  [20, 0, 1]]
        assert wtasks == rtasks

        set_number_of_threads(nthreads)

    def test_make_diag_work_group(self):

        nthreads = get_number_of_threads()
        mol_h2o, bas_svp = self.get_data()
        set_number_of_threads(3)

        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svp, mol_h2o, 'eri')

        wtasks = make_work_group(t4c_drv.gto_pair_blocks(), 12, 1)
        rtasks = [[0, 0, 0, 0, 0, 1, 0, 1], [0, 0, 0, 1, 0, 1, 0, 9],
                  [0, 0, 1, 1, 0, 9, 0, 9], [0, 1, 0, 1, 0, 1, 0, 8],
                  [0, 1, 1, 1, 0, 9, 0, 8], [0, 2, 0, 1, 0, 1, 0, 4],
                  [0, 2, 1, 1, 0, 9, 0, 4], [0, 3, 0, 1, 0, 1, 0, 12],
                  [0, 3, 1, 1, 0, 9, 0, 12], [0, 4, 0, 1, 0, 1, 0, 4],
                  [0, 4, 1, 1, 0, 9, 0, 4], [0, 5, 0, 1, 0, 1, 0, 4],
                  [0, 5, 1, 1, 0, 9, 0, 4], [0, 6, 0, 1, 0, 1, 0, 3],
                  [0, 6, 1, 1, 0, 9, 0, 3], [0, 7, 0, 2, 0, 1, 0, 2],
                  [0, 7, 1, 2, 0, 9, 0, 2], [0, 8, 0, 1, 0, 1, 0, 6],
                  [0, 8, 1, 1, 0, 9, 0, 6], [0, 9, 0, 1, 0, 1, 0, 2],
                  [0, 9, 1, 1, 0, 9, 0, 2], [0, 10, 0, 1, 0, 1, 0, 2],
                  [0, 10, 1, 1, 0, 9, 0, 2], [0, 11, 0, 0, 0, 1, 0, 1],
                  [0, 11, 1, 0, 0, 9, 0, 1], [0, 12, 0, 1, 0, 1, 0, 2],
                  [0, 12, 0, 2, 0, 1, 0, 1], [0, 12, 1, 1, 0, 9, 0, 2],
                  [0, 12, 1, 2, 0, 9, 0, 1], [0, 13, 0, 1, 0, 1, 0, 1],
                  [0, 13, 1, 1, 0, 9, 0, 1], [0, 14, 0, 2, 0, 1, 0, 1],
                  [0, 14, 1, 2, 0, 9, 0, 1], [0, 15, 0, 1, 0, 1, 0, 6],
                  [0, 15, 1, 1, 0, 9, 0, 6], [0, 16, 0, 1, 0, 1, 0, 3],
                  [0, 16, 1, 1, 0, 9, 0, 3], [0, 17, 0, 1, 0, 1, 0, 3],
                  [0, 17, 1, 1, 0, 9, 0, 3], [0, 18, 0, 0, 0, 1, 0, 1],
                  [0, 18, 1, 0, 0, 9, 0, 1], [0, 19, 0, 1, 0, 1, 0, 1],
                  [0, 19, 1, 1, 0, 9, 0, 1], [0, 20, 0, 1, 0, 1, 0, 1],
                  [0, 20, 1, 1, 0, 9, 0, 1], [1, 1, 1, 1, 0, 8, 0, 8],
                  [1, 2, 1, 1, 0, 8, 0, 4], [1, 3, 1, 1, 0, 8, 0, 12],
                  [1, 4, 1, 1, 0, 8, 0, 4], [1, 5, 1, 1, 0, 8, 0, 4],
                  [1, 6, 1, 1, 0, 8, 0, 3], [1, 7, 1, 2, 0, 8, 0, 2],
                  [1, 8, 1, 1, 0, 8, 0, 6], [1, 9, 1, 1, 0, 8, 0, 2],
                  [1, 10, 1, 1, 0, 8, 0, 2], [1, 11, 1, 0, 0, 8, 0, 1],
                  [1, 12, 1, 1, 0, 8, 0, 2], [1, 12, 1, 2, 0, 8, 0, 1],
                  [1, 13, 1, 1, 0, 8, 0, 1], [1, 14, 1, 2, 0, 8, 0, 1],
                  [1, 15, 1, 1, 0, 8, 0, 6], [1, 16, 1, 1, 0, 8, 0, 3],
                  [1, 17, 1, 1, 0, 8, 0, 3], [1, 18, 1, 0, 0, 8, 0, 1],
                  [1, 19, 1, 1, 0, 8, 0, 1], [1, 20, 1, 1, 0, 8, 0, 1],
                  [2, 2, 1, 1, 0, 4, 0, 4], [2, 3, 1, 1, 0, 4, 0, 12],
                  [2, 4, 1, 1, 0, 4, 0, 4], [2, 5, 1, 1, 0, 4, 0, 4],
                  [2, 6, 1, 1, 0, 4, 0, 3], [2, 7, 1, 2, 0, 4, 0, 2],
                  [2, 8, 1, 1, 0, 4, 0, 6], [2, 9, 1, 1, 0, 4, 0, 2],
                  [2, 10, 1, 1, 0, 4, 0, 2], [2, 11, 1, 0, 0, 4, 0, 1],
                  [2, 12, 1, 1, 0, 4, 0, 2], [2, 12, 1, 2, 0, 4, 0, 1],
                  [2, 13, 1, 1, 0, 4, 0, 1], [2, 14, 1, 2, 0, 4, 0, 1],
                  [2, 15, 1, 1, 0, 4, 0, 6], [2, 16, 1, 1, 0, 4, 0, 3],
                  [2, 17, 1, 1, 0, 4, 0, 3], [2, 18, 1, 0, 0, 4, 0, 1],
                  [2, 19, 1, 1, 0, 4, 0, 1], [2, 20, 1, 1, 0, 4, 0, 1],
                  [3, 3, 1, 1, 0, 12, 0, 12], [3, 4, 1, 1, 0, 12, 0, 4],
                  [3, 5, 1, 1, 0, 12, 0, 4], [3, 6, 1, 1, 0, 12, 0, 3],
                  [3, 7, 1, 2, 0, 12, 0, 2], [3, 8, 1, 1, 0, 12, 0, 6],
                  [3, 9, 1, 1, 0, 12, 0, 2], [3, 10, 1, 1, 0, 12, 0, 2],
                  [3, 11, 1, 0, 0, 12, 0, 1], [3, 12, 1, 1, 0, 12, 0, 2],
                  [3, 12, 1, 2, 0, 12, 0, 1], [3, 13, 1, 1, 0, 12, 0, 1],
                  [3, 14, 1, 2, 0, 12, 0, 1], [3, 15, 1, 1, 0, 12, 0, 6],
                  [3, 16, 1, 1, 0, 12, 0, 3], [3, 17, 1, 1, 0, 12, 0, 3],
                  [3, 18, 1, 0, 0, 12, 0, 1], [3, 19, 1, 1, 0, 12, 0, 1],
                  [3, 20, 1, 1, 0, 12, 0, 1], [4, 4, 1, 1, 0, 4, 0, 4],
                  [4, 5, 1, 1, 0, 4, 0, 4], [4, 6, 1, 1, 0, 4, 0, 3],
                  [4, 7, 1, 2, 0, 4, 0, 2], [4, 8, 1, 1, 0, 4, 0, 6],
                  [4, 9, 1, 1, 0, 4, 0, 2], [4, 10, 1, 1, 0, 4, 0, 2],
                  [4, 11, 1, 0, 0, 4, 0, 1], [4, 12, 1, 1, 0, 4, 0, 2],
                  [4, 12, 1, 2, 0, 4, 0, 1], [4, 13, 1, 1, 0, 4, 0, 1],
                  [4, 14, 1, 2, 0, 4, 0, 1], [4, 15, 1, 1, 0, 4, 0, 6],
                  [4, 16, 1, 1, 0, 4, 0, 3], [4, 17, 1, 1, 0, 4, 0, 3],
                  [4, 18, 1, 0, 0, 4, 0, 1], [4, 19, 1, 1, 0, 4, 0, 1],
                  [4, 20, 1, 1, 0, 4, 0, 1], [5, 5, 1, 1, 0, 4, 0, 4],
                  [5, 6, 1, 1, 0, 4, 0, 3], [5, 7, 1, 2, 0, 4, 0, 2],
                  [5, 8, 1, 1, 0, 4, 0, 6], [5, 9, 1, 1, 0, 4, 0, 2],
                  [5, 10, 1, 1, 0, 4, 0, 2], [5, 11, 1, 0, 0, 4, 0, 1],
                  [5, 12, 1, 1, 0, 4, 0, 2], [5, 12, 1, 2, 0, 4, 0, 1],
                  [5, 13, 1, 1, 0, 4, 0, 1], [5, 14, 1, 2, 0, 4, 0, 1],
                  [5, 15, 1, 1, 0, 4, 0, 6], [5, 16, 1, 1, 0, 4, 0, 3],
                  [5, 17, 1, 1, 0, 4, 0, 3], [5, 18, 1, 0, 0, 4, 0, 1],
                  [5, 19, 1, 1, 0, 4, 0, 1], [5, 20, 1, 1, 0, 4, 0, 1],
                  [6, 6, 1, 1, 0, 3, 0, 3], [6, 7, 1, 2, 0, 3, 0, 2],
                  [6, 8, 1, 1, 0, 3, 0, 6], [6, 9, 1, 1, 0, 3, 0, 2],
                  [6, 10, 1, 1, 0, 3, 0, 2], [6, 11, 1, 0, 0, 3, 0, 1],
                  [6, 12, 1, 1, 0, 3, 0, 2], [6, 12, 1, 2, 0, 3, 0, 1],
                  [6, 13, 1, 1, 0, 3, 0, 1], [6, 14, 1, 2, 0, 3, 0, 1],
                  [6, 15, 1, 1, 0, 3, 0, 6], [6, 16, 1, 1, 0, 3, 0, 3],
                  [6, 17, 1, 1, 0, 3, 0, 3], [6, 18, 1, 0, 0, 3, 0, 1],
                  [6, 19, 1, 1, 0, 3, 0, 1], [6, 20, 1, 1, 0, 3, 0, 1],
                  [7, 7, 2, 2, 0, 2, 0, 2], [7, 8, 2, 1, 0, 2, 0, 6],
                  [7, 9, 2, 1, 0, 2, 0, 2], [7, 10, 2, 1, 0, 2, 0, 2],
                  [7, 11, 2, 0, 0, 2, 0, 1], [7, 12, 2, 1, 0, 2, 0, 2],
                  [7, 12, 2, 2, 0, 2, 0, 1], [7, 13, 2, 1, 0, 2, 0, 1],
                  [7, 14, 2, 2, 0, 2, 0, 1], [7, 15, 2, 1, 0, 2, 0, 6],
                  [7, 16, 2, 1, 0, 2, 0, 3], [7, 17, 2, 1, 0, 2, 0, 3],
                  [7, 18, 2, 0, 0, 2, 0, 1], [7, 19, 2, 1, 0, 2, 0, 1],
                  [7, 20, 2, 1, 0, 2, 0, 1], [8, 8, 1, 1, 0, 6, 0, 6],
                  [8, 9, 1, 1, 0, 6, 0, 2], [8, 10, 1, 1, 0, 6, 0, 2],
                  [8, 11, 1, 0, 0, 6, 0, 1], [8, 12, 1, 1, 0, 6, 0, 2],
                  [8, 12, 1, 2, 0, 6, 0, 1], [8, 13, 1, 1, 0, 6, 0, 1],
                  [8, 14, 1, 2, 0, 6, 0, 1], [8, 15, 1, 1, 0, 6, 0, 6],
                  [8, 16, 1, 1, 0, 6, 0, 3], [8, 17, 1, 1, 0, 6, 0, 3],
                  [8, 18, 1, 0, 0, 6, 0, 1], [8, 19, 1, 1, 0, 6, 0, 1],
                  [8, 20, 1, 1, 0, 6, 0, 1], [9, 9, 1, 1, 0, 2, 0, 2],
                  [9, 10, 1, 1, 0, 2, 0, 2], [9, 11, 1, 0, 0, 2, 0, 1],
                  [9, 12, 1, 1, 0, 2, 0, 2], [9, 12, 1, 2, 0, 2, 0, 1],
                  [9, 13, 1, 1, 0, 2, 0, 1], [9, 14, 1, 2, 0, 2, 0, 1],
                  [9, 15, 1, 1, 0, 2, 0, 6], [9, 16, 1, 1, 0, 2, 0, 3],
                  [9, 17, 1, 1, 0, 2, 0, 3], [9, 18, 1, 0, 0, 2, 0, 1],
                  [9, 19, 1, 1, 0, 2, 0, 1], [9, 20, 1, 1, 0, 2, 0, 1],
                  [10, 10, 1, 1, 0, 2, 0, 2], [10, 11, 1, 0, 0, 2, 0, 1],
                  [10, 12, 1, 1, 0, 2, 0, 2], [10, 12, 1, 2, 0, 2, 0, 1],
                  [10, 13, 1, 1, 0, 2, 0, 1], [10, 14, 1, 2, 0, 2, 0, 1],
                  [10, 15, 1, 1, 0, 2, 0, 6], [10, 16, 1, 1, 0, 2, 0, 3],
                  [10, 17, 1, 1, 0, 2, 0, 3], [10, 18, 1, 0, 0, 2, 0, 1],
                  [10, 19, 1, 1, 0, 2, 0, 1], [10, 20, 1, 1, 0, 2, 0, 1],
                  [11, 11, 0, 0, 0, 1, 0, 1], [11, 12, 0, 1, 0, 1, 0, 2],
                  [11, 12, 0, 2, 0, 1, 0, 1], [11, 13, 0, 1, 0, 1, 0, 1],
                  [11, 14, 0, 2, 0, 1, 0, 1], [11, 15, 0, 1, 0, 1, 0, 6],
                  [11, 16, 0, 1, 0, 1, 0, 3], [11, 17, 0, 1, 0, 1, 0, 3],
                  [11, 18, 0, 0, 0, 1, 0, 1], [11, 19, 0, 1, 0, 1, 0, 1],
                  [11, 20, 0, 1, 0, 1, 0, 1], [12, 12, 1, 1, 0, 2, 0, 2],
                  [12, 12, 1, 2, 0, 2, 0, 1], [12, 12, 2, 2, 0, 1, 0, 1],
                  [12, 13, 1, 1, 0, 2, 0, 1], [12, 13, 2, 1, 0, 1, 0, 1],
                  [12, 14, 1, 2, 0, 2, 0, 1], [12, 14, 2, 2, 0, 1, 0, 1],
                  [12, 15, 1, 1, 0, 2, 0, 6], [12, 15, 2, 1, 0, 1, 0, 6],
                  [12, 16, 1, 1, 0, 2, 0, 3], [12, 16, 2, 1, 0, 1, 0, 3],
                  [12, 17, 1, 1, 0, 2, 0, 3], [12, 17, 2, 1, 0, 1, 0, 3],
                  [12, 18, 1, 0, 0, 2, 0, 1], [12, 18, 2, 0, 0, 1, 0, 1],
                  [12, 19, 1, 1, 0, 2, 0, 1], [12, 19, 2, 1, 0, 1, 0, 1],
                  [12, 20, 1, 1, 0, 2, 0, 1], [12, 20, 2, 1, 0, 1, 0, 1],
                  [13, 13, 1, 1, 0, 1, 0, 1], [13, 14, 1, 2, 0, 1, 0, 1],
                  [13, 15, 1, 1, 0, 1, 0, 6], [13, 16, 1, 1, 0, 1, 0, 3],
                  [13, 17, 1, 1, 0, 1, 0, 3], [13, 18, 1, 0, 0, 1, 0, 1],
                  [13, 19, 1, 1, 0, 1, 0, 1], [13, 20, 1, 1, 0, 1, 0, 1],
                  [14, 14, 2, 2, 0, 1, 0, 1], [14, 15, 2, 1, 0, 1, 0, 6],
                  [14, 16, 2, 1, 0, 1, 0, 3], [14, 17, 2, 1, 0, 1, 0, 3],
                  [14, 18, 2, 0, 0, 1, 0, 1], [14, 19, 2, 1, 0, 1, 0, 1],
                  [14, 20, 2, 1, 0, 1, 0, 1], [15, 15, 1, 1, 0, 6, 0, 6],
                  [15, 16, 1, 1, 0, 6, 0, 3], [15, 17, 1, 1, 0, 6, 0, 3],
                  [15, 18, 1, 0, 0, 6, 0, 1], [15, 19, 1, 1, 0, 6, 0, 1],
                  [15, 20, 1, 1, 0, 6, 0, 1], [16, 16, 1, 1, 0, 3, 0, 3],
                  [16, 17, 1, 1, 0, 3, 0, 3], [16, 18, 1, 0, 0, 3, 0, 1],
                  [16, 19, 1, 1, 0, 3, 0, 1], [16, 20, 1, 1, 0, 3, 0, 1],
                  [17, 17, 1, 1, 0, 3, 0, 3], [17, 18, 1, 0, 0, 3, 0, 1],
                  [17, 19, 1, 1, 0, 3, 0, 1], [17, 20, 1, 1, 0, 3, 0, 1],
                  [18, 18, 0, 0, 0, 1, 0, 1], [18, 19, 0, 1, 0, 1, 0, 1],
                  [18, 20, 0, 1, 0, 1, 0, 1], [19, 19, 1, 1, 0, 1, 0, 1],
                  [19, 20, 1, 1, 0, 1, 0, 1], [20, 20, 1, 1, 0, 1, 0, 1]]
        assert wtasks == rtasks

        set_number_of_threads(nthreads)
