from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import Molecule
from veloxchem.veloxchemlib import make_gto_blocks
from veloxchem.veloxchemlib import set_number_of_threads
from veloxchem.veloxchemlib import get_number_of_threads
from veloxchem.veloxchemlib import make_workgroup


class TestOpenMPFunc:
    """
    Implements tests for src/general/OpenMPFunc.hpp
    """

    def get_data(self):

        h2ostr = """O   0.000   0.000  -1.000
                    H   0.000   1.400  -2.100
                    H   0.000  -1.400  -2.100"""

        mol = Molecule.read_str(h2ostr, 'au')

        bas = MolecularBasis.read(mol, 'DEF2-SVP', 'basis')

        return (mol, bas)

    def test_number_of_threads(self):

        nthreads = get_number_of_threads()

        set_number_of_threads(128)
        assert get_number_of_threads() == 128

        set_number_of_threads(nthreads)

    def test_make_workgroup(self):

        nthreads = get_number_of_threads()

        mol_h2o, bas_svp = self.get_data()

        # check workgroup with single vector of basis function blocks
        set_number_of_threads(3)
        gblocks = make_gto_blocks(bas_svp, mol_h2o)
        wgroups = make_workgroup(gblocks)
        rgroups = [  # task 0
            [  # gtos: (0, 1)
                [0, 0, 0, 2],
                [0, 1, 0, 2],
                [0, 2, 0, 2],
                [0, 3, 0, 2],
                [0, 4, 0, 2],
                [0, 5, 0, 2],
                # gtos: (0, 3)
                [1, 1, 0, 1],
                [1, 2, 0, 1],
                [3, 1, 0, 1],
                [1, 4, 0, 1],
                [1, 5, 0, 1],
                # gtos: (0, 5)
                [2, 2, 0, 1],
                [3, 2, 0, 1],
                [4, 2, 0, 1],
                [5, 2, 0, 1],
                # gtos: (1, 1)
                [3, 3, 0, 1],
                [3, 4, 0, 1],
                [3, 5, 0, 1],
                # gtos: (1, 3)
                [4, 4, 0, 1],
                [5, 4, 0, 1],
                # gtos: (2, 1)
                [5, 5, 0, 1]
            ],
            # task 1
            [  # gtos: (0, 1)
                [0, 0, 2, 3],
                [0, 1, 2, 3],
                [0, 2, 2, 3],
                [0, 3, 2, 3],
                [0, 4, 2, 3],
                [0, 5, 2, 3],
                # gtos: (0, 3)
                [1, 1, 1, 2],
                [1, 2, 1, 2],
                [3, 1, 1, 2],
                [1, 4, 1, 2],
                [1, 5, 1, 2],
                # gtos: (0, 5)
                [3, 2, 1, 2],
                # gtos: (1, 1)
                [3, 3, 1, 2],
                [3, 4, 1, 2],
                [3, 5, 1, 2]
            ],
            # task 2
            [  # gtos: (0, 1)
                [0, 0, 3, 4],
                [0, 1, 3, 4],
                [0, 2, 3, 4],
                [0, 3, 3, 4],
                [0, 4, 3, 4],
                [0, 5, 3, 4],
                # gtos: (0, 3)
                [3, 1, 2, 3],
                # gtos: (0, 3)
                [3, 2, 2, 3],
                # gtos: (1, 1)
                [3, 3, 2, 3],
                [3, 4, 2, 3],
                [3, 5, 2, 3]
            ]
        ]
        assert wgroups == rgroups

        # check workgroup with two vectors of basis function blocks
        wgroups = make_workgroup(gblocks, gblocks)
        rgroups = [  # task 0
            [  # gtos: (0, 1)
                [0, 0, 0, 2],
                [0, 1, 0, 2],
                [0, 2, 0, 2],
                [0, 3, 0, 2],
                [0, 4, 0, 2],
                [0, 5, 0, 2],
                # gtos: (0, 3)
                [1, 0, 0, 1],
                [1, 1, 0, 1],
                [1, 2, 0, 1],
                [1, 3, 0, 1],
                [1, 4, 0, 1],
                [1, 5, 0, 1],
                # gtos: (0, 5)
                [2, 0, 0, 1],
                [2, 1, 0, 1],
                [2, 2, 0, 1],
                [2, 3, 0, 1],
                [2, 4, 0, 1],
                [2, 5, 0, 1],
                # gtos: (1, 1)
                [3, 0, 0, 1],
                [3, 1, 0, 1],
                [3, 2, 0, 1],
                [3, 3, 0, 1],
                [3, 4, 0, 1],
                [3, 5, 0, 1],
                # gtos: (1, 3)
                [4, 0, 0, 1],
                [4, 1, 0, 1],
                [4, 2, 0, 1],
                [4, 3, 0, 1],
                [4, 4, 0, 1],
                [4, 5, 0, 1],
                # gtos: (2, 1)
                [5, 0, 0, 1],
                [5, 1, 0, 1],
                [5, 2, 0, 1],
                [5, 3, 0, 1],
                [5, 4, 0, 1],
                [5, 5, 0, 1]
            ],
            # task 1
            [  # gtos: (0, 1)
                [0, 0, 2, 3],
                [0, 1, 2, 3],
                [0, 2, 2, 3],
                [0, 3, 2, 3],
                [0, 4, 2, 3],
                [0, 5, 2, 3],
                # gtos: (0, 3)
                [1, 0, 1, 2],
                [1, 1, 1, 2],
                [1, 2, 1, 2],
                [1, 3, 1, 2],
                [1, 4, 1, 2],
                [1, 5, 1, 2],
                # gtos: (1, 1)
                [3, 0, 1, 2],
                [3, 1, 1, 2],
                [3, 2, 1, 2],
                [3, 3, 1, 2],
                [3, 4, 1, 2],
                [3, 5, 1, 2]
            ],
            # task 2
            [  # gtos: (0, 1)
                [0, 0, 3, 4],
                [0, 1, 3, 4],
                [0, 2, 3, 4],
                [0, 3, 3, 4],
                [0, 4, 3, 4],
                [0, 5, 3, 4],
                # gtos: (1, 1)
                [3, 0, 2, 3],
                [3, 1, 2, 3],
                [3, 2, 2, 3],
                [3, 3, 2, 3],
                [3, 4, 2, 3],
                [3, 5, 2, 3]
            ]
        ]
        assert wgroups == rgroups

        set_number_of_threads(nthreads)
