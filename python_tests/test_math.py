from mpi4py import MPI
from veloxchem.veloxchemlib import DenseMatrix

import numpy as np
import unittest


class TestMath(unittest.TestCase):

    def test_numpy(self):

        data = [[ 1., 2., ], [ 3., 4., ]]

        array = np.array(data)
        matrix = DenseMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = DenseMatrix(array2)

        self.assertEqual(0, np.max(np.abs(array - array2)))
        self.assertEqual(matrix, matrix2)


if __name__ == "__main__":
    unittest.main()
