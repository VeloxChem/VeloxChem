from mpi4py import MPI
from veloxchem.veloxchemlib import DenseMatrix

import numpy as np
import unittest


class TestMath(unittest.TestCase):

    def test_rectanglar_matrix(self):

        data = [[ 1., 2., 3. ], [ 4., 5., 6. ]]

        array = np.array(data)
        matrix = DenseMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = DenseMatrix(array2)

        self.assertEqual(0, np.max(np.abs(array - array2)))
        self.assertEqual(matrix, matrix2)

        self.assertEqual(2, matrix.number_of_rows())
        self.assertEqual(3, matrix.number_of_columns())

        array_t = array.T
        matrix_t = DenseMatrix(array_t)
        array2_t = matrix_t.to_numpy()
        matrix2_t = DenseMatrix(array2_t)

        self.assertEqual(0, np.max(np.abs(array_t - array2_t)))
        self.assertEqual(matrix_t, matrix2_t)

        self.assertEqual(3, matrix_t.number_of_rows())
        self.assertEqual(2, matrix_t.number_of_columns())

    def test_square_matrix(self):

        matrix = DenseMatrix(np.array([[ 1., 2. ], [ 3., 4. ]]))

        matrix.symmetrize()
        self.assertEqual(matrix, DenseMatrix([[ 2., 5. ], [ 5., 8. ]]))

        matrix.zero()
        self.assertEqual(matrix, DenseMatrix([[ 0., 0. ], [ 0., 0. ]]))

    def test_slice_matrix(self):

        matrix = DenseMatrix([[ 1., 2., 3. ], [ 4., 5., 6. ]])

        self.assertEqual(matrix.slice(0, 0, 2, 2),
                         DenseMatrix([[ 1., 2. ], [ 4., 5. ]]))

        self.assertEqual(matrix.slice(0, 1, 2, 2),
                         DenseMatrix([[ 2., 3. ], [ 5., 6. ]]))


if __name__ == "__main__":
    unittest.main()
