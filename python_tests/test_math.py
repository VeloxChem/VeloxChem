import numpy as np
import unittest

from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import matmul


class TestMath(unittest.TestCase):

    def test_matrix_numpy(self):

        data = [[1., 2., 3.], [4., 5., 6.]]

        array = np.array(data)
        matrix = DenseMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = DenseMatrix(array2)

        self.assertTrue((array == array2).all())
        self.assertEqual(matrix, matrix2)

        self.assertEqual(2, matrix.number_of_rows())
        self.assertEqual(3, matrix.number_of_columns())

        array_t = array.T
        matrix_t = DenseMatrix(array_t)
        array2_t = matrix_t.to_numpy()
        matrix2_t = DenseMatrix(array2_t)

        self.assertTrue((array_t == array2_t).all())
        self.assertEqual(matrix_t, matrix2_t)

        self.assertEqual(3, matrix_t.number_of_rows())
        self.assertEqual(2, matrix_t.number_of_columns())

    def test_symmetrize(self):

        matrix = DenseMatrix(np.array([[1., 2.], [3., 4.]]))

        matrix.symmetrize()
        self.assertEqual(matrix, DenseMatrix([[2., 5.], [5., 8.]]))

    def test_slice_matrix(self):

        matrix = DenseMatrix([[1., 2., 3.], [4., 5., 6.]])

        self.assertEqual(matrix.slice(0, 0, 2, 2),
                         DenseMatrix([[1., 2.], [4., 5.]]))

        self.assertEqual(matrix.slice(0, 1, 2, 2),
                         DenseMatrix([[2., 3.], [5., 6.]]))

    def test_matmul(self):

        mat_A = np.arange(6.).reshape(2, 3)
        mat_B = np.arange(12.).reshape(3, 4)
        ref_C = np.matmul(mat_A, mat_B)

        mat_C = matmul(mat_A, mat_B)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_A = np.arange(6.).reshape(3, 2)
        mat_B = np.arange(12.).reshape(3, 4)
        ref_C = np.matmul(mat_A.T, mat_B)

        mat_C = matmul(mat_A.T, mat_B)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_C = matmul(mat_A.T.copy(), mat_B)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_A = np.arange(6.).reshape(2, 3)
        mat_B = np.arange(12.).reshape(4, 3)
        ref_C = np.matmul(mat_A, mat_B.T)

        mat_C = matmul(mat_A, mat_B.T)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_C = matmul(mat_A, mat_B.T.copy())
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_A = np.arange(6.).reshape(3, 2)
        mat_B = np.arange(12.).reshape(4, 3)
        ref_C = np.matmul(mat_A.T, mat_B.T)

        mat_C = matmul(mat_A.T, mat_B.T)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_C = matmul(mat_A.T.copy(), mat_B.T)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_C = matmul(mat_A.T, mat_B.T.copy())
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_C = matmul(mat_A.T.copy(), mat_B.T.copy())
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)


if __name__ == "__main__":
    unittest.main()
