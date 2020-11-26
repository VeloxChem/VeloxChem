import numpy as np
import unittest

from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import c_matmul
from veloxchem.veloxchemlib import c_dgemm
from veloxchem.veloxchemlib import c_multi_dot
from veloxchem.veloxchemlib import c_outer
from veloxchem.veloxchemlib import c_eigh


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

        mat_C = c_matmul(mat_A, mat_B)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_A = np.arange(6.).reshape(3, 2)
        mat_B = np.arange(12.).reshape(3, 4)
        ref_C = np.matmul(mat_A.T, mat_B)

        mat_C = c_matmul(mat_A.T, mat_B)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_C = c_matmul(mat_A.T.copy(), mat_B)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_A = np.arange(6.).reshape(2, 3)
        mat_B = np.arange(12.).reshape(4, 3)
        ref_C = np.matmul(mat_A, mat_B.T)

        mat_C = c_matmul(mat_A, mat_B.T)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_C = c_matmul(mat_A, mat_B.T.copy())
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_A = np.arange(6.).reshape(3, 2)
        mat_B = np.arange(12.).reshape(4, 3)
        ref_C = np.matmul(mat_A.T, mat_B.T)

        mat_C = c_matmul(mat_A.T, mat_B.T)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_C = c_matmul(mat_A.T.copy(), mat_B.T)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_C = c_matmul(mat_A.T, mat_B.T.copy())
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_C = c_matmul(mat_A.T.copy(), mat_B.T.copy())
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

    def test_dgemm_square_matrix(self):

        mat_A = np.arange(90.).reshape(9, 10)[:5, :5]
        mat_B = np.arange(56.).reshape(7, 8)[:5, :5]
        mat_C = np.zeros((5, 5))

        ref_C = np.matmul(mat_A, mat_B)
        c_dgemm('row-major', 'n', 'n', 5, 5, 5, 1.0, mat_A, 10, mat_B, 8, 0.0,
                mat_C, 5)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        ref_C += 2.0 * np.matmul(mat_A.T, mat_B)
        c_dgemm('row-major', 't', 'n', 5, 5, 5, 2.0, mat_A, 10, mat_B, 8, 1.0,
                mat_C, 5)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        ref_C *= 2.0
        ref_C += 3.0 * np.matmul(mat_A, mat_B.T)
        c_dgemm('row-major', 'n', 't', 5, 5, 5, 3.0, mat_A, 10, mat_B, 8, 2.0,
                mat_C, 5)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        ref_C += np.matmul(mat_A.T, mat_B.T)
        c_dgemm('row-major', 't', 't', 5, 5, 5, 1.0, mat_A, 10, mat_B, 8, 1.0,
                mat_C, 5)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

    def test_dgemm_rectanglar_matrix(self):

        mat_A = np.arange(90.).reshape(9, 10)[:3, :4]
        mat_B = np.arange(56.).reshape(7, 8)[:4, :5]
        mat_C = np.zeros((3, 5))

        ref_C = np.matmul(mat_A, mat_B)
        c_dgemm('row-major', 'n', 'n', 3, 5, 4, 1.0, mat_A, 10, mat_B, 8, 0.0,
                mat_C, 5)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_A = np.arange(90.).reshape(9, 10)[:3, :4]
        mat_B = np.arange(56.).reshape(7, 8)[:5, :4]

        ref_C += 2.0 * np.matmul(mat_A, mat_B.T)
        c_dgemm('row-major', 'n', 't', 3, 5, 4, 2.0, mat_A, 10, mat_B, 8, 1.0,
                mat_C, 5)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_A = np.arange(90.).reshape(9, 10)[:4, :3]
        mat_B = np.arange(56.).reshape(7, 8)[:4, :5]

        ref_C *= 2.0
        ref_C += 3.0 * np.matmul(mat_A.T, mat_B)
        c_dgemm('row-major', 't', 'n', 3, 5, 4, 3.0, mat_A, 10, mat_B, 8, 2.0,
                mat_C, 5)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        mat_A = np.arange(90.).reshape(9, 10)[:4, :3]
        mat_B = np.arange(56.).reshape(7, 8)[:5, :4]

        ref_C += np.matmul(mat_A.T, mat_B.T)
        c_dgemm('row-major', 't', 't', 3, 5, 4, 1.0, mat_A, 10, mat_B, 8, 1.0,
                mat_C, 5)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

    def test_multi_dot(self):

        mat_A = np.arange(6.).reshape(2, 3)
        mat_B = np.arange(12.).reshape(3, 4)
        mat_C = np.arange(16.).reshape(4, 4)

        ref_prod = np.linalg.multi_dot([mat_A, mat_B, mat_C.T])
        prod = c_multi_dot([mat_A, mat_B, mat_C.T])
        self.assertTrue(np.max(np.abs(prod - ref_prod)) < 1.0e-13)

        ref_prod = np.linalg.multi_dot([mat_A, mat_B, mat_C, mat_B.T, mat_A.T])
        prod = c_multi_dot([mat_A, mat_B, mat_C, mat_B.T, mat_A.T])
        self.assertTrue(np.max(np.abs(prod - ref_prod)) < 1.0e-13)

    def test_outer(self):

        vec_A = np.arange(1., 10.)
        vec_B = np.arange(10., 100.)

        ref_C = np.outer(vec_A, vec_B)
        mat_C = c_outer(vec_A, vec_B)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        ref_C = np.outer(vec_B, vec_A)
        mat_C = c_outer(vec_B, vec_A)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        ref_C = np.outer(vec_A, vec_A)
        mat_C = c_outer(vec_A, vec_A)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

        ref_C = np.outer(vec_B, vec_B)
        mat_C = c_outer(vec_B, vec_B)
        self.assertTrue(np.max(np.abs(mat_C - ref_C)) < 1.0e-13)

    def test_eigh(self):

        mat_A = np.arange(10., 20.)
        np.random.shuffle(mat_A)
        mat_A = np.diag(mat_A)
        mat_A += np.random.uniform(0.01, 0.99, 100).reshape(10, 10)
        mat_A = 0.5 * (mat_A + mat_A.T)
        ref_eigvals, ref_eigvecs = np.linalg.eigh(mat_A)

        eigvals, eigvecs = c_eigh(mat_A)
        self.assertTrue(np.max(np.abs(eigvals - ref_eigvals)) < 1.0e-13)
        for k in range(ref_eigvecs.shape[1]):
            vec = eigvecs[:, k].copy()
            ref_vec = ref_eigvecs[:, k].copy()
            if np.dot(vec, ref_vec) < 0.0:
                vec *= -1.0
            self.assertTrue(np.max(np.abs(vec - ref_vec)) < 1.0e-13)


if __name__ == "__main__":
    unittest.main()
