//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstring>
#include <memory>
#include <vector>

#ifdef ENABLE_MKL
#include "mkl.h"
#else
#include "cblas.h"
#include "lapacke.h"
#endif

#include "DenseMatrix.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "MathConst.hpp"
#include "TwoIndexes.hpp"

namespace py = pybind11;

namespace vlx_math {  // vlx_math namespace

// Helper function for printing CDenseMatrix

static std::string
CDenseMatrix_str(const CDenseMatrix& self)
{
    return self.getString();
}

// Helper function for converting CDenseMatrix to numpy array

static py::array_t<double>
CDenseMatrix_to_numpy(const CDenseMatrix& self)
{
    return vlx_general::pointer_to_numpy(self.values(), self.getNumberOfRows(), self.getNumberOfColumns());
}

// Helper function for CDenseMatrix constructor
// Not a static function; used in other files

std::shared_ptr<CDenseMatrix>
CDenseMatrix_from_numpy(const py::array_t<double>& arr)
{
    // check dimension

    std::string errdim("Matrix_from_numpy: need a 2D numpy array");

    errors::assertMsgCritical(arr.ndim() == 2, errdim);

    if (arr.data() == nullptr || arr.size() == 0)
    {
        return std::shared_ptr<CDenseMatrix>(new CDenseMatrix());
    }

    // check that the numpy array is c-style contiguous

    std::string errsrc("Matrix_from_numpy: need a contiguous numpy array");

    auto c_style = py::detail::check_flags(arr.ptr(), py::array::c_style);

    auto f_style = py::detail::check_flags(arr.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style | f_style, errsrc);

    // create CDenseMatrix from numpy array

    std::vector<double> vec(arr.size());

    if (c_style)
    {
        std::memcpy(vec.data(), arr.data(), arr.size() * sizeof(double));
    }
    else if (f_style)
    {
        for (ssize_t i = 0; i < arr.shape(0); i++)
        {
            for (ssize_t j = 0; j < arr.shape(1); j++)
            {
                vec.data()[i * arr.shape(1) + j] = arr.data()[j * arr.shape(0) + i];
            }
        }
    }

    int32_t nrows = static_cast<int32_t>(arr.shape(0));

    int32_t ncols = static_cast<int32_t>(arr.shape(1));

    return std::shared_ptr<CDenseMatrix>(new CDenseMatrix(vec, nrows, ncols));
}

static py::array_t<double>
c_matmul(const py::array_t<double>& A, const py::array_t<double>& B)
{
    // check dimension and shape

    errors::assertMsgCritical(A.ndim() == 2, "matmul - Invalid shape of matrix A");

    errors::assertMsgCritical(B.ndim() == 2, "matmul - Invalid shape of matrix B");

    auto nrow_A = A.shape(0);

    auto ncol_A = A.shape(1);

    auto nrow_B = B.shape(0);

    auto ncol_B = B.shape(1);

    errors::assertMsgCritical(ncol_A == nrow_B, "matmul - Inconsistent sizes");

    // check layout

    auto c_style_A = py::detail::check_flags(A.ptr(), py::array::c_style);

    auto f_style_A = py::detail::check_flags(A.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style_A | f_style_A, "matmul - Matrix A is noncontiguous");

    auto c_style_B = py::detail::check_flags(B.ptr(), py::array::c_style);

    auto f_style_B = py::detail::check_flags(B.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style_B | f_style_B, "matmul - Matrix B is noncontiguous");

    // check transpose

    auto trans_A = CblasNoTrans;

    auto lda_A = ncol_A;

    if (f_style_A)
    {
        trans_A = CblasTrans;

        lda_A = nrow_A;
    }

    auto trans_B = CblasNoTrans;

    auto lda_B = ncol_B;

    if (f_style_B)
    {
        trans_B = CblasTrans;

        lda_B = nrow_B;
    }

    // compute matrix-matrix multiplication

    py::array_t<double> C({nrow_A, ncol_B});

    auto lda_C = ncol_B;

    cblas_dgemm(CblasRowMajor,
                trans_A,
                trans_B,
                nrow_A,
                ncol_B,
                ncol_A,
                1.0,
                A.data(),
                lda_A,
                B.data(),
                lda_B,
                0.0,
                C.mutable_data(),
                lda_C);

    return C;
}

static py::array_t<double>
c_multi_dot(const std::vector<py::array_t<double>>& matrices)
{
    py::array_t<double> prod(matrices[0]);

    for (size_t i = 1; i < matrices.size(); i++)
    {
        auto prod_new = c_matmul(prod, matrices[i]);

        prod = py::array_t<double>(prod_new);
    }

    return prod;
}

static py::array_t<double>
c_outer(const py::array_t<double>& A, const py::array_t<double>& B)
{
    // check dimension and shape

    errors::assertMsgCritical(A.ndim() == 1, "outer - Invalid shape of vector A");

    errors::assertMsgCritical(B.ndim() == 1, "outer - Invalid shape of vector B");

    auto m = A.shape(0);

    auto n = B.shape(0);

    // check layout

    auto c_style_A = py::detail::check_flags(A.ptr(), py::array::c_style);

    auto f_style_A = py::detail::check_flags(A.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style_A & f_style_A, "outer - Vector A is noncontiguous");

    auto c_style_B = py::detail::check_flags(B.ptr(), py::array::c_style);

    auto f_style_B = py::detail::check_flags(B.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style_B & f_style_B, "outer - Vector B is noncontiguous");

    // compute outer

    py::array_t<double> C({m, n});

    std::memset(C.mutable_data(), 0, sizeof(double) * C.size());

    auto lda = n;

    cblas_dger(CblasRowMajor, m, n, 1.0, A.data(), 1, B.data(), 1, C.mutable_data(), lda);

    return C;
}

static py::list
c_eigh(const py::array_t<double>& A)
{
    // check dimension and shape

    errors::assertMsgCritical(A.ndim() == 2, "eigh - Invalid shape of matrix A");

    auto nrow_A = A.shape(0);

    auto ncol_A = A.shape(1);

    errors::assertMsgCritical(ncol_A == nrow_A, "eigh - Matrix A is not symmetric");

    // check layout

    auto c_style_A = py::detail::check_flags(A.ptr(), py::array::c_style);

    auto f_style_A = py::detail::check_flags(A.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style_A | f_style_A, "matmul - Matrix A is noncontiguous");

    auto layout_A = LAPACK_ROW_MAJOR;

    if (f_style_A) layout_A = LAPACK_COL_MAJOR;

    // copy matrix into temporary storage

    auto dim = nrow_A;

    py::array_t<double> tmp_A({dim, dim});

    std::memcpy(tmp_A.mutable_data(), A.data(), A.size() * sizeof(double));

    // initialize eigenvalues and eigenvectors

    py::array_t<double> eigenValues(dim);

    py::array_t<double> eigenVectors({dim, dim});

    // set up pointers to matrices and vectors

    auto mat = tmp_A.mutable_data();

    auto evecs = eigenVectors.mutable_data();

    auto evals = eigenValues.mutable_data();

    // temporary array for pivot data

    CMemBlock<int32_t> idx(2 * dim);

    // initialize number of eigenvalues

    int32_t nval = 0;

    // diagonalize matrix

    auto st = LAPACKE_dsyevr(
        layout_A, 'V', 'A', 'U', dim, mat, dim, 0.0, 0.0, 0, 0, 1.0e-13, &nval, evals, evecs, dim, idx.data());

    errors::assertMsgCritical(st == 0, "eigh - Diagonalization failed");

    py::list result;

    result.append(eigenValues);

    result.append(eigenVectors);

    return result;
}

// Exports classes/functions in src/math to python

void
export_math(py::module& m)
{
    // CDenseMatrix class

    py::class_<CDenseMatrix, std::shared_ptr<CDenseMatrix>>(m, "DenseMatrix")
        .def(py::init<>())
        .def(py::init<const int32_t>())
        .def(py::init<const int32_t, const int32_t>())
        .def(py::init<const CDenseMatrix&>())
        .def(py::init(&CDenseMatrix_from_numpy))
        .def("__str__", &CDenseMatrix_str)
        .def("to_numpy", &CDenseMatrix_to_numpy)
        .def("number_of_rows", &CDenseMatrix::getNumberOfRows)
        .def("number_of_columns", &CDenseMatrix::getNumberOfColumns)
        .def("symmetrize", &CDenseMatrix::symmetrize)
        .def("slice", &CDenseMatrix::slice)
        .def(py::self == py::self);

    // CTwoIndexes class

    py::class_<CTwoIndexes, std::shared_ptr<CTwoIndexes>>(m, "TwoIndexes")
        .def(py::init<>())
        .def(py::init<const int32_t, const int32_t>())
        .def("first", &CTwoIndexes::first)
        .def("second", &CTwoIndexes::second);

    // exposing functions

    m.def("mathconst_pi", &mathconst::getPiValue);

    m.def("c_matmul", &c_matmul);

    m.def("c_multi_dot", &c_multi_dot);

    m.def("c_outer", &c_outer);

    m.def("c_eigh", &c_eigh);
}

}  // namespace vlx_math
