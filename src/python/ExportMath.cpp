//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstring>
#include <memory>
#include <vector>

#ifdef ENABLE_MKL
#include <mkl.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif

#include "DenseMatrix.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "MathConst.hpp"
#include "StringFormat.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_math {  // vlx_math namespace

// Helper function for CDenseMatrix constructor
// Not a static function; used in other files

std::shared_ptr<CDenseMatrix>
CDenseMatrix_from_numpy(const py::array_t<double>& arr)
{
    // check dimension

    std::string errdim("DenseMatrix: Expecting a 2D numpy array");

    errors::assertMsgCritical(arr.ndim() == 2, errdim);

    if (arr.data() == nullptr || arr.size() == 0)
    {
        return std::make_shared<CDenseMatrix>();
    }

    // check that the numpy array is c-style contiguous

    std::string errsrc("DenseMatrix: Expecting a contiguous numpy array");

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
        for (py::ssize_t i = 0; i < arr.shape(0); i++)
        {
            for (py::ssize_t j = 0; j < arr.shape(1); j++)
            {
                vec.data()[i * arr.shape(1) + j] = arr.data()[j * arr.shape(0) + i];
            }
        }
    }

    int32_t nrows = static_cast<int32_t>(arr.shape(0));

    int32_t ncols = static_cast<int32_t>(arr.shape(1));

    return std::make_shared<CDenseMatrix>(vec, nrows, ncols);
}

static py::array_t<double>
c_matmul(const py::array_t<double>& A, const py::array_t<double>& B)
{
    // check dimension and shape

    errors::assertMsgCritical(A.ndim() == 2, "c_matmul: Invalid shape of matrix A");

    errors::assertMsgCritical(B.ndim() == 2, "c_matmul: Invalid shape of matrix B");

    auto nrow_A = A.shape(0);

    auto ncol_A = A.shape(1);

    auto nrow_B = B.shape(0);

    auto ncol_B = B.shape(1);

    errors::assertMsgCritical(ncol_A == nrow_B, "c_matmul: Inconsistent sizes");

    // check layout

    auto c_style_A = py::detail::check_flags(A.ptr(), py::array::c_style);

    auto f_style_A = py::detail::check_flags(A.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style_A | f_style_A, "c_matmul: Matrix A is noncontiguous");

    auto c_style_B = py::detail::check_flags(B.ptr(), py::array::c_style);

    auto f_style_B = py::detail::check_flags(B.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style_B | f_style_B, "c_matmul: Matrix B is noncontiguous");

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

    cblas_dgemm(CblasRowMajor, trans_A, trans_B, nrow_A, ncol_B, ncol_A, 1.0, A.data(), lda_A, B.data(), lda_B, 0.0, C.mutable_data(), lda_C);

    return C;
}

static void
c_dgemm(const std::string          layout_str,
        const std::string          trans_A_str,
        const std::string          trans_B_str,
        const int32_t              m,
        const int32_t              n,
        const int32_t              k,
        const double               alpha,
        const py::array_t<double>& A,
        const int32_t              lda_A,
        const py::array_t<double>& B,
        const int32_t              lda_B,
        const double               beta,
        py::array_t<double>&       C,
        const int32_t              lda_C)
{
    auto layout = CblasRowMajor;

    if (fstr::upcase(layout_str) == std::string("COL-MAJOR")) layout = CblasColMajor;

    auto trans_A = CblasNoTrans;

    if (fstr::upcase(trans_A_str) == std::string("T")) trans_A = CblasTrans;

    auto trans_B = CblasNoTrans;

    if (fstr::upcase(trans_B_str) == std::string("T")) trans_B = CblasTrans;

    cblas_dgemm(layout, trans_A, trans_B, m, n, k, alpha, A.data(), lda_A, B.data(), lda_B, beta, C.mutable_data(), lda_C);
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

    errors::assertMsgCritical(A.ndim() == 1, "c_outer: Invalid shape of vector A");

    errors::assertMsgCritical(B.ndim() == 1, "c_outer: Invalid shape of vector B");

    auto m = A.shape(0);

    auto n = B.shape(0);

    // check layout

    auto c_style_A = py::detail::check_flags(A.ptr(), py::array::c_style);

    auto f_style_A = py::detail::check_flags(A.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style_A & f_style_A, "c_outer: Vector A is noncontiguous");

    auto c_style_B = py::detail::check_flags(B.ptr(), py::array::c_style);

    auto f_style_B = py::detail::check_flags(B.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style_B & f_style_B, "c_outer: Vector B is noncontiguous");

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

    errors::assertMsgCritical(A.ndim() == 2, "c_eigh: Invalid shape of matrix A");

    auto nrow_A = A.shape(0);

    auto ncol_A = A.shape(1);

    errors::assertMsgCritical(ncol_A == nrow_A, "c_eigh: Matrix A is not symmetric");

    // check layout

    auto c_style_A = py::detail::check_flags(A.ptr(), py::array::c_style);

    auto f_style_A = py::detail::check_flags(A.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style_A | f_style_A, "c_eigh: Matrix A is noncontiguous");

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

    auto st = LAPACKE_dsyevr(layout_A, 'V', 'A', 'U', dim, mat, dim, 0.0, 0.0, 0, 0, 1.0e-13, &nval, evals, evecs, dim, idx.data());

    errors::assertMsgCritical(st == 0, "c_eigh: Diagonalization failed");

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

    PyClass<CDenseMatrix>(m, "DenseMatrix")
        .def(py::init<>())
        .def(py::init<const int32_t>())
        .def(py::init<const int32_t, const int32_t>())
        .def(py::init<const CDenseMatrix&>())
        .def(py::init(&CDenseMatrix_from_numpy))
        .def("__str__", &CDenseMatrix::getString)
        .def("number_of_rows", &CDenseMatrix::getNumberOfRows, "Gets number of rows in dense matrix.")
        .def("number_of_columns", &CDenseMatrix::getNumberOfColumns, "Gets number of columns in dense matrix.")
        .def("symmetrize", &CDenseMatrix::symmetrize, "Symmetrizes elements of square matrix: a_ij = a_ji = (a_ij + a_ji).")
        .def("slice",
             &CDenseMatrix::slice,
             "Creates dense matrix object by slicing specified size submatrix at selected position from this dense matrix object.",
             "iRow"_a,
             "iColumn"_a,
             "nRows"_a,
             "nColumns"_a)
        .def(
            "to_numpy",
            [](const CDenseMatrix& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.values(), self.getNumberOfRows(), self.getNumberOfColumns());
            },
            "Converts DenseMatrix to numpy array.")
        .def(py::self == py::self);

    // exposing functions

    m.def("mathconst_pi", &mathconst::getPiValue, "Gets value of PI constant.");

    m.def("c_matmul", &c_matmul, "matmul for testing purposes.");
    m.def("c_dgemm", &c_dgemm, "dgemm for testing purposes.");
    m.def("c_multi_dot", &c_multi_dot, "multi_dot for testing purposes");
    m.def("c_outer", &c_outer, "outer for testing purposes");
    m.def("c_eigh", &c_eigh, "eigh for testing purposes");
}

}  // namespace vlx_math
