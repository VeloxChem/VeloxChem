//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "PackedLinearAlgebra.hpp"

#include <memory>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"

#ifdef VLX_USE_MATHLIB
#include "MathLibrary.hpp"
#else
#include "Eigen/Dense"
#endif

namespace packlin {  // packlin namespace

/// @brief The message printed when the matrix turns out not to be positive
/// definite and the Cholesky factorization is abandoned.
static const std::string _indefinite_message =
    "PackedMatrix inversion: The matrix is not positive definite, and is inverted through its Bunch-Kaufman "
    "factorization. This is expected of a nearly linearly dependent basis, and is a sign of an ill conditioned "
    "matrix otherwise.";

#ifdef VLX_USE_MATHLIB

/// @brief Inverts the symmetric matrix held in the dense array, in place.
/// @param values The values of the dense matrix, as a row major array of ndim
/// rows and ndim columns, which is overwritten by the inverted matrix.
/// @param ndim The number of rows of matrix.
/// @return True if the matrix was positive definite and has been inverted,
/// false if it was not and has been left partly factorized.
/// @note The array is row major and the math library is column major, so the
/// upper triangle of the library is the lower triangle of the array. The matrix
/// is symmetric, so the two hold the same elements and the inverted matrix is
/// written into the lower triangle of the array, which is the triangle the
/// packed matrix stores.
static auto
_invert_dense(double *values, const size_t ndim) -> bool
{
    const char uplo = 'U';

    auto ndim_arg = static_cast<lapack_int_t>(ndim);

    lapack_int_t info = 0;

    dpotrf_(&uplo, &ndim_arg, values, &ndim_arg, &info);

    errors::assertMsgCritical(info >= 0, "PackedMatrix inversion: Invalid argument of the Cholesky factorization");

    if (info == 0)
    {
        dpotri_(&uplo, &ndim_arg, values, &ndim_arg, &info);

        errors::assertMsgCritical(info == 0, "PackedMatrix inversion: The matrix is singular");

        return true;
    }

    // NOTE: the Cholesky factorization has overwritten the leading part of the
    // matrix, which the Bunch-Kaufman factorization cannot pick up, so the
    // caller expands the matrix again before it factorizes it the other way.

    errors::msg(_indefinite_message, "Warning");

    return false;
}

/// @brief Inverts the symmetric indefinite matrix held in the dense array, in
/// place.
/// @param values The values of the dense matrix, as a row major array of ndim
/// rows and ndim columns, which is overwritten by the inverted matrix.
/// @param ndim The number of rows of matrix.
static auto
_invert_dense_indefinite(double *values, const size_t ndim) -> void
{
    const char uplo = 'U';

    auto ndim_arg = static_cast<lapack_int_t>(ndim);

    lapack_int_t info = 0;

    std::vector<lapack_int_t> ipiv(ndim);

    // NOTE: the size of the work array is the one the factorization asks for,
    // which it reports when it is called with a size of minus one.

    lapack_int_t lwork = -1;

    double work_size = 0.0;

    dsytrf_(&uplo, &ndim_arg, values, &ndim_arg, ipiv.data(), &work_size, &lwork, &info);

    errors::assertMsgCritical(info == 0, "PackedMatrix inversion: Failed to size the Bunch-Kaufman factorization");

    lwork = static_cast<lapack_int_t>(work_size);

    std::vector<double> work(static_cast<size_t>(lwork));

    dsytrf_(&uplo, &ndim_arg, values, &ndim_arg, ipiv.data(), work.data(), &lwork, &info);

    errors::assertMsgCritical(info == 0, "PackedMatrix inversion: The matrix is singular");

    // NOTE: the inversion asks for a work array of one element per row, which is
    // smaller than the one of the factorization and is reused rather than
    // allocated again.

    work.resize(ndim);

    dsytri_(&uplo, &ndim_arg, values, &ndim_arg, ipiv.data(), work.data(), &info);

    errors::assertMsgCritical(info == 0, "PackedMatrix inversion: The matrix is singular");
}

#else

/// @brief Inverts the symmetric matrix held in the dense array, in place.
/// @param values The values of the dense matrix, as a row major array of ndim
/// rows and ndim columns, which is overwritten by the inverted matrix.
/// @param ndim The number of rows of matrix.
/// @return True if the matrix was positive definite and has been inverted,
/// false if it was not and has been left partly factorized.
/// @note The matrix is mapped as column major, which is the transposed matrix
/// and is the matrix itself, as it is symmetric. The factorization is made in
/// place, through a reference to the mapped array, as a decomposition of its
/// own would hold a second copy of the matrix.
static auto
_invert_dense(double *values, const size_t ndim) -> bool
{
    Eigen::Map<Eigen::MatrixXd> matrix(values, static_cast<Eigen::Index>(ndim), static_cast<Eigen::Index>(ndim));

    Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt(matrix);

    if (llt.info() != Eigen::Success)
    {
        errors::msg(_indefinite_message, "Warning");

        return false;
    }

    auto inverse = Eigen::MatrixXd::Identity(static_cast<Eigen::Index>(ndim), static_cast<Eigen::Index>(ndim)).eval();

    llt.solveInPlace(inverse);

    std::copy(inverse.data(), inverse.data() + ndim * ndim, values);

    return true;
}

/// @brief Inverts the symmetric indefinite matrix held in the dense array, in
/// place.
/// @param values The values of the dense matrix, as a row major array of ndim
/// rows and ndim columns, which is overwritten by the inverted matrix.
/// @param ndim The number of rows of matrix.
static auto
_invert_dense_indefinite(double *values, const size_t ndim) -> void
{
    Eigen::Map<Eigen::MatrixXd> matrix(values, static_cast<Eigen::Index>(ndim), static_cast<Eigen::Index>(ndim));

    Eigen::LDLT<Eigen::Ref<Eigen::MatrixXd>> ldlt(matrix);

    errors::assertMsgCritical(ldlt.info() == Eigen::Success, "PackedMatrix inversion: The matrix is singular");

    auto inverse = Eigen::MatrixXd::Identity(static_cast<Eigen::Index>(ndim), static_cast<Eigen::Index>(ndim)).eval();

    ldlt.solveInPlace(inverse);

    std::copy(inverse.data(), inverse.data() + ndim * ndim, values);
}

#endif /* VLX_USE_MATHLIB */

auto
invert(const CPackedMatrix &matrix) -> CPackedMatrix
{
    errors::assertMsgCritical(matrix.get_type() == mat_t::symmetric, "PackedMatrix inversion: The matrix must be symmetric");

    const auto ndim = matrix.number_of_rows();

    errors::assertMsgCritical(ndim > 0, "PackedMatrix inversion: The matrix must not be empty");

    // NOTE: the dense matrix is not zeroed, as the expansion writes every one of
    // its elements, and it reaches gigabytes for the fitting bases the inversion
    // is used on.

    auto dense = std::make_unique_for_overwrite<double[]>(ndim * ndim);

    matrix.to_dense(dense.get());

    // NOTE: the Cholesky factorization leaves the matrix as it found it only
    // when it succeeds, so a matrix which turned out not to be positive definite
    // is expanded again before it is factorized the other way.

    if (!_invert_dense(dense.get(), ndim))
    {
        matrix.to_dense(dense.get());

        _invert_dense_indefinite(dense.get(), ndim);
    }

    auto inverse = CPackedMatrix(ndim, ndim, mat_t::symmetric);

    inverse.from_dense(dense.get());

    return inverse;
}

}  // namespace packlin

namespace packlin {  // packlin namespace

/// @brief Inverts the Cholesky factor of the symmetric matrix held in the dense
/// array, in place.
/// @param values The values of the dense matrix, as a row major array of ndim
/// rows and ndim columns, whose lower triangle is overwritten by the inverted
/// factor.
/// @param ndim The number of rows of matrix.
static auto
_invert_dense_factor(double *values, const size_t ndim) -> void
{
#ifdef VLX_USE_MATHLIB

    // NOTE: the array is row major and the library is column major, so the upper
    // triangle of the library is the lower triangle of the array. The
    // factorization of the upper triangle therefore leaves the lower triangle of
    // the array holding L, with the matrix equal to L L transposed, and the
    // inversion of that triangle leaves it holding L inverted.

    const char uplo = 'U';

    const char diag = 'N';

    auto ndim_arg = static_cast<lapack_int_t>(ndim);

    lapack_int_t info = 0;

    dpotrf_(&uplo, &ndim_arg, values, &ndim_arg, &info);

    errors::assertMsgCritical(info >= 0, "PackedMatrix Cholesky inversion: Invalid argument of the factorization");

    errors::assertMsgCritical(info == 0, "PackedMatrix Cholesky inversion: The matrix is not positive definite");

    dtrtri_(&uplo, &diag, &ndim_arg, values, &ndim_arg, &info);

    errors::assertMsgCritical(info == 0, "PackedMatrix Cholesky inversion: The factor is singular");

#else

    const auto nrows = static_cast<Eigen::Index>(ndim);

    // NOTE: the array is mapped as column major, which is the transposed matrix
    // and is the matrix itself, as it is symmetric. Eigen factorizes it as L L
    // transposed with L in the column major lower triangle, which is the upper
    // triangle of the row major array, so the factor is transposed into the lower
    // triangle the packed matrix stores.

    Eigen::Map<Eigen::MatrixXd> matrix(values, nrows, nrows);

    Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt(matrix);

    errors::assertMsgCritical(llt.info() == Eigen::Success,
                              "PackedMatrix Cholesky inversion: The matrix is not positive definite");

    auto factor = Eigen::MatrixXd::Identity(nrows, nrows).eval();

    llt.matrixL().solveInPlace(factor);

    // NOTE: the inverted factor is column major and lower triangular, and the
    // array is read as row major, so it is stored transposed.

    for (Eigen::Index i = 0; i < nrows; i++)
    {
        for (Eigen::Index j = 0; j <= i; j++)
        {
            values[static_cast<size_t>(i) * ndim + static_cast<size_t>(j)] = factor(i, j);
        }
    }

#endif /* VLX_USE_MATHLIB */
}

auto
cholesky_inverse(const CPackedMatrix &matrix) -> CPackedMatrix
{
    errors::assertMsgCritical(matrix.get_type() == mat_t::symmetric,
                              std::string("PackedMatrix Cholesky inversion: The matrix must be symmetric"));

    const auto ndim = matrix.number_of_rows();

    errors::assertMsgCritical(ndim > 0, std::string("PackedMatrix Cholesky inversion: The matrix must not be empty"));

    auto dense = std::make_unique_for_overwrite<double[]>(ndim * ndim);

    matrix.to_dense(dense.get());

    _invert_dense_factor(dense.get(), ndim);

    auto factor = CPackedMatrix(ndim, ndim, mat_t::lower_triangular);

    factor.from_dense(dense.get());

    return factor;
}

}  // namespace packlin
