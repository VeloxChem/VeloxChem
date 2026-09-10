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


#ifndef PackedLinearAlgebra_hpp
#define PackedLinearAlgebra_hpp

#include "PackedMatrix.hpp"

namespace packlin {  // packlin namespace

/// @brief Inverts a symmetric matrix stored in the packed format.
/// @param matrix The symmetric matrix to invert.
/// @return The inverted matrix, in the packed format.
/// @note The matrix is assumed to be invertible, which holds for the matrices
/// it is used on, such as the Coulomb metric of a fitting basis. A matrix which
/// is not invertible is a critical error rather than a value the caller checks.
/// @note The matrices are expected to be positive definite and are inverted
/// through their Cholesky factorization. A matrix which is not, as a metric of
/// a nearly linearly dependent fitting basis may turn out to be, is inverted
/// through its Bunch-Kaufman factorization instead and a warning is printed,
/// rather than the inversion failing.
/// @note The inversion expands the matrix into the dense matrix, as the packed
/// factorizations are level 2 BLAS and are not threaded, and their memory is
/// not worth the several times longer run. The peak memory is therefore about
/// twice the dense matrix.
auto invert(const CPackedMatrix &matrix) -> CPackedMatrix;

/// @brief Inverts the Cholesky factor of a symmetric positive definite matrix
/// stored in the packed format.
/// @param matrix The symmetric positive definite matrix to factorize.
/// @return The inverse of the lower triangular Cholesky factor, in the packed
/// format, as a lower triangular matrix.
/// @note The matrix is factorized as L L transposed, and the returned matrix is
/// L inverted. It is what the resolution of the identity needs: the B vectors
/// formed with it satisfy B transposed times B equal to the inverse of the
/// matrix, which is what makes the Coulomb matrix of the fitting close.
/// @note The factorization costs a third of the cube of the dimensions and the
/// inversion of the factor another third, so this is cheaper than invert, which
/// pays a further third to form the whole inverse.
/// @note The matrix must be positive definite, which the metrics of the fitting
/// bases are. There is no fallback of the kind invert has, as a matrix which is
/// not positive definite has no Cholesky factor to invert.
auto cholesky_inverse(const CPackedMatrix &matrix) -> CPackedMatrix;

/// @brief Inverts the square root of a symmetric positive semidefinite matrix
/// stored in the packed format, dropping the directions of its small eigenvalues.
/// @param matrix The symmetric matrix to invert the square root of.
/// @param threshold The eigenvalues at or below which a direction is dropped.
/// @return The inverted square root, in the packed format as a symmetric matrix.
/// @note The matrix returned is V times the inverted square root of the
/// eigenvalues times V transposed, with the inverted square root set to zero for
/// the eigenvalues at or below the threshold. It is therefore the inverted square
/// root of the matrix on the directions which are kept and zero on the rest, and
/// multiplied by itself it is the pseudo inverse of the matrix.
/// @note This is what the resolution of the identity needs of a fitting basis
/// which is close to linearly dependent, where the Cholesky factorization either
/// fails or returns a factor whose inverse is meaningless. The directions which
/// carry no information are projected out rather than inverted.
/// @note The matrix is square and symmetric, and its dimensions are those of the
/// matrix it inverts whatever is dropped, so an auxiliary basis function keeps the
/// index it had. A rectangular form of one column per direction kept would be
/// smaller, and would leave the index of the auxiliary side meaning nothing.
/// @note A warning names the number of directions dropped when any are, as a
/// fitting basis which loses them is worth knowing about.
auto inverse_square_root(const CPackedMatrix &matrix, const double threshold = 1.0e-12) -> CPackedMatrix;

}  // namespace packlin

#endif /* PackedLinearAlgebra_hpp */
