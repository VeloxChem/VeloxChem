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


#ifndef MathLibrary_hpp
#define MathLibrary_hpp

/// @file The declarations of the dense linear algebra of the hardware specific
/// math library.
///
/// @note This header is compiled only when VLX_USE_MATHLIB is defined, i.e.
/// when Makefile.setup has selected a math library for the platform. A file
/// including it must guard the include, and must provide the implementation
/// based on Eigen for the case where no math library is set.
///
/// @note Accelerate declares the routines itself, with its own integer type,
/// and its header is included rather than the declarations repeated, so that
/// ACCELERATE_NEW_LAPACK selects the modern LAPACK rather than the legacy
/// symbols of LAPACK 3.2.1. Elsewhere the routines are declared as the Fortran
/// symbols exported by MKL and by OpenBLAS alike.

#ifdef VLX_USE_MATHLIB

#ifdef __APPLE__

#include <Accelerate/Accelerate.h>

/// @brief The integer type of the arguments of the math library.
using lapack_int_t = __LAPACK_int;

#else

/// @brief The integer type of the arguments of the math library.
using lapack_int_t = int;

extern "C" {

/// @brief Computes the Cholesky factorization of a symmetric positive definite
/// matrix.
auto dpotrf_(const char *uplo, const lapack_int_t *n, double *a, const lapack_int_t *lda, lapack_int_t *info) -> void;

/// @brief Computes the inverse of a symmetric positive definite matrix from its
/// Cholesky factorization.
auto dpotri_(const char *uplo, const lapack_int_t *n, double *a, const lapack_int_t *lda, lapack_int_t *info) -> void;

/// @brief Computes the Bunch-Kaufman factorization of a symmetric matrix.
auto dsytrf_(const char          *uplo,
             const lapack_int_t  *n,
             double              *a,
             const lapack_int_t  *lda,
             lapack_int_t        *ipiv,
             double              *work,
             const lapack_int_t  *lwork,
             lapack_int_t        *info) -> void;

/// @brief Computes the inverse of a triangular matrix.
auto dtrtri_(const char         *uplo,
             const char         *diag,
             const lapack_int_t *n,
             double             *a,
             const lapack_int_t *lda,
             lapack_int_t       *info) -> void;

/// @brief Computes a general matrix product.
auto dgemm_(const char         *transa,
            const char         *transb,
            const lapack_int_t *m,
            const lapack_int_t *n,
            const lapack_int_t *k,
            const double       *alpha,
            const double       *a,
            const lapack_int_t *lda,
            const double       *b,
            const lapack_int_t *ldb,
            const double       *beta,
            double             *c,
            const lapack_int_t *ldc) -> void;

/// @brief Computes the inverse of a symmetric matrix from its Bunch-Kaufman
/// factorization.
auto dsytri_(const char         *uplo,
             const lapack_int_t *n,
             double             *a,
             const lapack_int_t *lda,
             const lapack_int_t *ipiv,
             double             *work,
             lapack_int_t       *info) -> void;

}  // extern "C"

#endif /* __APPLE__ */

#endif /* VLX_USE_MATHLIB */

#endif /* MathLibrary_hpp */
