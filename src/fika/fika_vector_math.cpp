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

#include "fika_vector_math.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#ifdef VLX_USE_MATHLIB
#include "MathLibrary.hpp"
#endif

namespace fika::detail {

void exp_scaled_negative(std::span<const double> x, double scale, std::span<double> out) {
  assert(x.size() >= out.size());
  for (std::size_t i = 0; i < out.size(); ++i) {
    out[i] = -scale * x[i];
  }
#if defined(VLX_USE_MATHLIB) && defined(__APPLE__)
  const int n = static_cast<int>(out.size());
  vvexp(out.data(), out.data(), &n);
#else
  for (double& value : out) {
    value = std::exp(value);
  }
#endif
}

void gemm(std::size_t m, std::size_t n, std::size_t k, const double* a, std::size_t lda,
          const double* b, std::size_t ldb, double* c, std::size_t ldc) {
  assert(lda >= k && ldb >= n && ldc >= n);
  for (std::size_t i = 0; i < m; ++i) {
    double* row = c + i * ldc;
    std::fill(row, row + n, 0.0);
    for (std::size_t p = 0; p < k; ++p) {
      const double factor = a[i * lda + p];
      if (factor == 0.0) {
        continue;  // zero contraction coefficients are common in general contractions
      }
      const double* source = b + p * ldb;
      for (std::size_t j = 0; j < n; ++j) {
        row[j] += factor * source[j];
      }
    }
  }
}

void blas_gemm_transposed(std::size_t m, std::size_t n, std::size_t k, const double* a,
                          std::size_t lda, const double* b, std::size_t ldb, double* c,
                          std::size_t ldc) {
  assert(lda >= k && ldb >= k && ldc >= n);
  if (m == 0 || n == 0) {
    return;
  }
#ifdef VLX_USE_MATHLIB
  // Row-major C = A B^T is column-major C^T = B A^T: B (k x n column-major) transposed.
  const auto mm = static_cast<lapack_int_t>(n), nn = static_cast<lapack_int_t>(m);
  const auto kk = static_cast<lapack_int_t>(k);
  const auto lda_ = static_cast<lapack_int_t>(ldb), ldb_ = static_cast<lapack_int_t>(lda);
  const auto ldc_ = static_cast<lapack_int_t>(ldc);
  const double one = 1.0, zero = 0.0;
  dgemm_("T", "N", &mm, &nn, &kk, &one, b, &lda_, a, &ldb_, &zero, c, &ldc_);
#else
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (std::size_t p = 0; p < k; ++p) {
        sum += a[i * lda + p] * b[j * ldb + p];
      }
      c[i * ldc + j] = sum;
    }
  }
#endif
}

}  // namespace fika::detail
