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

#ifndef fika_vector_math_hpp
#define fika_vector_math_hpp

// Internal: vectorized elementary functions over contiguous arrays.

#include <cstddef>
#include <span>

namespace fika::detail {

/// out[i] = exp(-scale * x[i]) for i < out.size() (x.size() >= out.size()). Uses Accelerate's
/// vectorized exponential on macOS, std::exp elsewhere.
void exp_scaled_negative(std::span<const double> x, double scale, std::span<double> out);

/// Row-major product c = a b of an m x k matrix a (row stride lda) and a k x n matrix b (row
/// stride ldb) into the m x n matrix c (row stride ldc). A plain loop: for the few-row products
/// of the overlap kernels it is faster than Accelerate's BLAS and, unlike it, gives results that
/// do not depend on n.
void gemm(std::size_t m, std::size_t n, std::size_t k, const double* a, std::size_t lda,
          const double* b, std::size_t ldb, double* c, std::size_t ldc);

}  // namespace fika::detail

#endif  // fika_vector_math_hpp
