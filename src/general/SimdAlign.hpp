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


#ifndef SimdAlign_hpp
#define SimdAlign_hpp

#include <cstddef>

namespace simd {  // simd namespace

/// @brief Gets size of cache line of the target architecture.
/// @return The size of cache line in bytes.
/// @note Apple silicon and other AArch64 targets use 128 byte cache lines, while
/// x86-64 targets use 64 byte cache lines.
constexpr auto
cache_line_size() -> size_t
{
#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
    return 128;
#else
    return 64;
#endif
}

/// @brief Checks if the exponential function is vectorized on the target.
/// @return True if a vector exponential is available, false otherwise.
/// @note The exponential is issued as a scalar call inside a vector loop, as no
/// vector math library is linked. Vectorizing a loop which contains it is then a
/// loss, since the values are unpacked and repacked around each call. Kernels use
/// this predicate to condition their vectorization, and turn it on once a vector
/// exponential is in place.
constexpr auto
has_vector_exp() -> bool
{
    return false;
}

/// @brief Gets number of values which fit into a cache line.
/// @return The number of values in a cache line.
constexpr auto
values_per_cache_line() -> size_t
{
    return cache_line_size() / sizeof(double);
}

/// @brief Computes the padded number of values in a row, so that the row which
/// follows it starts at a cache line boundary.
/// @param ncolumns The number of values in a row.
/// @return The padded number of values in a row.
constexpr auto
pitch_of(const size_t ncolumns) -> size_t
{
    const auto nvals = values_per_cache_line();

    return ((ncolumns + nvals - 1) / nvals) * nvals;
}

}  // namespace simd

#endif /* SimdAlign_hpp */
