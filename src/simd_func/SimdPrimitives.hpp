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


#ifndef SimdPrimitives_hpp
#define SimdPrimitives_hpp

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <string>
#include <vector>

#include "BasisFunction.hpp"
#include "ErrorHandler.hpp"
#include "SimdMatrix.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief The data of one pair of primitives of the basis functions on bra and
/// ket sides which reaches at least one atom pair.
/// @note The exponents and the normalization factors are passed by value, as the
/// accumulation of a pair of primitives reads each of them a handful of times and
/// forms its prefactors from them before the loop over the atom pairs.
struct CPrimitivePair
{
    /// @brief The exponent of the primitive on bra side.
    double aexp;

    /// @brief The exponent of the primitive on ket side.
    double bexp;

    /// @brief The normalization factor of the primitive on bra side.
    double anorm;

    /// @brief The normalization factor of the primitive on ket side.
    double bnorm;

    /// @brief The number of atom pairs the pair of primitives reaches, which is
    /// the number of leading columns of the accumulation buffer it contributes to.
    size_t ncols;
};

/// @brief Creates the buffer the pairs of primitives accumulate their
/// contributions in.
/// @param dimensions The number of atom pairs each pair of primitives reaches.
/// @param nrows The number of accumulators of the kernel.
/// @return The zeroed matrix of nrows rows spanning the atom pairs reached by the
/// pair of primitives reaching furthest, empty when no pair of primitives reaches
/// any atom pair.
/// @note The buffer spans the atom pairs reached by the pair of primitives
/// reaching furthest, which is searched for rather than assumed. The primitives
/// are sorted by descending exponent, but the bound of a pair of primitives
/// carries their prefactor as well as their decay, so a tighter pair with a larger
/// prefactor reaches further than a more diffuse pair with a smaller one, and the
/// last pair is not always the furthest reaching.
inline auto
make_primitive_buffer(const std::vector<size_t> &dimensions, const size_t nrows) -> CSimdMatrix
{
    errors::assertMsgCritical(nrows > 0, std::string("SimdPrimitives.make_primitive_buffer: Number of rows must be positive"));

    if (dimensions.empty()) return CSimdMatrix();

    const auto nmax = *std::ranges::max_element(dimensions);

    if (nmax == 0) return CSimdMatrix();

    auto matrix = CSimdMatrix(nrows, nmax);

    matrix.zero();

    return matrix;
}

/// @brief Accumulates the contribution of every pair of primitives which reaches
/// at least one atom pair.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param dimensions The number of atom pairs each pair of primitives reaches,
/// with the primitives on bra side as the slowest running index.
/// @param accumulate The accumulation of one pair of primitives, called as
/// accumulate(pair) with the data of the pair of primitives.
/// @note The pairs of primitives which reach no atom pair are skipped rather than
/// passed on with a length of zero, so the accumulation is never called with an
/// empty loop to run.
template <typename F>
inline auto
accumulate_primitives(const CBasisFunction &bra, const CBasisFunction &ket, const std::vector<size_t> &dimensions,
                      const F &accumulate) -> void
{
    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    errors::assertMsgCritical(dimensions.size() == nprim_a * nprim_b,
                              std::string("SimdPrimitives.accumulate_primitives: Dimensions do not match the pairs of primitives"));

    for (size_t i = 0; i < nprim_a; i++)
    {
        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = dimensions[i * nprim_b + j];

            if (ncols == 0) continue;

            accumulate(CPrimitivePair{a_exps[i], b_exps[j], a_norms[i], b_norms[j], ncols});
        }
    }
}

}  // namespace simdfunc

#endif /* SimdPrimitives_hpp */
