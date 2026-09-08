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

#include "ErrorHandler.hpp"
#include "SimdAlign.hpp"
#include "SimdMatrix.hpp"

namespace simdfunc {  // simdfunc namespace

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

/// @brief Computes the displacement of the Gaussian product center from the atom
/// on bra side, for one pair of primitives.
/// @param buffer The buffer of the pair of primitives.
/// @param coordinates The coordinates of the atom pairs, whose rows six to eight
/// hold the vector from the atom on bra side to the atom on ket side.
/// @param target The first of the three rows of the buffer to write.
/// @param ncols The number of atom pairs the pair of primitives reaches.
/// @param fa The displacement in units of the vector between the atoms, which is
/// -b / (a + b) for the atom on bra side.
/// @note The displacement is P - A = -b / (a + b) times the vector from the atom on
/// bra side to the atom on ket side, as the Gaussian product center P divides that
/// vector in the ratio of the exponents. The factor is formed by the caller, which
/// holds the exponents of the pair of primitives.
/// @note The rows of both matrices start at a cache line boundary, so the loop is
/// vectorized with aligned loads and stores.
inline auto
compute_pa(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target, const size_t ncols, const double fa) -> void
{
    auto *pa_x = buffer.data(target + 0);
    auto *pa_y = buffer.data(target + 1);
    auto *pa_z = buffer.data(target + 2);

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

#pragma omp simd aligned(pa_x, pa_y, pa_z, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        pa_x[k] = fa * ab_x[k];

        pa_y[k] = fa * ab_y[k];

        pa_z[k] = fa * ab_z[k];
    }
}

/// @brief Computes the displacement of the Gaussian product center from the atom
/// on ket side, for one pair of primitives.
/// @param buffer The buffer of the pair of primitives.
/// @param coordinates The coordinates of the atom pairs, whose rows six to eight
/// hold the vector from the atom on bra side to the atom on ket side.
/// @param target The first of the three rows of the buffer to write.
/// @param ncols The number of atom pairs the pair of primitives reaches.
/// @param fb The displacement in units of the vector between the atoms, which is
/// a / (a + b) for the atom on ket side.
/// @note The displacement is P - B = a / (a + b) times the vector from the atom on
/// bra side to the atom on ket side. It differs from the one on bra side in the
/// factor alone, which is why the two are separate routines rather than one with a
/// side to select.
/// @note The rows of both matrices start at a cache line boundary, so the loop is
/// vectorized with aligned loads and stores.
inline auto
compute_pb(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target, const size_t ncols, const double fb) -> void
{
    auto *pb_x = buffer.data(target + 0);
    auto *pb_y = buffer.data(target + 1);
    auto *pb_z = buffer.data(target + 2);

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

#pragma omp simd aligned(pb_x, pb_y, pb_z, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        pb_x[k] = fb * ab_x[k];

        pb_y[k] = fb * ab_y[k];

        pb_z[k] = fb * ab_z[k];
    }
}

}  // namespace simdfunc

#endif /* SimdPrimitives_hpp */
