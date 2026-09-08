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

// NOTE: the generated kernels reach the Boys function through this header, as
// they include it for the other helpers of simdfunc and do not include the Boys
// header of their own. The include belongs in the kernels which call it and can
// go once they carry it.
#include "SimdBoysFunc.hpp"

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

/// @brief Adds the contribution of one pair of primitives to the contracted rows
/// of the buffer.
/// @param buffer The buffer of the combination of basis functions.
/// @param target The first of the contracted rows to add to.
/// @param source The first of the rows of the pair of primitives to add.
/// @param nrows The number of rows to add.
/// @param ncols The number of atom pairs the pair of primitives reaches.
/// @note The rows of a pair of primitives are written over for every pair, while
/// the contracted rows accumulate across them, which is why make_primitive_buffer
/// zeroes the buffer: the first pair adds to rows nothing has written.
/// @note The rows of the buffer start at a cache line boundary, so the loop is
/// vectorized with aligned loads and stores.
inline auto
contract_primitives(CSimdMatrix &buffer, const size_t target, const size_t source, const size_t nrows, const size_t ncols) -> void
{
    for (size_t r = 0; r < nrows; r++)
    {
        auto *dst = buffer.data(target + r);

        const auto *src = buffer.data(source + r);

#pragma omp simd aligned(dst, src : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            dst[k] += src[k];
        }
    }
}

/// @brief Prepares a buffer a caller holds for one combination of basis functions.
/// @param buffer The view of the arena of the block, which is reshaped to this
/// combination.
/// @param nrows The number of rows this combination uses.
/// @param ncols The number of atom pairs this combination reaches.
/// @return The number of columns, which is ncols, so a caller reads it back in one
/// statement.
/// @note The arena belongs to the block and not to the combination, so nothing is
/// allocated here. The shape is taken per combination all the same: a combination
/// which reaches fewer atom pairs than the block holds wants its rows that much
/// closer together, and a view stretched to the pairs of the block would scatter
/// every row of it over a page of its own.
/// @note Only the rows this combination uses are zeroed. A combination which
/// accumulates with += into its contracted rows needs them to start at zero, which is
/// what this gives it; the rows beyond it hold whatever the combination before left
/// there and are none of its business.
/// @note Those rows are zeroed in one stroke, padding and all, and not row by row
/// over the columns alone. They are contiguous once the shape is taken, so this is a
/// single fill rather than one per row which stops short of the end of a cache line
/// and has the line read back to write its tail.
inline auto
prepare_buffer(CSimdMatrix &buffer, const size_t nrows, const size_t ncols) -> size_t
{
    buffer.reshape(nrows, ncols);

    auto *values = buffer.data();

    std::fill(values, values + nrows * buffer.pitch(), 0.0);

    return ncols;
}

/// @brief Prepares a buffer a caller holds for one combination of basis functions,
/// spanning the atom pairs the furthest reaching pair of primitives reaches.
/// @param buffer The buffer of the block.
/// @param nrows The number of rows this combination uses.
/// @param dimensions The number of atom pairs each pair of primitives reaches.
/// @return The number of columns to work over, zero when no pair of primitives
/// reaches any atom pair.
/// @note The columns are searched for rather than assumed, for the reason given for
/// make_primitive_buffer: the last pair of primitives is not always the furthest
/// reaching.
inline auto
prepare_buffer(CSimdMatrix &buffer, const size_t nrows, const std::vector<size_t> &dimensions) -> size_t
{
    if (dimensions.empty()) return 0;

    const auto nmax = *std::ranges::max_element(dimensions);

    if (nmax == 0) return 0;

    return prepare_buffer(buffer, nrows, nmax);
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

/// @brief Computes the displacement of the Gaussian product center of a pair of
/// primitives from one atom on the ket side of a three-center quantity.
/// @param buffer The buffer of the combination of basis functions.
/// @param coordinates The coordinates of the atom pairs, whose rows zero to two hold
/// the atom on bra side and whose rows six to eight hold the vector between the atoms.
/// @param c_coordinates The coordinates of the atoms on the ket side, as three rows.
/// @param target The first of the three rows of the buffer to write.
/// @param iatom The atom on the ket side, as its column of c_coordinates.
/// @param ncols The number of atom pairs the pair of primitives reaches.
/// @param fc The displacement of the product center from the atom on bra side in
/// units of the vector between the atoms, which is b / (a + b).
/// @note The Gaussian product center P of the pair sits at A - b / (a + b) times the
/// vector from A to B, so its displacement from the atom C on the ket side is
/// (A - C) - fc times that vector. The atom on the ket side is one point and the atom
/// pairs are many, so its coordinates are read once and the loop runs over the pairs.
/// @note The rows of the buffer and of the coordinates start at a cache line
/// boundary, so the loop is vectorized with aligned loads and stores. The coordinates
/// of the atom on the ket side are scalars and are not part of the clause.
inline auto
compute_pc(CSimdMatrix       &buffer,
           const CSimdMatrix &coordinates,
           const CSimdMatrix &c_coordinates,
           const size_t       target,
           const size_t       iatom,
           const size_t       ncols,
           const double       fc) -> void
{
    auto *pc_x = buffer.data(target + 0);
    auto *pc_y = buffer.data(target + 1);
    auto *pc_z = buffer.data(target + 2);

    const auto *a_x = coordinates.data(0);
    const auto *a_y = coordinates.data(1);
    const auto *a_z = coordinates.data(2);

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto c_x = c_coordinates.data(0)[iatom];
    const auto c_y = c_coordinates.data(1)[iatom];
    const auto c_z = c_coordinates.data(2)[iatom];

#pragma omp simd aligned(pc_x, pc_y, pc_z, a_x, a_y, a_z, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        pc_x[k] = (a_x[k] - c_x) - fc * ab_x[k];

        pc_y[k] = (a_y[k] - c_y) - fc * ab_y[k];

        pc_z[k] = (a_z[k] - c_z) - fc * ab_z[k];
    }
}

}  // namespace simdfunc

#endif /* SimdPrimitives_hpp */
