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

#ifndef SimdTypes_hpp
#define SimdTypes_hpp

#include <array>
#include <cstddef>
#include <vector>

#include "Point.hpp"

constexpr size_t simd_width = 32;

using TDoubleArray = std::array<double, simd_width>;

using TFloatArray = std::array<double, simd_width>;

template <size_t N>
using TDoubleArray2D = std::array<TDoubleArray, N>;

namespace simd {  // simd namespace

/**
 Loads coordinates into three arrays (x, y, z).

 @param coords_x the array with Cartesian X coordinates.
 @param coords_y the array with Cartesian Y coordinates.
 @param coords_z the array with Cartesian Z coordinates.
 @param coords_r the vector of Cartesian coordinates.
 @param first the index of the range [first, last) to load.
 @param last the index of the range [first, last) to load.
 */
inline auto
loadCoordinates(TDoubleArray&                coords_x,
                TDoubleArray&                coords_y,
                TDoubleArray&                coords_z,
                const std::vector<TPoint3D>& coords_r,
                const int64_t                first,
                const int64_t                last) -> void
{
    for (int64_t i = first; i < last; i++)
    {
        const auto [x, y, z] = coords_r[i];

        const auto index = i - first;

        coords_x[index] = x;

        coords_y[index] = y;

        coords_z[index] = z;
    }
}

/**
 Loads primitive GTOs data from vector of primtive GTOs.

 @param data the array loaded with primitive GTOs data.
 @param prim_data the vector of primitive GTOs data.
 @param ipgto the index of leading primitive GTO.
 @param ncgtos the number of contracted GTOs.
 @param first the index of the range [first, last) to load.
 @param last the index of the range [first, last) to load.
 */
inline auto
loadPrimitiveGTOsData(TDoubleArray&              data,
                      const std::vector<double>& prim_data,
                      const int64_t              ipgto,
                      const int64_t              ncgtos,
                      const int64_t              first,
                      const int64_t              last) -> void
{
    for (int64_t i = first; i < last; i++)
    {
        data[i - first] = prim_data[ipgto * ncgtos + i];
    }
}

/**
 Sets all elements of array to zero.

 @param data the array to set to zero.
 */
inline auto
zero(TDoubleArray& data) -> void
{
    auto ptr_data = data.data();

#pragma omp simd aligned(ptr_data : 64)
    for (size_t i = 0; i < simd_width; i++)
    {
        ptr_data[i] = 0.0;
    }
}

}  // namespace simd

#endif /* SimdTypes_hpp */
