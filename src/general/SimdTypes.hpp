//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

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
