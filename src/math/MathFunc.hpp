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

#ifndef MathFunc_hpp
#define MathFunc_hpp

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include <ranges>

#include "CustomConstrains.hpp"
#include "MathConst.hpp"

namespace mathfunc {  // mathfunc

/// Computes distance between two points in 3D space.
inline auto
distance(const double rax, const double ray, const double raz, const double rbx, const double rby, const double rbz) -> double
{
    const auto ab_x = rax - rbx;
    const auto ab_y = ray - rby;
    const auto ab_z = raz - rbz;

    return std::sqrt(ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);
}

/// Computes distance between two points in 3D space.
/// - Parameter r_a: the first point coordinates.
/// - Parameter r_b: the first point coordinates.
inline auto
distance(const std::array<double, 3>& ra, const std::array<double, 3>& rb) -> double
{
    const auto ab_x = ra[0] - rb[0];

    const auto ab_y = ra[1] - rb[1];

    const auto ab_z = ra[2] - rb[2];

    return std::sqrt(ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);
}

/// @brief Compares two floating point values.
/// @tparam T The floating point type.
/// @param flhs The left-hand-side floating point value to compare.
/// @param frhs The rigth-hand-side floating point value to compare.
/// @param rtol The relative tollerance.
/// @param atol The absolute tollerance.
/// @return True if floating point values are equal, False otherwise.
template <FloatingPoint T>
inline auto
equal(const T flhs, const T frhs, const T rtol, const T atol) -> bool
{
    return std::fabs(flhs - frhs) < std::max(atol, rtol * std::max(std::fabs(flhs), std::fabs(frhs)));
}

/// @brief Gets upper triangular matrix linearized index (column wise scheme).
/// @param i The index of row in matrix.
/// @param j The index of collumn in matrix.
/// @return The linearized index.
template <Integral T>
inline auto
uplo_index(const T i, const T j) -> T
{
    return i + j * (j + 1) / 2;
}

/// @brief Counts number of elements in vector matching given selector.
/// @param values  The vector of values.
/// @param selector  The selector to march vector values.
/// @return The number of elements in vector matching given selector.
template <Integral T>
inline auto
count_elements_by_values(const std::vector<T>& values,
                         const T               selector) -> T
{
    return static_cast<T>(std::ranges::count(values, selector));
}

/// Computes Chebtshev quadrature of second kind in [-1,1] interval.
/// @param coordinates the vector of quadature coordinates.
/// @param weights the vector of quadrature weights.
/// @param nPoints the number of points in quadrature.
inline auto
quadChebyshevOfKindTwo(double* coordinates, double* weights, const int nPoints) -> void
{
    // prefactor
    auto fstep = mathconst::pi_value() / (static_cast<double>(nPoints) + 1.0);

    // loop over grid points
    for (int i = 1; i < nPoints + 1; i++)
    {
        auto farg = static_cast<double>(i) * fstep;

        coordinates[i - 1] = std::cos(farg);

        auto warg = std::sin(farg);

        weights[i - 1] = fstep * warg * warg;
    }
}

/// Determines batch sizes for parallelization.
/// @param nElements the size of data vector.
/// @param nodes the number of processing elements.
/// @return the sizes of data batches.
inline auto
batch_sizes(const int nElements, const int nodes) -> std::vector<int>
{
    int ave = nElements / nodes;

    int rem = nElements % nodes;

    std::vector<int> counts;

    for (int p = 0; p < nodes; p++)
    {
        counts.push_back((p < rem) ? (ave + 1) : ave);
    }

    return counts;
}

/// Determines batch offsets for parallelization.
/// @param nElements the size of data vector.
/// @param nodes the number of processing elements.
/// @return the offsets of data batches.
inline auto
batch_offsets(const int nElements, const int nodes) -> std::vector<int>
{
    auto counts = mathfunc::batch_sizes(nElements, nodes);

    std::vector<int> displs;

    int index = 0;

    for (int p = 0; p < nodes; p++)
    {
        displs.push_back(index);

        index += counts[p];
    }

    return displs;
}

/// Determines batch size for parallelization.
/// @param nElements the size of data vector.
/// @param rank the rank of processing element.
/// @param nodes the number of processing elements.
/// @return the size of data batch.
inline auto
batch_size(const int nElements, const int rank, const int nodes) -> int
{
    auto counts = mathfunc::batch_sizes(nElements, nodes);

    return counts[rank];
}

/// Determines batch offset for parallelization.
/// @param nElements the size of data vector.
/// @param rank the rank of processing element.
/// @param nodes the number of processing elements.
/// @return the offset of data batch.
inline auto
batch_offset(const int nElements, const int rank, const int nodes) -> int
{
    auto displs = mathfunc::batch_offsets(nElements, nodes);

    return displs[rank];
}

/// @brief Computes linerized upper triangular index (row major arragment)
/// @param irow  The row of square matrix.
/// @param icolumn  The column of square matrix.
/// @return The linearized index.
template <Integral T>
inline auto
uplo_rm_index(const T irow,
              const T icolumn,
              const T nrows) -> T
{
    return icolumn + irow * nrows - irow * (irow + 1) / 2;
}

}  // namespace mathfunc

#endif /* MathFunc_hpp */
