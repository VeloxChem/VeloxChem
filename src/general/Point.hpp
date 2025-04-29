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

#ifndef Point_hpp
#define Point_hpp

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>

#include "CustomConstrains.hpp"
#include "MathFunc.hpp"

/// @brief Templated Cartesian point class.
/// @tparam T The Cartesian coordinates storage type.
template <FloatingPoint T>
class TPoint
{
   public:
    /// @brief The default constructor.
    TPoint() {};

    /// @brief The constructor with coordinates array.
    /// @param coordinates The coordinates array to create point.
    TPoint(const std::array<T, 3> &coordinates) : _coordinates(coordinates) {};

    /// @brief The default copy constructor.
    /// @param other The point to be copied.
    TPoint(const TPoint<T> &other) = default;

    /// @brief The default move constructor.
    /// @param other The point to be moved.
    TPoint(TPoint<T> &&other) noexcept = default;

    /// @brief The default destructor.
    ~TPoint() = default;

    /// @brief The default copy assignment operator.
    /// @param other The point to be copy assigned.
    /// @return The assigned point.
    auto operator=(const TPoint<T> &other) -> TPoint<T> & = default;

    /// @brief The default move assignment operator.
    /// @param other The point to be move assigned.
    /// @return The assigned point.
    auto operator=(TPoint<T> &&other) noexcept -> TPoint<T> & = default;

    /// @brief The equality operator.
    /// @param other The point to be compared.
    /// @return True if points are equal, False otherwise.
    auto
    operator==(const TPoint<T> &other) const -> bool
    {
        return std::ranges::equal(
            _coordinates, other._coordinates, [](auto lhs, auto rhs) -> bool { return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12); });
    }

    /// @brief The non-equality operator.
    /// @param other The point to be compared.
    /// @return True if points are not equal, False otherwise.
    auto
    operator!=(const TPoint<T> &other) const -> bool
    {
        return !(*this == other);
    }

    /// @brief The getter for Cartesian coordinates of point.
    /// @return The Cartesian coordinates of point.
    auto
    coordinates() const -> std::array<T, 3>
    {
        return _coordinates;
    };

    /// @brief Scales Cartesian coordinates of point by given factor.
    /// @param factor The factor to scale coordinates.
    auto
    scale(const T &factor) -> void
    {
        std::ranges::transform(_coordinates, _coordinates.begin(), [=](const T &val) { return factor * val; });
    }

    /// @brief Computes square of length for vector defined by this point.
    /// @return The square of length of vector given by point.
    auto
    length_square() const -> T
    {
        return std::accumulate(_coordinates.begin(), _coordinates.end(), T{0.0}, [](const T &sum, const T &val) { return sum + val * val; });
    }

    /// @brief Computes length for vector defined by this point.
    /// @return The length of vector given by point.
    auto
    length() const -> T
    {
        return std::sqrt(length_square());
    }

    /// @brief Computes square of distance between this point and given other
    /// point.
    /// @param other The other point to compute square of distance.
    /// @return The square of distance between two points.
    auto
    distance_square(const TPoint<T> &other) const -> T
    {
        return std::inner_product(
            _coordinates.begin(),
            _coordinates.end(),
            other._coordinates.begin(),
            T{0.0},
            [](const T &sum, const T &val) { return sum + val * val; },
            [](const T &val_a, const T &val_b) { return val_b - val_a; });
    }

    /// @brief Computes distance between this point and given other point.
    /// @param other The other point to compute distance.
    /// @return The distance between two points.
    auto
    distance(const TPoint<T> &other) const -> T
    {
        return std::sqrt(distance_square(other));
    }

   private:
    /// @brief The Cartesian coordinates of point.
    std::array<T, 3> _coordinates;
};

#endif /* Point_hpp */
