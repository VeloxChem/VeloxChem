//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <type_traits>


namespace mathfunc {  // mathfunc namespace

/**
 * Sets all elements of vector to zero.
 *
 * @tparam T scalar type
 * @param vector the vector.
 * @param nElements the number of elements in vector.
 */
template <typename T>
inline auto
zero(T* vector, const int32_t nElements) -> decltype((void)(std::is_arithmetic_v<T>), void())
{
#pragma omp simd aligned(vector : VLX_ALIGN)
    for (int32_t i = 0; i < nElements; i++)
        vector[i] = T{0};
}

/**
 * Sets all elements of vector to specific value.
 *
 * @tparam T scalar type
 * @param vector the vector.
 * @param value the value of element.
 * @param nElements the number of elements in vector.
 */
template <typename T>
inline auto
set_to(T* vector, const T value, const int32_t nElements) -> decltype((void)(std::is_arithmetic_v<T>), void())
{
#pragma omp simd aligned(vector : VLX_ALIGN)
    for (int32_t i = 0; i < nElements; i++)
        vector[i] = value;
}

/**
 * Computes sum of all elements vector.
 *
 * @tparam T scalar type
 * @param vector the vector.
 * @param nElements the number of elements in vector.
 * @return sum of all elements in vector.
 */
template <typename T>
inline auto
sum(const T* vector, const int32_t nElements) -> decltype((void)(std::is_arithmetic_v<T>), T())
{
    auto fsum = T{0};

#pragma omp simd aligned(vector : VLX_ALIGN)
    for (int32_t i = 0; i < nElements; i++)
        fsum += vector[i];

    return fsum;
}

/**
 Scales all elements of real numbers vector by specific factor.

 @param vector the vector of real numbers.
 @param factor the scaling factor.
 @param nElements the number of elements in vector.
 */
void scale(double* vector, const double factor, const int32_t nElements);

/**
 Adds vector scaled by factor to other vector i.e. va = va + f * vb.

 @param aVector the vector of real numbers.
 @param bVector the vector of real numbers.
 @param factor the scaling factor.
 @param nElements the number of elements in vector.
 */
void add_scaled(double* aVector, const double* bVector, const double factor, const int32_t nElements);

/**
 * Determines largest element in vector.
 *
 * @tparam T scalar type
 * @param vector the vector.
 * @param nElements the number of elements in vector.
 * @return the largest element in real numbers vector.
 */
template <typename T>
inline auto
max(const T* vector, const int32_t nElements) -> decltype((void)(std::is_arithmetic_v<T>), T())
{
    auto fmax = vector[0];

    for (int32_t i = 1; i < nElements; i++)
    {
        auto cmax = vector[i];

        if (cmax > fmax) fmax = cmax;
    }

    return fmax;
}

/**
 Normalizes vector of real numbers.

 @param vector the vector of real numbers.
 @param nElements the number of elements in vector.
 */
void normalize(double* vector, const int32_t nElements);

/**
 Sets indexes vector using size vector.

 @param aVector the indexes vector.
 @param bVector the sizes vector.
 @param nElements the number of elements in vectors.
 */
void indexes(int32_t* aVector, const int32_t* bVector, const int32_t nElements);

/**
 Sets indexes vector using size vector and offset.

 @param aVector the indexes vector.
 @param bVector the sizes vector.
 @param offset the offset of first index.
 @param nElements the number of elements in vectors.
 */
void indexes(int32_t* aVector, const int32_t* bVector, const int32_t offset, const int32_t nElements);

/**
 Sets ordering vector for given vector of binary values (0 or 1) by storing
 all indexes of binary vector elements equal to 1.

 @param aVector the ordering vector.
 @param bVector the binary vector.
 @param nElements the number of elements in vectors.
 */
void ordering(int32_t* aVector, const int32_t* bVector, const int32_t nElements);

/**
 Computes distance between two 3D vectors.

 @param aCoordX the Cartesian X coordinate of first vector.
 @param aCoordY the Cartesian Y coordinate of first vector.
 @param aCoordZ the Cartesian Z coordinate of first vector.
 @param bCoordX the Cartesian X coordinate of second vector.
 @param bCoordY the Cartesian Y coordinate of second vector.
 @param bCoordZ the Cartesian Z coordinate of second vector.
 @return the distance between vectors.
 */
inline double
distance(const double aCoordX, const double aCoordY, const double aCoordZ, const double bCoordX, const double bCoordY, const double bCoordZ)
{
    auto rx = aCoordX - bCoordX;

    auto ry = aCoordY - bCoordY;

    auto rz = aCoordZ - bCoordZ;

    return std::sqrt(rx * rx + ry * ry + rz * rz);
}

/**
 Computes distances between reference point A and vector of B points.

 @param abDistancesX the vector of distances R(AB)_x = A_x - B_x.
 @param abDistancesY the vector of distances R(AB)_y = A_y - B_y.
 @param abDistancesZ the vector of distances R(AB)_z = A_z - B_z.
 @param aCoordX the Cartesian X coordinate of point A.
 @param aCoordY the Cartesian Y coordinate of point A.
 @param aCoordZ the Cartesian Z coordinate of point A.
 @param bCoordsX the vector of Cartesian X coordinates of points B.
 @param bCoordsY the vector of Cartesian Y coordinates of points B.
 @param bCoordsZ the vector of Cartesian Z coordinates of points B.
 @param nElements the number of points B.
 */
void distances(double*       abDistancesX,
               double*       abDistancesY,
               double*       abDistancesZ,
               const double  aCoordX,
               const double  aCoordY,
               const double  aCoordZ,
               const double* bCoordsX,
               const double* bCoordsY,
               const double* bCoordsZ,
               const int32_t nElements);

/**
 Computes Chebtshev quadrature of second kind in [-1,1] interval.

 @param coordinates the vector of quadature coordinates.
 @param weights the vector of quadrature weights.
 @param nPoints the number of points in quadrature.
 */
void quadChebyshevOfKindTwo(double* coordinates, double* weights, const int32_t nPoints);

/**
 * Copies scalars from one vector to another vector.
 *
 * @tparam T scalar type.
 * @param aVector the destination vector.
 * @param aPosition the position of first copied element in destination vector.
 * @param bVector the source vector.
 * @param bPosition the position of first copied element in source vector.
 * @param nElements the number of elements.
 */
template <typename T>
inline auto
copy(T* aVector, const int32_t aPosition, const T* bVector, const int32_t bPosition, const int32_t nElements)
    -> decltype((void)(std::is_arithmetic_v<T>), void())
{
    std::copy_n(bVector + bPosition, nElements, aVector + aPosition);
}

/**
 Determines maximum number of components for tensor of given order.

 @param order the order of tensor.
 @return the number of components.
 */
int32_t maxTensorComponents(const int32_t order);

/** Fill raw array with random numbers in interval.
 *
 * @tparam T scalar type of raw array.
 * @param[in,out] dst raw array.
 * @param[in] lower lower bound of interval.
 * @param[in] upper upper bound of interval.
 * @param[in] sz number of elements in array.
 *
 * This function uses the C++ default random engine with random seeding.
 */
template <typename T>
auto
fill_random(T* dst, T lower, T upper, size_t sz) -> void
{
    static_assert(std::is_arithmetic_v<T>, "Scalar type must be arithmetic.");

    // random number generator
    auto gen = std::default_random_engine(std::random_device()());

    // distribution (use IIFE idiom to get the right distribution at compile-time)
    auto dist = [lower, upper] {
        if constexpr (std::is_floating_point_v<T>)
        {
            return std::uniform_real_distribution<T>(lower, upper);
        }
        else
        {
            return std::uniform_int_distribution<T>(lower, upper);
        }
    }();

    std::generate(dst, dst + sz, [&dist, &gen]() { return dist(gen); });
}
}  // namespace mathfunc

#endif /* MathFunc_hpp */
