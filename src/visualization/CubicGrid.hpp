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

#ifndef CubicGrid_hpp
#define CubicGrid_hpp

#include <array>
#include <cstdint>
#include <vector>

/**
 Class CCubicGrid contains cubic grid points.
 */
class CCubicGrid
{
    /**
     Origin of the cubic grid points.
     */
    std::array<double, 3> _origin;

    /**
     Step size in three dimensions.
     */
    std::array<double, 3> _stepSize;

    /**
     Number of points in three dimensions.
     */
    std::array<int, 3> _numPoints;

    /**
     Values at the grid points.
     */
    std::vector<double> _values;

   public:
    /**
     Creates an empty cubic grid object.
     */
    CCubicGrid();

    /**
     Creates a cubic grid object from origin, step size, and
     number of points.
     */
    CCubicGrid(const std::array<double, 3>& origin, const std::array<double, 3>& stepSize, const std::array<int, 3>& numPoints);

    /**
     Creates a cubic grid object by copying other cubic grid object.

     @param source the cubic grid object.
     */
    CCubicGrid(const CCubicGrid& source);

    /**
     Creates a cubic grid object by moving other cubic grid object.

     @param source the cubic grid object.
     */
    CCubicGrid(CCubicGrid&& source) noexcept;

    /**
     Destroys a cubic grid object.
     */
    ~CCubicGrid();

    /**
     Assigns a cubic grid object by copying other cubic grid object.

     @param source the cubic grid object.
     */
    CCubicGrid& operator=(const CCubicGrid& source);

    /**
     Assigns a cubic grid object by moving other cubic grid object.

     @param source the cubic grid object.
     */
    CCubicGrid& operator=(CCubicGrid&& source) noexcept;

    /**
     Compares cubic grid object with other cubic grid object.

     @param other the cubic grid object.
     @return true if cubic grid objects are equal, false otherwise.
     */
    bool operator==(const CCubicGrid& other) const;

    /**
     Compares cubic grid object with other cubic grid object.

     @param other the cubic grid object.
     @return true if cubic grid objects are not equal, false otherwise.
     */
    bool operator!=(const CCubicGrid& other) const;

    /**
     Gets coordinate of the origin.

     @return coordinate of the origin.
     */
    std::array<double, 3> getOrigin() const;

    /**
     Gets step size in X, Y and Z direction.

     @return step size in X, Y and Z direction.
     */
    std::array<double, 3> getStepSize() const;

    /**
     Gets number of points in X, Y and Z direction.

     @return number of points in X, Y and Z direction.
     */
    std::array<int, 3> getNumPoints() const;

    /**
     Gets constant pointer to first element of cubic grid values.

     @return the constant pointer to first element of cubic grid values.
     */
    const double* values() const;

    /**
     Gets pointer to first element of cubic grid values.

     @return the pointer to first element of cubic grid values.
     */
    double* values();

    /**
     Sets the cubic grid values.

     @param vals the values on the grid points.
     */
    void setValues(const std::vector<double> vals);
};

#endif /* CubicGrid_hpp */
