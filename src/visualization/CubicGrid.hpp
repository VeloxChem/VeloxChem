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

#ifndef CubicGrid_hpp
#define CubicGrid_hpp

#include <cstdint>
#include <vector>

#include "MemBlock.hpp"

/**
 Class CCubicGrid contains cubic grid points.

 @author X. Li
 */
class CCubicGrid
{
    /**
     Origin of the cubic grid points.
     */
    std::vector<double> _origin;

    /**
     Step size in three dimensions.
     */
    std::vector<double> _stepSize;

    /**
     Number of points in three dimensions.
     */
    std::vector<int32_t> _numPoints;

    /**
     Values at the grid points.
     */
    CMemBlock<double> _values;

   public:
    /**
     Creates an empty cubic grid object.
     */
    CCubicGrid();

    /**
     Creates a cubic grid object from origin, step size, and
     number of points.
     */
    CCubicGrid(const std::vector<double>& origin, const std::vector<double>& stepSize, const std::vector<int32_t>& numPoints);

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
    std::vector<double> getOrigin() const;

    /**
     Gets step size in X, Y and Z direction.

     @return step size in X, Y and Z direction.
     */
    std::vector<double> getStepSize() const;

    /**
     Gets number of points in X, Y and Z direction.

     @return number of points in X, Y and Z direction.
     */
    std::vector<int32_t> getNumPoints() const;

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
