//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
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

#ifndef GridBox_hpp
#define GridBox_hpp

#include <array>

#include "DenseMatrix.hpp"

/**
 Class CGridBox contains DFT grid points.

 @author X. Li
 */
class CGridBox
{
    /**
     Dimension of the box (xmin,ymin,zmin,xmax,ymax,zmax).
     */
    std::array<double, 6> _dimension;

    /**
     DFT grid points (x,y,z,weight).
     */
    CDenseMatrix _points;

   public:
    /**
     Creates an empty grid box object.
     */
    CGridBox();

    /**
     Creates a grid box object from dimension and points.
     */
    CGridBox(const std::array<double, 6>& dimension, const CDenseMatrix& points);

    /**
     Creates a grid box object by copying other grid box object.

     @param source the grid box object.
     */
    CGridBox(const CGridBox& source);

    /**
     Creates a grid box object by moving other grid box object.

     @param source the grid box object.
     */
    CGridBox(CGridBox&& source) noexcept;

    /**
     Destroys a grid box object.
     */
    ~CGridBox();

    /**
     Assigns a grid box object by copying other grid box object.

     @param source the grid box object.
     */
    auto operator=(const CGridBox& source) -> CGridBox&;

    /**
     Assigns a grid box object by moving other grid box object.

     @param source the grid box object.
     */
    auto operator=(CGridBox&& source) noexcept -> CGridBox&;

    /**
     Compares grid box object with other grid box object.

     @param other the grid box object.
     @return true if grid box objects are equal, false otherwise.
     */
    auto operator==(const CGridBox& other) const -> bool;

    /**
     Compares grid box object with other grid box object.

     @param other the grid box object.
     @return true if grid box objects are not equal, false otherwise.
     */
    auto operator!=(const CGridBox& other) const -> bool;

    /**
     Gets the box dimension.

     @return the box dimension
     */
    auto getBoxDimension() const -> std::array<double, 6>;

    /**
     Gets the grid points.

     @return the grid points.
     */
    auto getGridPoints() const -> CDenseMatrix;

    /**
     Gets number of grid points.

     @return the number of grid points.
     */
    auto getNumberOfGridPoints() const -> int32_t;

    /**
     Gets Cartesian X coordinates of grid points.

     @return the constant pointer to Cartesian X coordinates of grid points.
     */
    auto getCoordinatesX() const -> const double*;

    /**
     Gets Cartesian Y coordinates of grid points.

     @return the constant  pointer to Cartesian Y coordinates of grid points.
     */
    auto getCoordinatesY() const -> const double*;

    /**
     Gets Cartesian Z coordinates of grid points.

     @return the constant pointer to Cartesian Z coordinates of grid points.
     */
    auto getCoordinatesZ() const -> const double*;

    /**
     Gets weights of grid points.

     @return the constant pointer to weights of grid points.
     */
    auto getWeights() const -> const double*;
};

#endif /* GridBox_hpp */
