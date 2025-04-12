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
    auto getNumberOfGridPoints() const -> int64_t;

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
