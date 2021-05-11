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

#ifndef LebedevLaikovQuadrature_hpp
#define LebedevLaikovQuadrature_hpp

//===================================================//
// Adaptation of the original C code written         //
// by Dmitri Laikov.                                 //
//===================================================//
// Reference:                                        //
//===================================================//
// Vyacheslav Lebedev, Dmitri Laikov,                //
// A quadrature formula for the sphere of the 131st  //
// algebraic order of accuracy,                      //
// Russian Academy of Sciences Doklady Mathematics,  //
// Volume 59, Number 3, 1999, pages 477-481.         //
//===================================================//

#include <cstdint>
#include <cstdlib>

#include "MemBlock2D.hpp"

/**
 Class CLebedevLaikovQuadrature class generates Lebedev-Laikov quadrature.

 @author Z. Rinkevicius
 */
class CLebedevLaikovQuadrature
{
    /**
     The number of angular points.
     */
    int32_t _nAngularPoints;

    /**
     Generates angular Lebedev-Laikov quadrature with 6 points.

     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> _generateQuadratureWith6Points() const;

    /**
     Generates angular Lebedev-Laikov quadrature with 50 points.

     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> _generateQuadratureWith50Points() const;

    /**
     Generates angular Lebedev-Laikov quadrature with 110 points.

     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> _generateQuadratureWith110Points() const;

    /**
     Generates angular Lebedev-Laikov quadrature with 194 points.

     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> _generateQuadratureWith194Points() const;

    /**
     Generates angular Lebedev-Laikov quadrature with 302 points.

     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> _generateQuadratureWith302Points() const;

    /**
     Generates angular Lebedev-Laikov quadrature with 434 points.

     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> _generateQuadratureWith434Points() const;

    /**
     Generates angular Lebedev-Laikov quadrature with 590 points.

     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> _generateQuadratureWith590Points() const;

    /**
     Generates angular Lebedev-Laikov quadrature with 770 points.

     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> _generateQuadratureWith770Points() const;

    /**
     Generates angular Lebedev-Laikov quadrature with 974 points.

     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> _generateQuadratureWith974Points() const;

    /**
     Generates angular Lebedev-Laikov quadrature with 2030 points.

     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> _generateQuadratureWith2030Points() const;

    /**
     Generates 6 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param weight the weight assigned to grid weights.
     */
    void _generateCaseOne(CMemBlock2D<double>& gridPoints, const int32_t offset, const double weight) const;

    /**
     Generates 12 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param weight the weight assigned to grid weights.
     */
    void _generateCaseTwo(CMemBlock2D<double>& gridPoints, const int32_t offset, const double weight) const;
    /**
     Generates 8 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param weight the weight assigned to grid weights.
     */
    void _generateCaseThree(CMemBlock2D<double>& gridPoints, const int32_t offset, const double weight) const;

    /**
     Generates first 24 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param factor the scaling factor of grid coordinates.
     @param weight the weight assigned to grid weights.
     */
    void _generateCaseFour(CMemBlock2D<double>& gridPoints, const int32_t offset, const double factor, const double weight) const;

    /**
     Generates second 24 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param factor the scaling factor of grid coordinates.
     @param weight the weight assigned to grid weights.
     */
    void _generateCaseFive(CMemBlock2D<double>& gridPoints, const int32_t offset, const double factor, const double weight) const;

    /**
     Generates 48 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param factorA the first scaling factor of grid coordinates.
     @param factorB the second scaling factor of grid coordinates.
     @param weight the weight assigned to grid weights.
     */
    void _generateCaseSix(CMemBlock2D<double>& gridPoints,
                          const int32_t        offset,
                          const double         factorA,
                          const double         factorB,
                          const double         weight) const;

   public:
    /**
     Creates a Lebedev-Laikov quadrature object.

     @param nAngularPoints the number of angular points.
     */
    CLebedevLaikovQuadrature(const int32_t nAngularPoints);

    /**
     Destroys a Lebedev-Laikov quadrature object.
     */
    ~CLebedevLaikovQuadrature();

    /**
     Generates quadrature points for angular Lebedev-Laikov quadrature.

     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> generate() const;
};

#endif /* LebedevLaikovQuadrature_hpp */
