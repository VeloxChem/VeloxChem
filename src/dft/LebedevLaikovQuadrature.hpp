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

#include "DenseMatrix.hpp"

/**
 Class CLebedevLaikovQuadrature class generates Lebedev-Laikov quadrature.

 @author Z. Rinkevicius
 */
class CLebedevLaikovQuadrature
{
    /**
     The number of angular points.
     */
    int64_t _nAngularPoints;

    /**
     Generates angular Lebedev-Laikov quadrature with 6 points.

     @return the quadrature points (coordinates, weights).
     */
    auto _generateQuadratureWith6Points() const -> CDenseMatrix;

    /**
     Generates angular Lebedev-Laikov quadrature with 50 points.

     @return the quadrature points (coordinates, weights).
     */
    auto _generateQuadratureWith50Points() const -> CDenseMatrix;

    /**
     Generates angular Lebedev-Laikov quadrature with 110 points.

     @return the quadrature points (coordinates, weights).
     */
    auto _generateQuadratureWith110Points() const -> CDenseMatrix;

    /**
     Generates angular Lebedev-Laikov quadrature with 194 points.

     @return the quadrature points (coordinates, weights).
     */
    auto _generateQuadratureWith194Points() const -> CDenseMatrix;

    /**
     Generates angular Lebedev-Laikov quadrature with 302 points.

     @return the quadrature points (coordinates, weights).
     */
    auto _generateQuadratureWith302Points() const -> CDenseMatrix;

    /**
     Generates angular Lebedev-Laikov quadrature with 434 points.

     @return the quadrature points (coordinates, weights).
     */
    auto _generateQuadratureWith434Points() const -> CDenseMatrix;

    /**
     Generates angular Lebedev-Laikov quadrature with 590 points.

     @return the quadrature points (coordinates, weights).
     */
    auto _generateQuadratureWith590Points() const -> CDenseMatrix;

    /**
     Generates angular Lebedev-Laikov quadrature with 770 points.

     @return the quadrature points (coordinates, weights).
     */
    auto _generateQuadratureWith770Points() const -> CDenseMatrix;

    /**
     Generates angular Lebedev-Laikov quadrature with 974 points.

     @return the quadrature points (coordinates, weights).
     */
    auto _generateQuadratureWith974Points() const -> CDenseMatrix;

    /**
     Generates angular Lebedev-Laikov quadrature with 2030 points.

     @return the quadrature points (coordinates, weights).
     */
    auto _generateQuadratureWith2030Points() const -> CDenseMatrix;

    /**
     Generates 6 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param weight the weight assigned to grid weights.
     */
    auto _generateCaseOne(CDenseMatrix& gridPoints, const int64_t offset, const double weight) const -> void;

    /**
     Generates 12 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param weight the weight assigned to grid weights.
     */
    auto _generateCaseTwo(CDenseMatrix& gridPoints, const int64_t offset, const double weight) const -> void;
    /**
     Generates 8 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param weight the weight assigned to grid weights.
     */
    auto _generateCaseThree(CDenseMatrix& gridPoints, const int64_t offset, const double weight) const -> void;

    /**
     Generates first 24 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param factor the scaling factor of grid coordinates.
     @param weight the weight assigned to grid weights.
     */
    auto _generateCaseFour(CDenseMatrix& gridPoints, const int64_t offset, const double factor, const double weight) const -> void;

    /**
     Generates second 24 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param factor the scaling factor of grid coordinates.
     @param weight the weight assigned to grid weights.
     */
    auto _generateCaseFive(CDenseMatrix& gridPoints, const int64_t offset, const double factor, const double weight) const -> void;

    /**
     Generates 48 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param factorA the first scaling factor of grid coordinates.
     @param factorB the second scaling factor of grid coordinates.
     @param weight the weight assigned to grid weights.
     */
    auto _generateCaseSix(CDenseMatrix& gridPoints, const int64_t offset, const double factorA, const double factorB, const double weight) const
        -> void;

   public:
    /**
     Creates a Lebedev-Laikov quadrature object.

     @param nAngularPoints the number of angular points.
     */
    CLebedevLaikovQuadrature(const int64_t nAngularPoints);

    /**
     Destroys a Lebedev-Laikov quadrature object.
     */
    ~CLebedevLaikovQuadrature();

    /**
     Generates quadrature points for angular Lebedev-Laikov quadrature.

     @return the quadrature points (coordinates, weights).
     */
    auto generate() const -> CDenseMatrix;
};

#endif /* LebedevLaikovQuadrature_hpp */
