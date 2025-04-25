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
 */
class CLebedevLaikovQuadrature
{
    /**
     The number of angular points.
     */
    int _nAngularPoints;

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
    auto _generateCaseOne(CDenseMatrix& gridPoints, const int offset, const double weight) const -> void;

    /**
     Generates 12 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param weight the weight assigned to grid weights.
     */
    auto _generateCaseTwo(CDenseMatrix& gridPoints, const int offset, const double weight) const -> void;
    /**
     Generates 8 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param weight the weight assigned to grid weights.
     */
    auto _generateCaseThree(CDenseMatrix& gridPoints, const int offset, const double weight) const -> void;

    /**
     Generates first 24 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param factor the scaling factor of grid coordinates.
     @param weight the weight assigned to grid weights.
     */
    auto _generateCaseFour(CDenseMatrix& gridPoints, const int offset, const double factor, const double weight) const -> void;

    /**
     Generates second 24 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param factor the scaling factor of grid coordinates.
     @param weight the weight assigned to grid weights.
     */
    auto _generateCaseFive(CDenseMatrix& gridPoints, const int offset, const double factor, const double weight) const -> void;

    /**
     Generates 48 points term of Lebedev-Laikov quadrature.

     @param gridPoints the quadrature points.
     @param offset the offset of grid points in quadrature points.
     @param factorA the first scaling factor of grid coordinates.
     @param factorB the second scaling factor of grid coordinates.
     @param weight the weight assigned to grid weights.
     */
    auto _generateCaseSix(CDenseMatrix& gridPoints, const int offset, const double factorA, const double factorB, const double weight) const
        -> void;

   public:
    /**
     Creates a Lebedev-Laikov quadrature object.

     @param nAngularPoints the number of angular points.
     */
    CLebedevLaikovQuadrature(const int nAngularPoints);

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
