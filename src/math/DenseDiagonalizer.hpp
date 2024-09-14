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

#ifndef DenseDiagonalizer_hpp
#define DenseDiagonalizer_hpp

#include <vector>

#include "DenseMatrix.hpp"

/**
 Class CDenseDiagonalizer provides methods for diagonalization of dense real
 symmetric matrices and for manipulation with computed eigenvalues and
 eigenvectors.

 @author Z. Rinkevicius
 */
class CDenseDiagonalizer
{
    /**
     The state of dense diagonalizer object : true - no errors,
     false - otherwise.
     */
    bool _state;

    /**
     The availability of eigenvalues and eigenvectors: true - available,
     false - otherwise.
     */
    bool _isSolved;

    /**
     The temporary dense matrix used by diagonalization routine.
     */
    CDenseMatrix _matrix;

    /**
     The eigenvectors of dense matrix.
     */
    CDenseMatrix _eigenVectors;

    /**
     The eigenvalues of dense matrix.
     */
    std::vector<double> _eigenValues;

   public:
    /**
     Creates an dense diagonalizer object.
     */
    CDenseDiagonalizer();

    /**
     Diagonalizes dense matrix and stores eigenvalues/eigenvectors into internal
     data buffers.

     @param matrix the dense matrix.
     */
    auto diagonalize(const CDenseMatrix& matrix) -> void;

    /**
     Gets state of dense diagonalizer object.

     @return true if no errors, false otherwise.
     */
    auto getState() const -> bool;

    /**
     Gets eigenvectors of dense matrix.

     @return the eigenvectors of matrix.
     */
    auto getEigenVectors() const -> CDenseMatrix;

    /**
     Gets eigenvalues of dense matrix.

     @return the eigenvalues of matrix.
     */
    auto getEigenValues() const -> std::vector<double>;

    /**
     Computes A^-1/2 matrix using eigenvalues and eigenvectors of A matrix.

     @return the A^-1/2 matrix.
     */
    auto getInvertedSqrtMatrix() const -> CDenseMatrix;

    /**
     Computes A^-1 matrix using eigenvalues and eigenvectors of A matrix.

     @return the A^-1 matrix.
     */
    auto getInvertedMatrix() const -> CDenseMatrix;
};

#endif /* DenseDiagonalizer_hpp */
