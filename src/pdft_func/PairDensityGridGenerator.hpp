//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#ifndef PairDensityGridGenerator_hpp
#define PairDensityGridGenerator_hpp

#include "DenseMatrix.hpp"
#include "MultiTimer.hpp"

namespace pairdengridgen {  // pairdengridgen namespace

/**
 Generates density and on-top pair-density for LDA.

 @param rho the pointer to density.
 @param gtoValues the GTO values on grid points.
 @param densityMatrix the total density matrix.
 @param activeMOs the MO coefficients.
 @param twoBodyDensityMatrix the MO two-body density matrix.
 @param timer the timer.
 */
void generatePairDensityForLDA(double*               rho,
                               const CDenseMatrix&   gtoValues,
                               const CDenseMatrix&   densityMatrix,
                               const CDenseMatrix&   activeMOs,
                               const CDenseMatrix&   twoBodyDensityMatrix,
                               CMultiTimer&          timer);

/**
 Generates density and on-top pair-density for GGA.

 @param rho the pointer to density.
 @param rhograd the pointer to density gradient.
 @param sigma the pointer to dot product of density gradient.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO X derivative values on grid points.
 @param gtoValuesY the GTO Y derivative values on grid points.
 @param gtoValuesZ the GTO Z derivative values on grid points.
 @param densityMatrix the total density matrix.
 @param activeMOs the MO coefficients.
 @param twoBodyDensityMatrix the MO two-body density matrix.
 @param timer the timer.
 */

void generatePairDensityForGGA(double*               rho,
                               double*               rhograd,
                               double*               sigma,
                               const CDenseMatrix&   gtoValues,
                               const CDenseMatrix&   gtoValuesX,
                               const CDenseMatrix&   gtoValuesY,
                               const CDenseMatrix&   gtoValuesZ,
                               const CDenseMatrix&   densityMatrix,
                               const CDenseMatrix&   activeMOs,
                               const CDenseMatrix&   twoBodyDensityMatrix,
                               CMultiTimer&          timer);

}  // namespace pairdengridgen

#endif /* PairDensityGridGenerator_hpp */
