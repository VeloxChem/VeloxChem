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

#ifndef DensityGridGenerator_hpp
#define DensityGridGenerator_hpp

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityGrid.hpp"
#include "MultiTimer.hpp"

namespace dengridgen {  // dengridgen namespace

/**
 Generates density for LDA.

 @param rho the pointer to density.
 @param gtoValues the GTO values on grid points.
 @param densityMatrix the density matrix.
 @param timer the timer.
 */
auto generateDensityForLDA(double* rho, const CDenseMatrix& gtoValues, const CDenseMatrix& densityMatrix, CMultiTimer& timer) -> void;

/**
 Generates density for LDA.

 @param rho the pointer to density.
 @param gtoValues the GTO values on grid points.
 @param densityMatrixAlpha the alpha density matrix.
 @param densityMatrixBeta the beta density matrix.
 @param timer the timer.
 */
auto generateDensityForLDA(double*             rho,
                           const CDenseMatrix& gtoValues,
                           const CDenseMatrix& densityMatrixAlpha,
                           const CDenseMatrix& densityMatrixBeta,
                           CMultiTimer&        timer) -> void;

/**
 Generates density for GGA.

 @param rho the pointer to density.
 @param rhograd the pointer to density gradient.
 @param sigma the pointer to dot product of density gradient.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO X derivative values on grid points.
 @param gtoValuesY the GTO Y derivative values on grid points.
 @param gtoValuesZ the GTO Z derivative values on grid points.
 @param densityMatrix the density matrix.
 @param timer the timer.
 */
auto generateDensityForGGA(double*             rho,
                           double*             rhograd,
                           double*             sigma,
                           const CDenseMatrix& gtoValues,
                           const CDenseMatrix& gtoValuesX,
                           const CDenseMatrix& gtoValuesY,
                           const CDenseMatrix& gtoValuesZ,
                           const CDenseMatrix& densityMatrix,
                           CMultiTimer&        timer) -> void;

/**
 Generates density for GGA.

 @param rho the pointer to density.
 @param rhograd the pointer to density gradient.
 @param sigma the pointer to dot product of density gradient.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO X derivative values on grid points.
 @param gtoValuesY the GTO Y derivative values on grid points.
 @param gtoValuesZ the GTO Z derivative values on grid points.
 @param densityMatrixAlpha the alpha density matrix.
 @param densityMatrixBeta the beta density matrix.
 @param timer the timer.
 */
auto generateDensityForGGA(double*             rho,
                           double*             rhograd,
                           double*             sigma,
                           const CDenseMatrix& gtoValues,
                           const CDenseMatrix& gtoValuesX,
                           const CDenseMatrix& gtoValuesY,
                           const CDenseMatrix& gtoValuesZ,
                           const CDenseMatrix& densityMatrixAlpha,
                           const CDenseMatrix& densityMatrixBeta,
                           CMultiTimer&        timer) -> void;

/**
 Generates density for meta-GGA.

 @param rho the pointer to density.
 @param rhograd the pointer to density gradient.
 @param sigma the pointer to dot product of density gradient.
 @param lapl the pointer to laplacian.
 @param tau the pointer to tau.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO X derivative values on grid points.
 @param gtoValuesY the GTO Y derivative values on grid points.
 @param gtoValuesZ the GTO Z derivative values on grid points.
 @param densityMatrix the density matrix.
 @param timer the timer.
 */
auto generateDensityForMGGA(double*             rho,
                            double*             rhograd,
                            double*             sigma,
                            double*             lapl,
                            double*             tau,
                            const CDenseMatrix& gtoValues,
                            const CDenseMatrix& gtoValuesX,
                            const CDenseMatrix& gtoValuesY,
                            const CDenseMatrix& gtoValuesZ,
                            const CDenseMatrix& densityMatrix,
                            CMultiTimer&        timer) -> void;

/**
 Generates density for meta-GGA.

 @param rho the pointer to density.
 @param rhograd the pointer to density gradient.
 @param sigma the pointer to dot product of density gradient.
 @param lapl the pointer to laplacian.
 @param tau the pointer to tau.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO X derivative values on grid points.
 @param gtoValuesY the GTO Y derivative values on grid points.
 @param gtoValuesZ the GTO Z derivative values on grid points.
 @param densityMatrixAlpha the alpha density matrix.
 @param densityMatrixBeta the beta density matrix.
 @param timer the timer.
 */
auto generateDensityForMGGA(double*             rho,
                            double*             rhograd,
                            double*             sigma,
                            double*             lapl,
                            double*             tau,
                            const CDenseMatrix& gtoValues,
                            const CDenseMatrix& gtoValuesX,
                            const CDenseMatrix& gtoValuesY,
                            const CDenseMatrix& gtoValuesZ,
                            const CDenseMatrix& densityMatrixAlpha,
                            const CDenseMatrix& densityMatrixBeta,
                            CMultiTimer&        timer) -> void;

}  // namespace dengridgen

#endif /* DensityGridGenerator_hpp */
