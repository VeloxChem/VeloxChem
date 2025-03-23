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

#ifndef SerialDensityGridGenerator_hpp
#define SerialDensityGridGenerator_hpp

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityGrid.hpp"

namespace sdengridgen {  // sdengridgen namespace

/**
 Generates density for LDA.

 @param rho the pointer to density.
 @param gtoValues the GTO values on grid points.
 @param densityMatrix the density matrix.
 */
auto serialGenerateDensityForLDA(double* rho, const CDenseMatrix& gtoValues, const CDenseMatrix& densityMatrix) -> void;

/**
 Generates density for LDA.

 @param rho the pointer to density.
 @param gtoValues the GTO values on grid points.
 @param densityMatrixAlpha the alpha density matrix.
 @param densityMatrixBeta the beta density matrix.
 */
auto serialGenerateDensityForLDA(double*             rho,
                                 const CDenseMatrix& gtoValues,
                                 const CDenseMatrix& densityMatrixAlpha,
                                 const CDenseMatrix& densityMatrixBeta) -> void;

/**
 Generates density grid (multiple density matrices) for LDA.

 @param gtoValues the GTO values on grid points.
 @param densityMatrix the AO density matrices.
 @param xcFunType the type of exchange-correlation functional.
 @return the density grid.
 */
auto serialGenerateDensityGridForLDA(const CDenseMatrix&     gtoValues,
                                     const CAODensityMatrix& densityMatrix,
                                     const xcfun             xcFunType) -> CDensityGrid;

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
 */
auto serialGenerateDensityForGGA(double*             rho,
                                 double*             rhograd,
                                 double*             sigma,
                                 const CDenseMatrix& gtoValues,
                                 const CDenseMatrix& gtoValuesX,
                                 const CDenseMatrix& gtoValuesY,
                                 const CDenseMatrix& gtoValuesZ,
                                 const CDenseMatrix& densityMatrix) -> void;

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
 */
auto serialGenerateDensityForGGA(double*             rho,
                                 double*             rhograd,
                                 double*             sigma,
                                 const CDenseMatrix& gtoValues,
                                 const CDenseMatrix& gtoValuesX,
                                 const CDenseMatrix& gtoValuesY,
                                 const CDenseMatrix& gtoValuesZ,
                                 const CDenseMatrix& densityMatrixAlpha,
                                 const CDenseMatrix& densityMatrixBeta) -> void;

/**
 Generates density grid (multiple density matrices) for GGA.

 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO gradient X values on grid points.
 @param gtoValuesY the GTO gradient Y values on grid points.
 @param gtoValuesZ the GTO gradient Z values on grid points.
 @param densityMatrix the AO density matrices.
 @param xcFunType the type of exchange-correlation functional.
 @return the density grid.
 */
auto serialGenerateDensityGridForGGA(const CDenseMatrix&     gtoValues,
                                     const CDenseMatrix&     gtoValuesX,
                                     const CDenseMatrix&     gtoValuesY,
                                     const CDenseMatrix&     gtoValuesZ,
                                     const CAODensityMatrix& densityMatrix,
                                     const xcfun             xcFunType) -> CDensityGrid;

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
 */
auto serialGenerateDensityForMGGA(double*             rho,
                                  double*             rhograd,
                                  double*             sigma,
                                  double*             lapl,
                                  double*             tau,
                                  const CDenseMatrix& gtoValues,
                                  const CDenseMatrix& gtoValuesX,
                                  const CDenseMatrix& gtoValuesY,
                                  const CDenseMatrix& gtoValuesZ,
                                  const CDenseMatrix& densityMatrix) -> void;

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
 */
auto serialGenerateDensityForMGGA(double*             rho,
                                  double*             rhograd,
                                  double*             sigma,
                                  double*             lapl,
                                  double*             tau,
                                  const CDenseMatrix& gtoValues,
                                  const CDenseMatrix& gtoValuesX,
                                  const CDenseMatrix& gtoValuesY,
                                  const CDenseMatrix& gtoValuesZ,
                                  const CDenseMatrix& densityMatrixAlpha,
                                  const CDenseMatrix& densityMatrixBeta) -> void;

/**
 Generates density grid (multiple density matrices) for MGGA.

 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO gradient X values on grid points.
 @param gtoValuesY the GTO gradient Y values on grid points.
 @param gtoValuesZ the GTO gradient Z values on grid points.
 @param densityMatrix the AO density matrices.
 @param xcFunType the type of exchange-correlation functional.
 @return the density grid.
 */
auto serialGenerateDensityGridForMGGA(const CDenseMatrix&     gtoValues,
                                      const CDenseMatrix&     gtoValuesX,
                                      const CDenseMatrix&     gtoValuesY,
                                      const CDenseMatrix&     gtoValuesZ,
                                      const CAODensityMatrix& densityMatrix,
                                      const xcfun             xcFunType) -> CDensityGrid;

}  // namespace sdengridgen

#endif /* SerialDensityGridGenerator_hpp */
