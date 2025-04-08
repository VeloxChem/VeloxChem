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
