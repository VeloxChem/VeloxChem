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

#ifndef GridScreener_hpp
#define GridScreener_hpp

#include <cstdint>
#include <vector>

#include "DensityGrid.hpp"

namespace gridscreen {  // gridscreen namespace

/**
 Screens LDA density grid to get rid of invalid grid points.

 @param gridPointInds mapping between grid points before and after screening
 @param destDensityGrid the density grid to store the result (destination).
 @param srcDensityGrid the density grid to provide the grid points (source).
 @param densityThreshold the threshold for density grid screening.
 */
void screenDensityGridForLDA(std::vector<int32_t>& gridPointInds,
                             CDensityGrid&         destDensityGrid,
                             const CDensityGrid&   srcDensityGrid,
                             const double          densityThreshold);

/**
 Screens GGA density grid to get rid of invalid grid points.

 @param gridPointInds mapping between grid points before and after screening
 @param destDensityGrid the density grid to store the result (destination).
 @param srcDensityGrid the density grid to provide the grid points (source).
 @param densityThreshold the threshold for density grid screening.
 */
void screenDensityGridForGGA(std::vector<int32_t>& gridPointInds,
                             CDensityGrid&         destDensityGrid,
                             const CDensityGrid&   srcDensityGrid,
                             const double          densityThreshold);

/**
 Screens weights.

 @param screenedWeights pointer to the screened weights.
 @param gridBlockPosition the starting position of the grid box.
 @param weights pointer to the original weights of all grid points.
 @param gridPointInds mapping between grid points before and after screening
 @param nScreenedGridPoints the number of grid points after screening.
 */
void screenWeights(double*                     screenedWeights,
                   const int32_t               gridBlockPosition,
                   const double*               weights,
                   const std::vector<int32_t>& gridPointInds,
                   const int32_t               nScreenedGridPoints);

/**
 Screens GTO matrix for LDA.

 @param screenedGtoValues the matrix containing screened GTO values.
 @param originalGtoValues the matrix containing original GTO values.
 @param gridPointInds mapping between grid points before and after screening
 @param nScreenedGridPoints the number of grid points after screening.
 */
void screenGtoMatrixForLDA(CDenseMatrix&               screenedGtoValues,
                           const CDenseMatrix&         originalGtoValues,
                           const std::vector<int32_t>& gridPointInds,
                           const int32_t               nScreenedGridPoints);

/**
 Screens GTO matrix for GGA.

 @param screenedGtoValues the matrix containing screened GTO values.
 @param screenedGtoValuesX the matrix containing screened GTO X derivatives.
 @param screenedGtoValuesY the matrix containing screened GTO Y derivatives.
 @param screenedGtoValuesZ the matrix containing screened GTO Z derivatives.
 @param originalGtoValues the matrix containing original GTO values.
 @param originalGtoValuesX the matrix containing original GTO X derivatives.
 @param originalGtoValuesY the matrix containing original GTO Y derivatives.
 @param originalGtoValuesZ the matrix containing original GTO Z derivatives.
 @param gridPointInds mapping between grid points before and after screening
 @param nScreenedGridPoints the number of grid points after screening.
 */
void screenGtoMatrixForGGA(CDenseMatrix&               screenedGtoValues,
                           CDenseMatrix&               screenedGtoValuesX,
                           CDenseMatrix&               screenedGtoValuesY,
                           CDenseMatrix&               screenedGtoValuesZ,
                           const CDenseMatrix&         originalGtoValues,
                           const CDenseMatrix&         originalGtoValuesX,
                           const CDenseMatrix&         originalGtoValuesY,
                           const CDenseMatrix&         originalGtoValuesZ,
                           const std::vector<int32_t>& gridPointInds,
                           const int32_t               nScreenedGridPoints);

}  // namespace gridscreen

#endif /* GridScreener_hpp */
