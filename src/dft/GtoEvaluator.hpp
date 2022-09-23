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

#ifndef GtoEvaluator_hpp
#define GtoEvaluator_hpp

#include <array>

#include "GtoContainer.hpp"
#include "MemBlock.hpp"
#include "MemBlock2D.hpp"
#include "XCFuncType.hpp"

namespace gtoeval {  // gtoeval namespace

/**
 Gets grid box dimension.

 @param gridBlockPosition the displacement of grid points in this box.
 @param nGridPoints the number of grid points in this box.
 @param xcoords the X coordinates of grid points.
 @param ycoords the Y coordinates of grid points.
 @param zcoords the Z coordinates of grid points.
 @return grid box dimension as (xmin, ymin, zmin, xmax, ymax, zmax).
 */
std::array<double, 6> getGridBoxDimension(const int32_t gridBlockPosition,
                                          const int32_t nGridPoints,
                                          const double* xcoords,
                                          const double* ycoords,
                                          const double* zcoords);

/**
 Prescreens GTOs for a grid box.

 @param skipCgtoIds the array to store whether a CGTO should be skipped.
 @param skipAOIds the array to store whether an AO should be skipped.
 @param gtoContainer the pointer to the GTO container.
 @param gtoDeriv the level of GTO derivative.
 @param gtoThreshold the screening threshold for GTO.
 @param boxDimension the dimension of the grid box.
 */
void preScreenGtos(CMemBlock<int32_t>&          skipCgtoIds,
                   CMemBlock<int32_t>&          skipAOIds,
                   const CGtoContainer*         gtoContainer,
                   const int32_t                gtoDeriv,
                   const double                 gtoThreshold,
                   const std::array<double, 6>& boxDimension);

/**
 Computes GTOs values for batch of grid points.

 @param gtoValues the GTOs values buffer.
 @param gtoContainer the poitner to GTO container.
 @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
 @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
 @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
 @param gridBlockPosition the starting position of grid box.
 @param gridOffset the offset of grid points in grid box.
 @param nGridPoints the number of grid points in grid points batch.
 @param skipCgtoids whether a CGTO should be skipped.
 */
void computeGtosValuesForLDA(CMemBlock2D<double>&      gtoValues,
                             const CGtoContainer*      gtoContainer,
                             const double*             gridCoordinatesX,
                             const double*             gridCoordinatesY,
                             const double*             gridCoordinatesZ,
                             const int32_t             gridBlockPosition,
                             const int32_t             gridOffset,
                             const int32_t             nGridPoints,
                             const CMemBlock<int32_t>& skipCgtoIds);

/**
 Computes GTOs values and derivatives for batch of grid points.

 @param gtoValues the GTOs values buffer.
 @param gtoValueX the GTOs X derivative values buffer.
 @param gtoValueY the GTOs Y derivative values buffer.
 @param gtoValueZ the GTOs Z derivative values buffer.
 @param gtoContainer the poitner to GTO container.
 @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
 @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
 @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
 @param gridBlockPosition the starting position of grid box.
 @param gridOffset the offset of grid points in grid box.
 @param nGridPoints the number of grid points in grid points batch.
 @param skipCgtoids whether a CGTO should be skipped.
 */
void computeGtosValuesForGGA(CMemBlock2D<double>&      gtoValues,
                             CMemBlock2D<double>&      gtoValuesX,
                             CMemBlock2D<double>&      gtoValuesY,
                             CMemBlock2D<double>&      gtoValuesZ,
                             const CGtoContainer*      gtoContainer,
                             const double*             gridCoordinatesX,
                             const double*             gridCoordinatesY,
                             const double*             gridCoordinatesZ,
                             const int32_t             gridBlockPosition,
                             const int32_t             gridOffset,
                             const int32_t             nGridPoints,
                             const CMemBlock<int32_t>& skipCgtoIds);

/**
 Computes GTOs values, 1st- and 2nd-order derivatives for batch of grid points.

 @param gtoValues the GTOs values buffer.
 @param gtoValueX the GTOs X derivative values buffer.
 @param gtoValueY the GTOs Y derivative values buffer.
 @param gtoValueZ the GTOs Z derivative values buffer.
 @param gtoValueXX the GTOs XX derivative values buffer.
 @param gtoValueXY the GTOs XY derivative values buffer.
 @param gtoValueXZ the GTOs XZ derivative values buffer.
 @param gtoValueYY the GTOs YY derivative values buffer.
 @param gtoValueYZ the GTOs YZ derivative values buffer.
 @param gtoValueZZ the GTOs ZZ derivative values buffer.
 @param gtoContainer the poitner to GTO container.
 @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
 @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
 @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
 @param gridBlockPosition the starting position of grid box.
 @param gridOffset the offset of grid points in grid box.
 @param nGridPoints the number of grid points in grid points batch.
 @param skipCgtoids whether a CGTO should be skipped.
 */
void computeGtosValuesForMetaGGA(CMemBlock2D<double>&      gtoValues,
                                 CMemBlock2D<double>&      gtoValuesX,
                                 CMemBlock2D<double>&      gtoValuesY,
                                 CMemBlock2D<double>&      gtoValuesZ,
                                 CMemBlock2D<double>&      gtoValuesXX,
                                 CMemBlock2D<double>&      gtoValuesXY,
                                 CMemBlock2D<double>&      gtoValuesXZ,
                                 CMemBlock2D<double>&      gtoValuesYY,
                                 CMemBlock2D<double>&      gtoValuesYZ,
                                 CMemBlock2D<double>&      gtoValuesZZ,
                                 const CGtoContainer*      gtoContainer,
                                 const double*             gridCoordinatesX,
                                 const double*             gridCoordinatesY,
                                 const double*             gridCoordinatesZ,
                                 const int32_t             gridBlockPosition,
                                 const int32_t             gridOffset,
                                 const int32_t             nGridPoints,
                                 const CMemBlock<int32_t>& skipCgtoIds);

}  // namespace gtoeval

#endif /* GtoEvaluator_hpp */
