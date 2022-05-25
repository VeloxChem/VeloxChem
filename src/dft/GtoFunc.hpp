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

#ifndef GtoFunc_hpp
#define GtoFunc_hpp

#include <cstdint>

#include "DenseMatrix.hpp"
#include "GtoBlock.hpp"
#include "GtoContainer.hpp"
#include "MemBlock2D.hpp"
#include "MolecularGrid.hpp"
#include "XCFuncType.hpp"

namespace gtorec {  // gtorec namespace

    /**
     Computes GTOs values for batch of grid points.

     @param gtoValues the GTOs values buffer.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void computeGtosValuesForLDA(CMemBlock2D<double>& gtoValues,
                                 const CGtoContainer* gtoContainer,
                                 const double*        gridCoordinatesX,
                                 const double*        gridCoordinatesY,
                                 const double*        gridCoordinatesZ,
                                 const int32_t        gridOffset,
                                 const int32_t        gridBlockPosition,
                                 const int32_t        nGridPoints);

    /**
     Computes GTOs values for batch of grid points.

     @param gtoValues the GTOs values buffer.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void computeGtosValuesForGGA(CMemBlock2D<double>& gtoValues,
                                 CMemBlock2D<double>& gtoValuesX,
                                 CMemBlock2D<double>& gtoValuesY,
                                 CMemBlock2D<double>& gtoValuesZ,
                                 const CGtoContainer* gtoContainer,
                                 const double*        gridCoordinatesX,
                                 const double*        gridCoordinatesY,
                                 const double*        gridCoordinatesZ,
                                 const int32_t        gridOffset,
                                 const int32_t        gridBlockPosition,
                                 const int32_t        nGridPoints);

    /**
     Computes GTOs values for batch of grid points.

     @param gtoValues the GTOs values buffer.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void computeGtosValuesForMGGA(CMemBlock2D<double>& gtoValues,
                                  CMemBlock2D<double>& gtoValuesX,
                                  CMemBlock2D<double>& gtoValuesY,
                                  CMemBlock2D<double>& gtoValuesZ,
                                  const CGtoContainer* gtoContainer,
                                  const double*        gridCoordinatesX,
                                  const double*        gridCoordinatesY,
                                  const double*        gridCoordinatesZ,
                                  const int32_t        gridOffset,
                                  const int32_t        gridBlockPosition,
                                  const int32_t        nGridPoints);

    /**
     Computes GTOs values for batch of grid points.

     @param gtoValues the GTOs values buffer.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gtoBlock the GTOs block.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void computeGtosValuesForGGA(CMemBlock2D<double>& gtoValues,
                                 CMemBlock2D<double>& gtoValuesX,
                                 CMemBlock2D<double>& gtoValuesY,
                                 CMemBlock2D<double>& gtoValuesZ,
                                 const CGtoBlock&     gtoBlock,
                                 const double*        gridCoordinatesX,
                                 const double*        gridCoordinatesY,
                                 const double*        gridCoordinatesZ,
                                 const int32_t        gridOffset,
                                 const int32_t        gridBlockPosition,
                                 const int32_t        nGridPoints);

    /**
     Computes GTOs values for batch of grid points.

     @param gtoValues the GTOs values buffer.
     @param gtoValuesX the GTOs gradient along X axis values buffer.
     @param gtoValuesY the GTOs gradient along Y axis values buffer.
     @param gtoValuesZ the GTOs gradient along Z axis values buffer.
     @param gtoBlock the GTOs block.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void computeGtosValuesForMGGA(CMemBlock2D<double>& gtoValues,
                                  CMemBlock2D<double>& gtoValuesX,
                                  CMemBlock2D<double>& gtoValuesY,
                                  CMemBlock2D<double>& gtoValuesZ,
                                  const CGtoBlock&     gtoBlock,
                                  const double*        gridCoordinatesX,
                                  const double*        gridCoordinatesY,
                                  const double*        gridCoordinatesZ,
                                  const int32_t        gridOffset,
                                  const int32_t        gridBlockPosition,
                                  const int32_t        nGridPoints);

    /**
     Computes GTOs values matrix for given block of grid points.

     @param gtoMatrix the GTOs values matrix.
     @param gtoContainer the container of GTOs blocks.
     @param molecularGrid the molecular grid.
     @param gridOffset the offset in molecular grid.
     @param nGridPoints the number of grid points.
     */
    void computeGtosMatrixForLDA(CDenseMatrix&         gtoMatrix,
                                 const CGtoContainer*  gtoContainer,
                                 const CMolecularGrid& molecularGrid,
                                 const int32_t         gridOffset,
                                 const int32_t         nGridPoints);

    /**
     Computes GTOs values matrix for given block of grid points.

     @param gtoMatrix the GTOs values matrix.
     @param gtoMatrixX the GTOs gradient along X coordinate values matrix.
     @param gtoMatrixY the GTOs gradient along Y coordinate values matrix.
     @param gtoMatrixZ the GTOs gradient along Z coordinate values matrix.
     @param gtoContainer the container of GTOs blocks.
     @param molecularGrid the molecular grid.
     @param gridOffset the offset in molecular grid.
     @param nGridPoints the number of grid points.
     */
    void computeGtosMatrixForGGA(CDenseMatrix&         gtoMatrix,
                                 CDenseMatrix&         gtoMatrixX,
                                 CDenseMatrix&         gtoMatrixY,
                                 CDenseMatrix&         gtoMatrixZ,
                                 const CGtoContainer*  gtoContainer,
                                 const CMolecularGrid& molecularGrid,
                                 const int32_t         gridOffset,
                                 const int32_t         nGridPoints);

    /**
     Computes GTOs values for batch of grid points.

     @param gtoMatrix the pointer to GTOs values matrix.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void computeGtosValuesForLDA(double*              gtoMatrix,
                                 const CGtoContainer* gtoContainer,
                                 const double*        gridCoordinatesX,
                                 const double*        gridCoordinatesY,
                                 const double*        gridCoordinatesZ,
                                 const int32_t        gridOffset,
                                 const int32_t        gridBlockPosition,
                                 const int32_t        nGridPoints);

    /**
     Computes GTOs values for batch of grid points.

     @param gtoMatrix the pointer to GTOs values matrix.
     @param gtoMatrixX the pointer to GTOs gradient along X coordinate values matrix.
     @param gtoMatrixY the pointer to GTOs gradient along Y coordinate values matrix.
     @param gtoMatrixZ the pointer to GTOs gradient along Z coordinate values matrix.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void computeGtosValuesForGGA(double*              gtoMatrix,
                                 double*              gtoMatrixX,
                                 double*              gtoMatrixY,
                                 double*              gtoMatrixZ,
                                 const CGtoContainer* gtoContainer,
                                 const double*        gridCoordinatesX,
                                 const double*        gridCoordinatesY,
                                 const double*        gridCoordinatesZ,
                                 const int32_t        gridOffset,
                                 const int32_t        gridBlockPosition,
                                 const int32_t        nGridPoints);

    /**
     Computes GTOs values for batch of grid points.

     @param gtoMatrix the pointer to GTOs values matrix.
     @param gtoMatrixX the pointer to GTOs gradient along X coordinate values matrix.
     @param gtoMatrixY the pointer to GTOs gradient along Y coordinate values matrix.
     @param gtoMatrixZ the pointer to GTOs gradient along Z coordinate values matrix.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridOffset the offset of grid points batch in molecular grid.
     @param gridBlockPosition the position of grid block in GTOs values grid.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void computeGtosValuesForMGGA(double*              gtoMatrix,
                                  double*              gtoMatrixX,
                                  double*              gtoMatrixY,
                                  double*              gtoMatrixZ,
                                  const CGtoContainer* gtoContainer,
                                  const double*        gridCoordinatesX,
                                  const double*        gridCoordinatesY,
                                  const double*        gridCoordinatesZ,
                                  const int32_t        gridOffset,
                                  const int32_t        gridBlockPosition,
                                  const int32_t        nGridPoints);

    /**
     Computes contracted GTOs values at grid points for specific type of functional.

     @param spherGtoGridBuffer the buffer for storing contracted spherical GTOs values on the grid.
     @param cartGtoGridBuffer the buffer for storing primitive Cartesian GTOs values on the grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid
     points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid
     points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid
     points.
     @param gridOffset the batch offset in vector grid points.
     @param gtoBlock the GTOs block.
     @param iContrGto the index of contracted GTO is GTOs block.
     @param xcFunctional the exchange-correlations functional type.
     */
    void computeGtoValuesOnGrid(CMemBlock2D<double>& spherGtoGridBuffer,
                                CMemBlock2D<double>& cartGtoGridBuffer,
                                const double*        gridCoordinatesX,
                                const double*        gridCoordinatesY,
                                const double*        gridCoordinatesZ,
                                const int32_t        gridOffset,
                                const CGtoBlock&     gtoBlock,
                                const int32_t        iContrGto,
                                const xcfun          xcFunctional);

}  // namespace gtorec

#endif /* GtoFunc_hpp */
