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

#ifndef GtoFuncForMGGA_hpp
#define GtoFuncForMGGA_hpp

#include <cstdint>

#include "GtoContainer.hpp"
#include "MemBlock2D.hpp"

namespace mggarec {  // mggarec namespace
    
    /**
     Computes S-type GTOs values on given grid.
     
     @param spherGtoGridBuffer the buffer for storing contracted spherical GTOs values on the grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridOffset the batch offset in vector grid points.
     @param gtoBlock the GTOs block.
     @param iContrGto the index of contracted GTO is GTOs block.
     */
    void compGtoValuesForS(      CMemBlock2D<double>& spherGtoGridBuffer,
                           const double*              gridCoordinatesX,
                           const double*              gridCoordinatesY,
                           const double*              gridCoordinatesZ,
                           const int32_t              gridOffset,
                           const CGtoBlock&           gtoBlock,
                           const int32_t              iContrGto);
    
    
    /**
     Computes P-type GTOs values on given grid.
     
     @param spherGtoGridBuffer the buffer for storing contracted spherical GTOs values on the grid.
     @param cartGtoGridBuffer the buffer for storing primitive Cartesian GTOs values on the grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridOffset the batch offset in vector grid points.
     @param gtoBlock the GTOs block.
     @param iContrGto the index of contracted GTO is GTOs block.
     */
    void compGtoValuesForP(      CMemBlock2D<double>& spherGtoGridBuffer,
                                 CMemBlock2D<double>& cartGtoGridBuffer,
                           const double*              gridCoordinatesX,
                           const double*              gridCoordinatesY,
                           const double*              gridCoordinatesZ,
                           const int32_t              gridOffset,
                           const CGtoBlock&           gtoBlock,
                           const int32_t              iContrGto);
    
    /**
     Computes D-type GTOs values on given grid.
     
     @param spherGtoGridBuffer the buffer for storing contracted spherical GTOs values on the grid.
     @param cartGtoGridBuffer the buffer for storing primitive Cartesian GTOs values on the grid.
     @param gridCoordinatesX the vector of Cartesian X coordinates of grid points.
     @param gridCoordinatesY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordinatesZ the vector of Cartesian Y coordinates of grid points.
     @param gridOffset the batch offset in vector grid points.
     @param gtoBlock the GTOs block.
     @param iContrGto the index of contracted GTO is GTOs block.
     */
    void compGtoValuesForD(      CMemBlock2D<double>& spherGtoGridBuffer,
                                 CMemBlock2D<double>& cartGtoGridBuffer,
                           const double*              gridCoordinatesX,
                           const double*              gridCoordinatesY,
                           const double*              gridCoordinatesZ,
                           const int32_t              gridOffset,
                           const CGtoBlock&           gtoBlock,
                           const int32_t              iContrGto);
    
}  // namespace mggarec

#endif /* GtoFuncForMGGA_hpp */
