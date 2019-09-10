//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef GtoFunc_hpp
#define GtoFunc_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MemBlock2D.hpp"
#include "XCFuncType.hpp"

namespace gtorec {  // gtorec namespace
    
   
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
    void computeGtoValuesOnGrid(      CMemBlock2D<double>& spherGtoGridBuffer,
                                      CMemBlock2D<double>& cartGtoGridBuffer,
                                const double*              gridCoordinatesX,
                                const double*              gridCoordinatesY,
                                const double*              gridCoordinatesZ,
                                const int32_t              gridOffset,
                                const CGtoBlock&           gtoBlock,
                                const int32_t              iContrGto,
                                const xcfun                xcFunctional);
    
}  // namespace gtorec

#endif /* GtoFunc_hpp */
