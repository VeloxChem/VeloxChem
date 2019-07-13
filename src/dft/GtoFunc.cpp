//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "GtoFunc.hpp"

#include "GtoFuncForLDA.hpp"

namespace gtorec {  // gtorec namespace
    
    void
    computeGtoValuesOnGrid(      CMemBlock2D<double>& spherGtoGridBuffer,
                                 CMemBlock2D<double>& cartGtoGridBuffer,
                           const double*              gridCoordinatesX,
                           const double*              gridCoordinatesY,
                           const double*              gridCoordinatesZ,
                           const int32_t              gridOffset,
                           const CGtoBlock&           gtoBlock,
                           const int32_t              iContrGto,
                           const xcfun                xcFunctional)
    {
        // set up angular momentum
        
        auto mang = gtoBlock.getAngularMomentum();
        
        // local density approximation
        
        if (xcFunctional == xcfun::lda)
        {
            switch (mang)
            {
                    // s-type GTOs on grid
                    
                case 0:
                    ldarec::compGtoValuesForS(spherGtoGridBuffer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, gtoBlock, iContrGto);
                    
                    break;
                    
                    // p-type GTOs on grid
                    
                case 1:
                    ldarec::compGtoValuesForP(spherGtoGridBuffer, cartGtoGridBuffer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, gtoBlock, iContrGto);
                    
                    break;
                    
                    // d-type GTOs on grid
                    
                case 2:
                    ldarec::compGtoValuesForD(spherGtoGridBuffer, cartGtoGridBuffer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, gtoBlock, iContrGto);
                    
                    break;
                    
                    // FIX ME: implement l > 2 cases
                    
                default:
                    break;
            }
            
            return;
        }
        
        // generalized gradient approximation
        
        if (xcFunctional == xcfun::gga)
        {
            // FIX ME: GGA case
            
            return;
        }
        
        // meta generalized gradient approximation
        
        if (xcFunctional == xcfun::mgga)
        {
            // FIX ME: MGGA case
            
            return;
        }
    }
    
}  // namespace gtorec
