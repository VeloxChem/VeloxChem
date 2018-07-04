//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GenFunc.hpp"

namespace genfunc { // genfunc namespace

void
contract(      CMemBlock2D<double>& contrData,
         const int32_t              contrIndex,
         const CMemBlock2D<double>& primData,
         const int32_t              primIndex,
         const int32_t*             startPositions,
         const int32_t*             endPositions,
         const int32_t              nBlocks,
         const int32_t              nElements)
{
    for (int32_t i = 0; i < nBlocks; i++)
    {
        // set up data vectors
        
        auto pdat = primData.data(primIndex + i);
        
        auto cdat = contrData.data(i);
        
        for (int32_t j = 0; j < nElements; j++)
        {
            double fsum = 0.0;
            
            // contract data vector components 
            
            for (int32_t k = startPositions[j]; k < endPositions[j]; k++)
            {
                fsum += pdat[k];
            }
            
            cdat[j] = fsum;
        }
    }
}
    
} // genfunc namespace
