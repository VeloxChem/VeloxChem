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
         const int32_t              nElements, 
         const int32_t              nBlocks)
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
    
void
transform(      CMemBlock2D<double>& spherData,
          const CMemBlock2D<double>& cartData,
          const CSphericalMomentum&  spherMomentum,
          const int32_t              nElements,
          const int32_t              nBlocks)
{
    auto ncomp = spherMomentum.getNumberOfComponents();
    
    for (int32_t i = 0; i < ncomp; i++)
    {
        // set up Cartesian to spherical transformation data
        
        auto nfact = spherMomentum.getNumberOfFactors(i);
        
        auto tidx = spherMomentum.getIndexes(i);
        
        auto tfact = spherMomentum.getFactors(i);
        
        for (int32_t j = 0; j < nBlocks; j++)
        {
            // set up spherical data vector
            
            auto sphervec = spherData.data(i * nBlocks + j);
            
            // first term: assignment
            
            auto cartvec = cartData.data(tidx[0] * nBlocks + j);
            
            auto cfact = tfact[0];
            
            #pragma omp simd aligned(sphervec, cartvec: VLX_ALIGN)
            for (int32_t k = 0; k < nElements; k++)
            {
                sphervec[k] = cfact * cartvec[k];
            }
            
            // remaining terms: addition
            
            for (int32_t k = 1; k < nfact; k++)
            {
                cartvec = cartData.data(tidx[k] * nBlocks + j);
                
                cfact = tfact[k];
                
                #pragma omp simd aligned(sphervec, cartvec: VLX_ALIGN)
                for (int32_t l = 0; l < nElements; l++)
                {
                    sphervec[l] += cfact * cartvec[l];
                }
            }
        }
    }
}
    
} // genfunc namespace
