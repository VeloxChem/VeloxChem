//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ThreeCenterEriFunc.hpp"

#include <cmath>

#include "MathConst.hpp"
#include "GenFunc.hpp"

namespace t3erifunc { // t3erifunc namespace
    
    void
    compElectronicPotentialForSSS(      CMemBlock2D<double>&  primBuffer,
                                  const CVecThreeIndexes&     recPattern,
                                  const std::vector<int32_t>& recIndexes,
                                  const CBoysFunction&        bfTable,
                                        CMemBlock<double>&    bfArguments,
                                        CMemBlock2D<double>&  bfValues,
                                  const int32_t               bfOrder,
                                  const CMemBlock2D<double>&  osFactors,
                                  const CMemBlock2D<double>&  aqDistances,
                                  const CGtoBlock&            braGtoBlock,
                                  const CGtoPairsBlock&       ketGtoPairsBlock,
                                  const int32_t               iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto bnorm = braGtoBlock.getNormFactors();
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto kfss = ketGtoPairsBlock.getOverlaps();
        
        auto nprim = ketGtoPairsBlock.getNumberOfScreenedPrimPairs();
        
        // set up pointers to R(AQ) distances
        
        auto aqx = aqDistances.data(0);
        
        auto aqy = aqDistances.data(1);
        
        auto aqz = aqDistances.data(2);
        
        // set up up pi prefactor
        
        auto fpi = 2.0 * mathconst::getPiValue();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up Obara-Saika prefactors
            
            auto fz = osFactors.data(5 * idx + 1);
            
            // compute Boys function argument
            
            auto fargs = bfArguments.data();
            
            #pragma omp simd aligned(fargs, fz, aqx, aqy, aqz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fargs[j] = fz[j] * (aqx[j] * aqx[j] + aqy[j] * aqy[j] +
                                    
                                    aqz[j] * aqz[j]);
            }
            
            // evaluate Boys function values
            
            bfTable.compute(bfValues, bfArguments, bfOrder);
            
            // set up pointers to Obara-Saika factors
            
            auto fia = osFactors.data(5 * idx + 2);
            
            auto fb = bnorm[i];
            
            // compute overlap scaling factor
            
            #pragma omp simd aligned(fz, fia, kfss, fargs: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fargs[j] = fb * kfss[j] * fpi * std::sqrt(fz[j]) * std::pow(fia[j], 1.50);
            }
            
            // distribute (s|g(r,r')|ss) integrals
            
            for (int32_t j = 0; j <= bfOrder; j++)
            {
                auto pidx = genfunc::findTripleIndex(recIndexes, recPattern,
                                                     {0, 0, j});
                
                if (pidx != -1)
                {
                    auto g_0_00 = primBuffer.data(pidx + idx);
                    
                    auto bvals = bfValues.data(j);
                    
                    #pragma omp simd aligned(g_0_00, bvals, fargs: VLX_ALIGN)
                    for (int32_t k = 0; k < nprim; k++)
                    {
                        g_0_00[k] = bvals[k] * fargs[k];
                    }
                }
            }
            
            idx++;
        }
    }
    
}  // t3erifunc namespace
