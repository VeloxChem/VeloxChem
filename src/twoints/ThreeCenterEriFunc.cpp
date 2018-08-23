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
        if (iContrGto == 0) printf("-> computing VRR(0|00)\n");
        
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
            
            // distribute (S|g(r,r')|SS) integrals
            
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
    
    void
    compElectronicPotentialForSSP(      CMemBlock2D<double>&  primBuffer,
                                  const CVecThreeIndexes&     recPattern,
                                  const std::vector<int32_t>& recIndexes,
                                  const CMemBlock2D<double>&  wdDistances,
                                  const CGtoBlock&            braGtoBlock,
                                  const CGtoPairsBlock&       ketGtoPairsBlock,
                                  const int32_t               iContrGto)
    {
        auto bord = genfunc::maxOrderOfPair(recPattern, 0, 1);
        
        // skip integrals if not included in recursion pattern
        
        if (bord < 0) return;
        
        if (iContrGto == 0) printf("-> computing VRR(0|01)\n");
        
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoPairsBlock.getNumberOfScreenedPrimPairs();
        
        // set up pointers to distances R(QD)
        
        auto qdx = ketGtoPairsBlock.getDistancesPBX();
        
        auto qdy = ketGtoPairsBlock.getDistancesPBY();
        
        auto qdz = ketGtoPairsBlock.getDistancesPBZ();
        
        // compute primitive integrals up to required order
        
        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer
            
            auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                 {0, 1, i});
            
            // skip integrals if this order is not required
            
            if (goff == -1) continue;
            
            // get position of integrals in primitves buffer
            
            auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 0, i});
            
            auto g2off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 0, i + 1});
            
            // loop over contracted GTO on bra side
            
            int32_t idx = 0;
            
            for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
            {
                // set up pointers to distances R(WD)
                
                auto wdx = wdDistances.data(3 * idx);
                
                auto wdy = wdDistances.data(3 * idx + 1);
                
                auto wdz = wdDistances.data(3 * idx + 2);
                
                // set up pointers to (S|g(r,r')|SS)^(m) integrals
                
                auto g0_0_0 = primBuffer.data(g1off + idx);
                
                // set up pointers to (S|g(r,r')|SS)^(m+1) integrals
                
                auto g1_0_0 = primBuffer.data(g2off + idx);
                
                // set up pointers to (S|g(r,r')|SP)^(m) integrals
                
                auto g_0_x = primBuffer.data(goff + 3 * idx);
                
                auto g_0_y = primBuffer.data(goff + 3 * idx + 1);
                
                auto g_0_z = primBuffer.data(goff + 3 * idx + 2);
                
                #pragma omp simd aligned(qdx, qdy, qdz, wdx, wdy, wdz, g_0_x,\
                                         g_0_y, g_0_z, g0_0_0, g1_0_0: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    double fact0 = g0_0_0[k];
                    
                    double fact1 = g1_0_0[k];
                    
                    g_0_x[k] = fact0 * qdx[k] + fact1 * wdx[k];
                    
                    g_0_y[k] = fact0 * qdy[k] + fact1 * wdy[k];
                    
                    g_0_z[k] = fact0 * qdz[k] + fact1 * wdz[k];
                }
                
                idx++;
            }
        }
    }
    
    void
    compElectronicPotentialForPSS(      CMemBlock2D<double>&  primBuffer,
                                  const CVecThreeIndexes&     recPattern,
                                  const std::vector<int32_t>& recIndexes,
                                  const CMemBlock2D<double>&  waDistances,
                                  const CGtoBlock&            braGtoBlock,
                                  const CGtoPairsBlock&       ketGtoPairsBlock,
                                  const int32_t               iContrGto)
    {
        auto bord = genfunc::maxOrderOfPair(recPattern, 1, 0);
        
        // skip integrals if not included in recursion pattern
        
        if (bord < 0) return;
        
        if (iContrGto == 0) printf("-> computing VRR(1|00)\n");
        
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoPairsBlock.getNumberOfScreenedPrimPairs();
        
        // compute primitive integrals up to required order
        
        for (int32_t i = 0; i <= bord; i++)
        {
            // get position of integrals in primitves buffer
            
            auto goff = genfunc::findTripleIndex(recIndexes, recPattern,
                                                 {1, 0, i});
            
            // skip integrals if this order is not required
            
            if (goff == -1) continue;
            
            // get position of integrals in primitves buffer
            
            auto g1off = genfunc::findTripleIndex(recIndexes, recPattern,
                                                  {0, 0, i + 1});
            
            // loop over contracted GTO on bra side
            
            int32_t idx = 0;
            
            for (int32_t j = spos[iContrGto]; j < epos[iContrGto]; j++)
            {
                // set up pointers to distances R(WA)
                
                auto wax = waDistances.data(3 * idx);
                
                auto way = waDistances.data(3 * idx + 1);
                
                auto waz = waDistances.data(3 * idx + 2);
                
                // set up pointers to (S|g(r,r')|SS)^(m+1) integrals
                
                auto g0_0_0 = primBuffer.data(g1off + idx);
                
                // set up pointers to (P|g(r,r')|SS)^(m) integrals
                
                auto g_0_x = primBuffer.data(goff + 3 * idx);
                
                auto g_0_y = primBuffer.data(goff + 3 * idx + 1);
                
                auto g_0_z = primBuffer.data(goff + 3 * idx + 2);
                
                #pragma omp simd aligned(wax, way, waz, g_0_x, g_0_y, g_0_z,\
                                         g0_0_0: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    double fact = g0_0_0[k];
                    
                    g_0_x[k] = fact * wax[k];
                    
                    g_0_y[k] = fact * way[k];
                    
                    g_0_z[k] = fact * waz[k];
                }
                
                idx++;
            }
        }
    }
    
}  // t3erifunc namespace
