//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricFieldRecFunc.hpp"

#include <cmath>

#include "MathConst.hpp"
#include "GenFunc.hpp"

namespace efieldrecfunc { // npotrecfunc namespace
    
    void
    compElectricFieldForSS(      CMemBlock2D<double>&  primBuffer,
                           const CVecFourIndexes&      recPattern,
                           const std::vector<int32_t>& recIndexes,
                           const CBoysFunction&        bfTable,
                                 CMemBlock<double>&    bfArguments,
                                 CMemBlock2D<double>&  bfValues,
                           const int32_t               bfOrder,
                           const CMemBlock2D<double>&  osFactors,
                           const CMemBlock2D<double>&  abDistances,
                           const CMemBlock2D<double>&  pcDistances,
                           const CGtoBlock&            braGtoBlock,
                           const CGtoBlock&            ketGtoBlock,
                           const int32_t               iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto bnorm = braGtoBlock.getNormFactors();
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        auto bdim = epos[iContrGto] - spos[iContrGto];
        
        // set up pointers to primitives data on ket side
        
        auto knorm = ketGtoBlock.getNormFactors();
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // set up pointers to R(AB) distances
        
        auto abx = abDistances.data(0);
        
        auto aby = abDistances.data(1);
        
        auto abz = abDistances.data(2);
        
        // compute all (s|A(0)|s)^(m) terms
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up Obara-Saika prefactors
            
            auto fg = osFactors.data(3 * idx + 2);
            
            // set up pointers to ditances R(PC)
            
            auto pcx = pcDistances.data(3 * idx);
            
            auto pcy = pcDistances.data(3 * idx + 1);
            
            auto pcz = pcDistances.data(3 * idx + 2);
            
            // compute Boys function argument
            
            auto fargs = bfArguments.data();
            
            #pragma omp simd aligned(fargs, fg, pcx, pcy, pcz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fargs[j] = fg[j] * (pcx[j] * pcx[j] + pcy[j] * pcy[j] +
                                    
                                    pcz[j] * pcz[j]);
            }
            
            // evaluate Boys function values
            
            bfTable.compute(bfValues, bfArguments, bfOrder);
            
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(3 * idx);
            
            auto fz = osFactors.data(3 * idx + 1);
            
            auto fb = bnorm[i];
            
            // fetch up pi values
            
            auto fpi = 2.0 * mathconst::getPiValue();
            
            // compute overlap scaling factor
            
            #pragma omp simd aligned(fx, fz, knorm, abx, aby, abz,\
                                     fargs: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fargs[j] = fb * knorm[j] * fpi * fx[j]
                
                         * std::exp(-fz[j] * (abx[j] * abx[j] + aby[j] * aby[j] +
                                     
                                     abz[j] * abz[j]));
            }
            
            // distribute (s|A(0)|s) integrals
            
            for (int32_t j = 0; j <= bfOrder; j++)
            {
                auto pidx = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {0, 0, j, 1});
                
                if (pidx == -1) continue;
                
                auto t_0_0 = primBuffer.data(pidx + idx);
                
                auto bvals = bfValues.data(j);
                
                #pragma omp simd aligned(t_0_0, bvals, fargs: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    t_0_0[k] = bvals[k] * fargs[k];
                }
            }
            
            idx++;
        }
        
        // compute all (s|A(1)|s)^(m) terms
        
        auto maxord = genfunc::maxOrderOfTriple(recPattern, 0, 0, 0);
        
        idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to ditances R(PC)
            
            auto pcx = pcDistances.data(3 * idx);
            
            auto pcy = pcDistances.data(3 * idx + 1);
            
            auto pcz = pcDistances.data(3 * idx + 2);
            
            // loop over requested orders
            
            for (int32_t j = 0; j <= maxord; j++)
            {
                auto pidx = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {0, 0, j, 0});
                
                if (pidx == -1) continue;
                
                auto sidx = genfunc::findQuadrupleIndex(recIndexes, recPattern,
                                                        {0, 0, j + 1, 1});
                
                auto frep = primBuffer.data(sidx + idx);
                
                auto fex = primBuffer.data(pidx + idx);
                
                auto fey = primBuffer.data(pidx + bdim + idx);
                
                auto fez = primBuffer.data(pidx + 2 * bdim + idx);
                
                #pragma omp simd aligned(fex, fey, fez, frep: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    fex[k] = pcx[k] * frep[k];
                    
                    fey[k] = pcy[k] * frep[k];
                    
                    fez[k] = pcz[k] * frep[k];
                }
            }
           
            idx++;
        }
    }
    
} // efieldrecfunc namespace
