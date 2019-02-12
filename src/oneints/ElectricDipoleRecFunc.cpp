//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricDipoleRecFunc.hpp"

#include "MathConst.hpp"
#include "GenFunc.hpp"

namespace ediprecfunc { // ediprecfunc namespace
    
    void
    compElectricDipoleForSS(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
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
        
        // fetch up pi values
        
        auto fpi = mathconst::getPiValue();
        
        // get position of integrals in primitves buffer
        
        auto soff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 1});
        
        auto doff = genfunc::findTripleIndex(recIndexes, recPattern, {0, 0, 0});
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(2 * idx);
            
            auto fz = osFactors.data(2 * idx + 1);
            
            auto fb = bnorm[i];
            
            // set up pointers to ditances R(PC)
            
            auto pcx = pcDistances.data(3 * idx);
            
            auto pcy = pcDistances.data(3 * idx + 1);
            
            auto pcz = pcDistances.data(3 * idx + 2);
            
            // set up primitives buffer data
            
            auto fovl = primBuffer.data(soff + idx);
            
            auto fdipx = primBuffer.data(doff + idx);
            
            auto fdipy = primBuffer.data(doff + bdim + idx);
            
            auto fdipz = primBuffer.data(doff + 2 * bdim + idx);
            
            #pragma omp simd aligned(fovl, fx, fz, knorm, abx, aby,\
                                     abz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fovl[j] = fb * knorm[j] * std::pow(fpi * fx[j], 1.5)
                
                        * std::exp(-fz[j] * (abx[j] * abx[j] + aby[j] * aby[j] +
                                     
                                             abz[j] * abz[j]));
                
                fdipx[j] = pcx[j] * fovl[j];
                
                fdipy[j] = pcy[j] * fovl[j];
                
                fdipz[j] = pcz[j] * fovl[j];
            }
            
            idx++;
        }
    }
    
} // ediprecfunc namespace
