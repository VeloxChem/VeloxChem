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

#ifndef DiagEriRecForSSSS_hpp
#define DiagEriRecForSSSS_hpp

#include <cstdint>

#include "Buffer.hpp"
#include "BinnedGtoPairBlock.hpp"
#include "DiagEriRecFacts.hpp"
#include "BoysFunc.hpp"

namespace derirec { // derirec  namespace

/**
Computes diagonal (SS|SS) integrals batch and distributes them into given
integrals buffer.

@param intsBuffer The integrals buffer.
@param gtoPairBlock The pointer to GTOs pairs block.
@param bPosition The start position of contracted GTOs batch.
@param ePosition The endposition of contracted GTOs batch.
*/
template <typename T>
auto
compHostSSSS(      T*                                 intsBuffer,
             const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock,
             const int32_t                            bPosition,
             const int32_t                            ePosition) -> void
{
    // set up dimensions
    
    const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();
    
    const auto ncpairs = ePosition - bPosition;
    
    // allocate recursion buffers
    
    BufferHostMY<T, 3> rpq(ncpairs);
    
    BufferHostMY<T, 2> osfacts(ncpairs);
    
    // allocate Boys function data
    
    BufferHostX<T> bargs(ncpairs);
    
    BufferHostXY<T> bvals(1, ncpairs);
    
    CBoysFunc<T, 0> bftable;
    
    // set up scaling factor
        
    const auto fpi = static_cast<T>(4.0) / static_cast<T>(mathconst::getPiValue());
        
    // set up pointers to Obara-Saika factors
        
    auto frho = osfacts.data(0);
        
    auto fnorm = osfacts.data(1);
        
    // loop over primitive integrals
    
    for (int32_t i = 0; i < nppairs; i++)
    {
        for (int j = i; j < nppairs; j++)
        {
            // compute recursion data
            
            derirec::compHostDistancesPQ(rpq, gtoPairBlock,
                                         bPosition, ePosition, i, j);
            
            derirec::compHostFactorRho(frho, gtoPairBlock,
                                       bPosition, ePosition, i, j);
            
            derirec::compHostFactorNorm(fnorm, gtoPairBlock,
                                        bPosition, ePosition, i, j);
            
            // compute Boys function values
            
            derirec::compHostBoysArguments(bargs, rpq, frho, ncpairs);
            
            bftable.compute(bvals, bargs);
            
            // set up pointers to Boys function values
            
            auto bf0 = bvals.data();
            
            // compute (ss|ss) integrals
            
            if ((i + j) == 0)
            {
                #pragma omp simd aligned(frho, fnorm, bvals: VLX_ALIGN)
                for (int32_t k = 0; k < ncpairs; k++)
                {
                    intsBuffer[k] = std::sqrt(fpi * frho[k]) * bf0[k] * fnorm[k];
                }
            }
            else
            {
                #pragma omp simd aligned(frho, fnorm, bvals: VLX_ALIGN)
                for (int32_t k = 0; k < ncpairs; k++)
                {
                    intsBuffer[k] += std::sqrt(fpi * frho[k]) * bf0[k] * fnorm[k];
                }
            }
        }
    }
}

}  // derirec  namespace

#endif /* DiagEriRecForSSSS_hpp */
