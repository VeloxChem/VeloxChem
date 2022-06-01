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

#ifndef DiagEriRecForSPSP_hpp
#define DiagEriRecForSPSP_hpp

#include <cstdint>

#include "Buffer.hpp"
#include "BinnedGtoPairBlock.hpp"
#include "DiagEriRecFacts.hpp"
#include "BoysFunc.hpp"

namespace derirec { // derirec  namespace

/**
Computes diagonal (SP|SP) integrals batch and distributes them into given
integrals buffer.

@param intsBuffer The integrals buffer.
@param gtoPairBlock The pointer to GTOs pairs block.
@param bPosition The start position of contracted GTOs batch.
@param ePosition The endposition of contracted GTOs batch.
*/
template <typename T>
auto
compHostSPSP(      T*                                 intsBuffer,
             const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock,
             const int32_t                            bPosition,
             const int32_t                            ePosition) -> void
{
    // set up dimensions
    
    const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();
    
    const auto ncpairs = ePosition - bPosition;
    
    // allocate recursion buffers
    
    BufferHostMY<T, 3> rpq(ncpairs);
    
    BufferHostMY<T, 3> osfacts(ncpairs);
    
    BufferHostMY<T, 3> rpb(ncpairs);
    
    BufferHostMY<T, 3> rqd(ncpairs);
    
    BufferHostMY<T, 3> rw(ncpairs);
    
    BufferHostMY<T, 3> rwp(ncpairs);
    
    BufferHostMY<T, 3> rwq(ncpairs);
    
    // allocate Boys function data
    
    BufferHostX<T> bargs(ncpairs);
    
    BufferHostXY<T> bvals(3, ncpairs);
    
    CBoysFunc<T, 2> bftable;
    
    // allocate contracted integral buffers
    
    auto t_spsp = BufferHostMY<T, 3>::Zero(ncpairs);

    // allocate primitive integral buffers
    
    BufferHostMY<T, 3> t0_sssp(ncpairs);
    
    BufferHostMY<T, 3> t1_sssp(ncpairs);
    
    BufferHostXY<T> t_ssss(3, ncpairs);
    
    // set up pointers to Obara-Saika factors
        
    auto frho = osfacts.data(0);
        
    auto fnorm = osfacts.data(1);
    
    auto fzeta = osfacts.data(2);
        
    // loop over primitive integrals
    
    for (int32_t i = 0; i < nppairs; i++)
    {
        derirec::compHostDistancesPT(rpb, gtoPairBlock, bPosition, ePosition, i);
        
        for (int j = i; j < nppairs; j++)
        {
            // compute recursion data
            
            derirec::compHostDistancesPQ(rpq, gtoPairBlock,
                                         bPosition, ePosition, i, j);
            
            derirec::compHostFactorRho(frho, gtoPairBlock,
                                       bPosition, ePosition, i, j);
            
            derirec::compHostFactorNorm(fnorm, gtoPairBlock,
                                        bPosition, ePosition, i, j);
            
            derirec::compHostFactorZeta(fzeta, gtoPairBlock,
                                        bPosition, ePosition, i, j);
            
            derirec::compHostDistancesPT(rqd, gtoPairBlock, bPosition, ePosition, j);
            
            derirec::compHostCoordinatesW(rw, gtoPairBlock, bPosition, ePosition, i, j);
            
            if (i == j)
            {
                rwp.setZero();
                
                rwq.setZero();
            }
            else
            {
                derirec::compHostDistancesWT(rwp, rw, gtoPairBlock, bPosition, ePosition, i);
                
                derirec::compHostDistancesWT(rwq, rw, gtoPairBlock, bPosition, ePosition, j);
            }
            
            // compute Boys function values
            
            derirec::compHostBoysArguments(bargs, rpq, frho, ncpairs);
            
            bftable.compute(bvals, bargs);
            
            // compute vertical recursion data
            
            // accumulate integrals
        }
    }
    
    derirec::selectHostMaxValues<T, 3>(intsBuffer, t_spsp, bPosition, ncpairs); 
}

} // derirec  namespace

#endif /* DiagEriRecForSPSP_hpp */
