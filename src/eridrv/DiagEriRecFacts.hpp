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

#ifndef DiagEriRecFacts_hpp
#define DiagEriRecFacts_hpp

#include <cstdint>

#include "Buffer.hpp"
#include "BinnedGtoPairBlock.hpp"

namespace derirec { // derirec  namespace

/**
Computes R(PQ) = P - Q distances for contracted GTOs batch.

@param rDistancesPQ The R(PQ) = P - Q distances.
@param gtoPairBlock The pointer to GTOs pairs block.
@param bPosition The start position of contracted GTOs batch.
@param ePosition The endposition of contracted GTOs batch.
@param braPrimGto The primitive GTO on bra side.
@param ketPrimGto The primitive GTO on ket side.
*/
template <typename T>
auto
compHostDistancesPQ(      BufferHostMY<T, 3>&                rDistancesPQ,
                    const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock,
                    const int32_t                            bPosition,
                    const int32_t                            ePosition,
                    const int32_t                            braPrimGto,
                    const int32_t                            ketPrimGto) -> void
{
    // set up dimentsions
    
    const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();
    
    const auto ncpairs = ePosition - bPosition;
    
    // set up pointers to P center coordinates
        
    auto rpx = gtoPairBlock->getCoordinatesPX();
        
    auto rpy = gtoPairBlock->getCoordinatesPY();
        
    auto rpz = gtoPairBlock->getCoordinatesPZ();
        
    // set up pointers to R(PQ) = P - Q distances
           
    auto rpqx = rDistancesPQ.data(0);
           
    auto rpqy = rDistancesPQ.data(1);
           
    auto rpqz = rDistancesPQ.data(2);

    // compute R(PQ) distances
    
    if (braPrimGto == ketPrimGto)
    {
        rDistancesPQ.setZero();
    }
    else
    {
        for (int32_t i = 0; i < ncpairs; i++)
        {
            const auto ioff = (bPosition + i) * nppairs;
            
            rpqx[i] = rpx[ioff + braPrimGto] - rpx[ioff + ketPrimGto];
                               
            rpqy[i] = rpy[ioff + braPrimGto] - rpy[ioff + ketPrimGto];
                               
            rpqz[i] = rpz[ioff + braPrimGto] - rpz[ioff + ketPrimGto];
        }
    }
}

/**
Computes rho = xi * eta / (xi + eta) factor for contracted GTOs batch.

@param osFactorRho The Obara-Saika factor rho = xi * eta / (xi + eta).
@param gtoPairBlock The pointer to GTOs pairs block.
@param bPosition The start position of contracted GTOs batch.
@param ePosition The endposition of contracted GTOs batch.
@param braPrimGto The primitive GTO on bra side.
@param ketPrimGto The primitive GTO on ket side.
*/
template <typename T>
auto
compHostFactorRho(      T*                                 osFactorRho,
                  const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock,
                  const int32_t                            bPosition,
                  const int32_t                            ePosition,
                  const int32_t                            braPrimGto,
                  const int32_t                            ketPrimGto) -> void
{
    // set up dimentsions
    
    const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();
    
    const auto ncpairs = ePosition - bPosition;
    
    // set up pointers to Xi values
        
    auto fxi = gtoPairBlock->getFactorsXi();
    
    // compute rho factors

    for (int32_t i = 0; i < ncpairs; i++)
    {
        const auto ioff = (bPosition + i) * nppairs;
            
        const auto bxi = fxi[ioff + braPrimGto];
        
        const auto kxi = fxi[ioff + ketPrimGto];
        
        osFactorRho[i] = bxi * kxi / (bxi + kxi);
    }
}

/**
Computes normalization factors for contracted GTOs batch.

@param osFactorNorm The Obara-Saika normalization factors.
@param gtoPairBlock The pointer to GTOs pairs block.
@param bPosition The start position of contracted GTOs batch.
@param ePosition The endposition of contracted GTOs batch.
@param braPrimGto The primitive GTO on bra side.
@param ketPrimGto The primitive GTO on ket side.
*/
template <typename T>
auto
compHostFactorNorm(      T*                                 osFactorNorm,
                   const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock,
                   const int32_t                            bPosition,
                   const int32_t                            ePosition,
                   const int32_t                            braPrimGto,
                   const int32_t                            ketPrimGto) -> void
{
    // set up dimentsions
    
    const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();
    
    const auto ncpairs = ePosition - bPosition;
    
    // set up pointers to Xi values
        
    auto fovl = gtoPairBlock->getOverlaps();
    
    // compute normalization factors
    
    auto fact = static_cast<T>((braPrimGto == ketPrimGto) ? 1.0 : 2.0);

    for (int32_t i = 0; i < ncpairs; i++)
    {
        const auto ioff = (bPosition + i) * nppairs;
        
        osFactorNorm[i] = fact * fovl[ioff + braPrimGto] * fovl[ioff + ketPrimGto];
    }
}

/**
Computes Boys function arguments for contracted GTOs batch.

@param bfArguments The Boys function arguments.
@param rDistancesPQ The R(PQ) = P - Q distances.
@param osFactorRho The Obara-Saika factor rho = xi * eta / (xi + eta).
@param nBatchPairs The number of pairs in batch.
*/
template <typename T>
auto
compHostBoysArguments(      BufferHostX<T>&     bfArguments,
                      const BufferHostMY<T, 3>& rDistancesPQ,
                      const T*                  osFactorRho,
                      const int32_t             nBatchPairs) -> void
{
    // set up pointers to R(PQ) = P - Q distances
           
    auto rpqx = rDistancesPQ.data(0);
           
    auto rpqy = rDistancesPQ.data(1);
           
    auto rpqz = rDistancesPQ.data(2);
    
    // compute Boys function arguments
    
    auto bargs = bfArguments.data();
    
    #pragma omp simd aligned(bargs, osFactorRho, rpqx, rpqy, rpqz: VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        bargs[i] = osFactorRho[i] * (rpqx[i] * rpqx[i] + rpqy[i] * rpqy[i] + rpqz[i] * rpqz[i]);
    }
}

}  // derirec  namespace

#endif /* DiagEriRecFacts_hpp */
