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
    // compute R(PQ) distances
    
    if (braPrimGto == ketPrimGto)
    {
        rDistancesPQ.setZero();
    }
    else
    {
        // set up dimentsions
        
        const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();
        
        const auto ncpairs = ePosition - bPosition;
        
        // set up pointers to R(PQ) = P - Q distances
               
        auto pqx = rDistancesPQ.data(0);
               
        auto pqy = rDistancesPQ.data(1);
               
        auto pqz = rDistancesPQ.data(2);
        
        // set up pointers to P center coordinates
            
        auto rpx = gtoPairBlock->getCoordinatesPX();
            
        auto rpy = gtoPairBlock->getCoordinatesPY();
            
        auto rpz = gtoPairBlock->getCoordinatesPZ();
        
        for (int32_t i = 0; i < ncpairs; i++)
        {
            const auto ioff = (bPosition + i) * nppairs;
            
            pqx[i] = rpx[ioff + braPrimGto] - rpx[ioff + ketPrimGto];
                               
            pqy[i] = rpy[ioff + braPrimGto] - rpy[ioff + ketPrimGto];
                               
            pqz[i] = rpz[ioff + braPrimGto] - rpz[ioff + ketPrimGto];
        }
    }
}

/**
Computes rho = xi * eta / (xi + eta) factors for contracted GTOs batch.

@param osFactorsRho The Obara-Saika factors rho = xi * eta / (xi + eta).
@param gtoPairBlock The pointer to GTOs pairs block.
@param bPosition The start position of contracted GTOs batch.
@param ePosition The endposition of contracted GTOs batch.
@param braPrimGto The primitive GTO on bra side.
@param ketPrimGto The primitive GTO on ket side.
*/
template <typename T>
auto
compHostFactorsRho(      T*                                 osFactorsRho,
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
    
    if (braPrimGto == ketPrimGto)
    {
        const auto fact = static_cast<T>(0.5);
        
        for (int32_t i = 0; i < ncpairs; i++)
        {
            const auto ioff = (bPosition + i) * nppairs + braPrimGto;

            osFactorsRho[i] = fact * fxi[ioff];
        }
    }
    else
    {
        for (int32_t i = 0; i < ncpairs; i++)
        {
            const auto ioff = (bPosition + i) * nppairs;
                
            const auto bxi = fxi[ioff + braPrimGto];
            
            const auto kxi = fxi[ioff + ketPrimGto];
            
            osFactorsRho[i] = bxi * kxi / (bxi + kxi);
        }
    }
}

/**
Computes normalization factors for contracted GTOs batch.

@param osFactorsNorm The Obara-Saika normalization factors.
@param gtoPairBlock The pointer to GTOs pairs block.
@param bPosition The start position of contracted GTOs batch.
@param ePosition The endposition of contracted GTOs batch.
@param braPrimGto The primitive GTO on bra side.
@param ketPrimGto The primitive GTO on ket side.
*/
template <typename T>
auto
compHostFactorsNorm(      T*                                 osFactorsNorm,
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
    
    const auto fact = static_cast<T>((braPrimGto == ketPrimGto) ? 1.0 : 2.0);

    for (int32_t i = 0; i < ncpairs; i++)
    {
        const auto ioff = (bPosition + i) * nppairs;
        
        osFactorsNorm[i] = fact * fovl[ioff + braPrimGto] * fovl[ioff + ketPrimGto];
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
           
    auto pqx = rDistancesPQ.data(0);
           
    auto pqy = rDistancesPQ.data(1);
           
    auto pqz = rDistancesPQ.data(2);
    
    // compute Boys function arguments
    
    auto bargs = bfArguments.data();
    
    #pragma omp simd aligned(bargs, osFactorRho, pqx, pqy, pqz: VLX_ALIGN)
    for (int32_t i = 0; i < nBatchPairs; i++)
    {
        bargs[i] = osFactorRho[i] * (pqx[i] * pqx[i] + pqy[i] * pqy[i] + pqz[i] * pqz[i]);
    }
}

/**
Computes zeta = 1.0 / (xi + eta) factors for contracted GTOs batch.

@param osFactorsZeta The Obara-Saika factors zeta = 1.0 / (xi + eta).
@param gtoPairBlock The pointer to GTOs pairs block.
@param bPosition The start position of contracted GTOs batch.
@param ePosition The endposition of contracted GTOs batch.
@param braPrimGto The primitive GTO on bra side.
@param ketPrimGto The primitive GTO on ket side.
*/
template <typename T>
auto
compHostFactorsZeta(     T*                                 osFactorsZeta,
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
    
    // compute zeta factors
    
    if (braPrimGto == ketPrimGto)
    {
        const auto fact = static_cast<T>(0.5);
        
        for (int32_t i = 0; i < ncpairs; i++)
        {
            const auto ioff = (bPosition + i) * nppairs + braPrimGto;
            
            osFactorsZeta[i] = fact / fxi[ioff];
        }
    }
    else
    {
        const auto fact = static_cast<T>(1.0);
        
        for (int32_t i = 0; i < ncpairs; i++)
        {
            const auto ioff = (bPosition + i) * nppairs;
            
            osFactorsZeta[i] = fact / (fxi[ioff + braPrimGto] + fxi[ioff + ketPrimGto]);
        }
    }
}

/**
Computes partial zeta = 1.0 / xi factors for contracted GTOs batch.

@param osFactorsPartZeta The Obara-Saika factors partial zeta = 1.0 / xi.
@param gtoPairBlock The pointer to GTOs pairs block.
@param bPosition The start position of contracted GTOs batch.
@param ePosition The endposition of contracted GTOs batch.
@param partPrimGto The primitive GTO.
*/
template <typename T>
auto
compHostFactorsPartialZeta(      T*                                 osFactorsPartZeta,
                           const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock,
                           const int32_t                            bPosition,
                           const int32_t                            ePosition,
                           const int32_t                            partPrimGto) -> void
{
    // set up dimentsions
    
    const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();
    
    const auto ncpairs = ePosition - bPosition;
    
    // set up pointers to Xi values
        
    auto fxi = gtoPairBlock->getFactorsXi();
    
    // compute zeta factors

    const auto fact = static_cast<T>(1.0);
    
    for (int32_t i = 0; i < ncpairs; i++)
    {
        const auto ioff = (bPosition + i) * nppairs;
        
        osFactorsPartZeta[i] = fact / fxi[ioff + partPrimGto];
    }
}

/**
Computes W = (P * xi  + Q * eta) / (xi + eta)  coordinates for contracted GTOs batch.

@param rCoordinatesW The W = (P * xi  + Q * eta) / (xi + eta)  coordinates.
@param gtoPairBlock The pointer to GTOs pairs block.
@param bPosition The start position of contracted GTOs batch.
@param ePosition The endposition of contracted GTOs batch.
@param braPrimGto The primitive GTO on bra side.
@param ketPrimGto The primitive GTO on ket side.
*/
template <typename T>
auto
compHostCoordinatesW(      BufferHostMY<T, 3>&                rCoordinatesW,
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
        
    // set up pointers to W coordinates
           
    auto rwx = rCoordinatesW.data(0);
           
    auto rwy = rCoordinatesW.data(1);
           
    auto rwz = rCoordinatesW.data(2);

    // compute W coordinates
    
    if (braPrimGto == ketPrimGto)
    {
        for (int32_t i = 0; i < ncpairs; i++)
        {
            const auto ioff = (bPosition + i) * nppairs + braPrimGto;
            
            rwx[i] = rpx[ioff];
                               
            rwy[i] = rpy[ioff];
                               
            rwz[i] = rpz[ioff];
        }
    }
    else
    {
        // set up pointers to Xi values
            
        auto fxi = gtoPairBlock->getFactorsXi();
        
        const auto fact = static_cast<T>(1.0);
        
        for (int32_t i = 0; i < ncpairs; i++)
        {
            const auto ioff = (bPosition + i) * nppairs;
            
            const auto bxi = fxi[ioff + braPrimGto];
            
            const auto kxi = fxi[ioff + ketPrimGto];
            
            const auto bkxi = fact / (bxi + kxi);
            
            rwx[i] = bkxi * (rpx[ioff + braPrimGto] * bxi + rpx[ioff + ketPrimGto] * kxi);
                               
            rwy[i] = bkxi * (rpy[ioff + braPrimGto] * bxi + rpy[ioff + ketPrimGto] * kxi);
                               
            rwz[i] = bkxi * (rpz[ioff + braPrimGto] * bxi + rpz[ioff + ketPrimGto] * kxi);
        }
    }
}

/**
Computes R(WT) = W - T distances for contracted GTOs batch.

@param rDistancesWT The R(WT) = W - T distances.
@param rCoordinatesW The W coordinates.
@param gtoPairBlock The pointer to GTOs pairs block.
@param bPosition The start position of contracted GTOs batch.
@param ePosition The endposition of contracted GTOs batch.
@param partPrimGto The primitive GTO.
*/
template <typename T>
auto
compHostDistancesWT(      BufferHostMY<T, 3>&                rDistancesWT,
                    const BufferHostMY<T, 3>&                rCoordinatesW,
                    const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock,
                    const int32_t                            bPosition,
                    const int32_t                            ePosition,
                    const int32_t                            partPrimGto) -> void
{
    // set up dimentsions
    
    const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();
    
    const auto ncpairs = ePosition - bPosition;
    
    // set up pointers to R(WT) = W - T distances
           
    auto wtx = rDistancesWT.data(0);
           
    auto wty = rDistancesWT.data(1);
           
    auto wtz = rDistancesWT.data(2);
    
    // set up pointers to W coordinates
           
    auto rwx = rCoordinatesW.data(0);
           
    auto rwy = rCoordinatesW.data(1);
           
    auto rwz = rCoordinatesW.data(2);
    
    // set up pointers to P center coordinates
        
    auto rpx = gtoPairBlock->getCoordinatesPX();
        
    auto rpy = gtoPairBlock->getCoordinatesPY();
        
    auto rpz = gtoPairBlock->getCoordinatesPZ();

    // compute R(WX) = W - T distances
    
    for (int32_t i = 0; i < ncpairs; i++)
    {
        const auto ioff = (bPosition + i) * nppairs + partPrimGto;
            
        wtx[i] = rwx[i] - rpx[ioff];
                               
        wty[i] = rwy[i] - rpy[ioff];
                               
        wtz[i] = rwz[i] - rpz[ioff];
    }
}

/**
Computes R(PT) = P - T distances for contracted GTOs batch.

@param rDistancesPT The R(PT) = P - T distances.
@param gtoPairBlock The pointer to GTOs pairs block.
@param bPosition The start position of contracted GTOs batch.
@param ePosition The endposition of contracted GTOs batch.
@param partPrimGto The primitive GTO.
*/
template <typename T>
auto
compHostDistancesPT(      BufferHostMY<T, 3>&                rDistancesPT,
                    const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock,
                    const int32_t                            bPosition,
                    const int32_t                            ePosition,
                    const int32_t                            partPrimGto) -> void
{
    // set up dimentsions
    
    const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();
    
    const auto ncpairs = ePosition - bPosition;
    
    // set up pointers to R(PT) = P - T distances
           
    auto ptx = rDistancesPT.data(0);
           
    auto pty = rDistancesPT.data(1);
           
    auto ptz = rDistancesPT.data(2);
    
    // set up pointers to R(PB) = P - B distances
        
    auto pbx = gtoPairBlock->getDistancesPBX();
        
    auto pby = gtoPairBlock->getDistancesPBY();
        
    auto pbz = gtoPairBlock->getDistancesPBZ();

    // compute R(PT) = P - T distances
    
    for (int32_t i = 0; i < ncpairs; i++)
    {
        const auto ioff = (bPosition + i) * nppairs + partPrimGto;
            
        ptx[i] = pbx[ioff];
                               
        pty[i] = pby[ioff];
                               
        ptz[i] = pbz[ioff];
    }
}

/**
Select maximum values from integrals buffer and stores into given array.

@param maxValues The array of maximum values.
@param intsBuffer The integrals buffer.
@param bPosition The start position of contracted GTOs batch.
@param nBatchPairs The number of pairs in batch.
*/
template <typename T, int32_t N>
auto
selectHostMaxValues(      T*                  maxValues,
                    const BufferHostMY<T, N>& intsBuffer,
                    const int32_t             bPosition,
                    const int32_t             nBatchPairs) -> void
{
    for (int32_t i = 0; i < N; i++)
    {
        auto tints = intsBuffer.data(i);
        
        for (int32_t j = 0; j < nBatchPairs; j++)
        {
            if (tints[j] > maxValues[bPosition + j])
            {
                maxValues[bPosition + j] = tints[j];
            }
        }
    }
}
    
}  // derirec  namespace

#endif /* DiagEriRecFacts_hpp */
