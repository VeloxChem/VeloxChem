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

#ifndef BinnedGtoPairBlock_hpp
#define BinnedGtoPairBlock_hpp

#include <cstdint>

#include "Buffer.hpp"
#include "BinnedGtoBlock.hpp"
#include "AngularMomentum.hpp"
#include "MathConst.hpp"

/**
 Class CBinnedGtoPairBlock stores data about binned basis function pairs and provides set of
 methods for manipulating with binned basis function pairs.

 @author Z. Rinkevicius
 */
template <typename T, typename B = mem::Host>
class CBinnedGtoPairBlock
{
    /**
     The angular momentum of bra side.
     */
    int32_t _braAngularMomentum;

    /**
     The angular momentum of ket side.
     */
    int32_t _ketAngularMomentum;
    
    /**
     The number of primitive GTO pairs in single contracted pair.
     */
    int32_t _nPrimitivePairs;

    /**
     The angular indexes of contracted pairs.
     */
    BufferXY<int32_t, B> _pairIndexes;
    
    /**
     The atomic indexes of contracted pairs.
     */
    BufferMY<int32_t, B, 2> _atomicIndexes;

    /**
     The primitives Gaussian function pair factors data (overlap, P coordinates,
     PA and PB distances, etc).
     */
    BufferMY<T, B, 24> _pairFactors;

   public:
    /**
     Creates an empty binned GTOs pair block object.
     */
    CBinnedGtoPairBlock()
    
        : _braAngularMomentum(-1)

        , _ketAngularMomentum(-1)

        , _nPrimitivePairs(0)

        , _pairIndexes(BufferXY<int32_t, B>())

        , _atomicIndexes(BufferMY<int32_t, B, 2>())

        , _pairFactors(BufferMY<T, B, 24>())
    {
        
    }
    
    /**
    Creates a binned GTOs pairs block object.

    @param braAngularMomentum the angular momentum on bra side.
    @param ketAngularMomentum the angular momentum on ket side.
    @param nPrimitivePairs the number of primitive GTO pairs in contracted GTO pair.
    @param pairIndexes the angular indexes of GTO pairs.
    @param atomicIndexes the  atomic indexes of GTO pairs.
    @param pairFactors the various primitive GTO pair factors (overlap, P coordinates,
           PA and PB distances, etc)
    */
    CBinnedGtoPairBlock(const int32_t                  braAngularMomentum,
                        const int32_t                  ketAngularMomentum,
                        const int32_t                  nPrimitivePairs,
                        const BufferXY<int32_t, B>&    pairIndexes,
                        const BufferMY<int32_t, B, 2>& atomicIndexes,
                        const BufferMY<T, B, 24>&      pairFactors)
    
        : _braAngularMomentum(braAngularMomentum)

        , _ketAngularMomentum(ketAngularMomentum)

        , _nPrimitivePairs(nPrimitivePairs)

        , _pairIndexes(pairIndexes)

        , _atomicIndexes(atomicIndexes)

        , _pairFactors(pairFactors)
    {
        
    }
    
    /**
     Creates a binned GTOs pairs block object.

     @param braGtoBlock the binned GTO block object for bra side.
     @param ketGtoBlock the binned GTO block object for ket side.
     */
    CBinnedGtoPairBlock(const CBinnedGtoBlock<T, B>& braGtoBlock,
                        const CBinnedGtoBlock<T, B>& ketGtoBlock)
    {
        // set up dimensions of GTO pairs
            
        _braAngularMomentum = braGtoBlock.getAngularMomentum();
            
        _ketAngularMomentum = ketGtoBlock.getAngularMomentum();
            
        const auto bpgtos = braGtoBlock.getNumberOfPrimGtos();
            
        const auto kpgtos = ketGtoBlock.getNumberOfPrimGtos();
            
        _nPrimitivePairs = bpgtos * kpgtos;
            
        // allocate GTO pairs data
            
        const auto bcgtos = braGtoBlock.getNumberOfContrGtos();
            
        const auto kcgtos = ketGtoBlock.getNumberOfContrGtos();
            
        const auto bcomps = angmom::to_SphericalComponents(_braAngularMomentum);
            
        const auto kcomps = angmom::to_SphericalComponents(_ketAngularMomentum);
           
        _pairIndexes = BufferXY<int32_t, B>(bcomps + kcomps, bcgtos * kcgtos);
            
        _atomicIndexes = BufferMY<int32_t, B, 2>(bcgtos * kcgtos);
            
        _pairFactors = BufferMY<T, B, 24>(bcgtos * kcgtos * _nPrimitivePairs);
            
        // set up pointers to GTOs block data on bra side

        auto bexp = braGtoBlock.getExponents();

        auto bnorm = braGtoBlock.getNormFactors();

        auto brx = braGtoBlock.getCoordinatesX();

        auto bry = braGtoBlock.getCoordinatesY();

        auto brz = braGtoBlock.getCoordinatesZ();
            
        auto batmids = braGtoBlock.getAtomicIdentifiers();
            
        // set up pointers to GTOs block data on ket side

        auto kexp = ketGtoBlock.getExponents();

        auto knorm = ketGtoBlock.getNormFactors();

        auto krx = ketGtoBlock.getCoordinatesX();

        auto kry = ketGtoBlock.getCoordinatesY();

        auto krz = ketGtoBlock.getCoordinatesZ();
            
        auto katmids = ketGtoBlock.getAtomicIdentifiers();
            
        // set up pointers to Obara-Saika prefactors

        auto ppfx = _pairFactors.data(0);

        auto ppfi = _pairFactors.data(1);

        auto ppfz = _pairFactors.data(2);

        auto ppss = _pairFactors.data(3);

        // set up pointers to P centers

        auto pprpx = _pairFactors.data(4);

        auto pprpy = _pairFactors.data(5);

        auto pprpz = _pairFactors.data(6);

        // set up pointers to R(PA) distances

        auto pprpax = _pairFactors.data(7);

        auto pprpay = _pairFactors.data(8);

        auto pprpaz = _pairFactors.data(9);

        // set up pointers to R(PB) distances

        auto pprpbx = _pairFactors.data(10);

        auto pprpby = _pairFactors.data(11);

        auto pprpbz = _pairFactors.data(12);

        // set up pointers to A centers

        auto pprax = _pairFactors.data(13);

        auto ppray = _pairFactors.data(14);

        auto ppraz = _pairFactors.data(15);

        // set up pointers to B centers

        auto pprbx = _pairFactors.data(16);

        auto pprby = _pairFactors.data(17);

        auto pprbz = _pairFactors.data(18);

        // set up pointers to R(AB) distances

        auto pprabx = _pairFactors.data(19);

        auto ppraby = _pairFactors.data(20);

        auto pprabz = _pairFactors.data(21);
            
        // set up pointers to exponents
            
        auto ppbexp = _pairFactors.data(22);
            
        auto ppkexp = _pairFactors.data(23);
            
        // set up pointers to atomic identifiers
            
        auto ppbaids = _atomicIndexes.data(0);
            
        auto ppkaids = _atomicIndexes.data(1);
            
        // loop over bra and ket GTOs blocks
            
        int32_t icgto = 0;
            
        for (int32_t i = 0; i < bcgtos; i++)
        {
            const auto koff = i * bpgtos;
                
            for (int32_t j = 0; j < kcgtos; j++)
            {
                const auto loff = j * kpgtos;
                    
                // loop over primitive GTO pairs
                    
                for (int32_t k = 0; k < bpgtos; k++)
                {
                    // A center coordinates

                    const auto rax = brx[koff + k];

                    const auto ray = bry[koff + k];

                    const auto raz = brz[koff + k];
                        
                    // A center primitive GTO

                    const auto fae = bexp[koff + k];

                    const auto fan = bnorm[koff + k];
                        
                    // primitive index
                        
                    const auto pgoff = (icgto * bpgtos + k) * kpgtos;

                    #pragma omp simd
                    for (int32_t l = 0; l < kpgtos; l++)
                    {
                        // store A center coordinates

                        pprax[pgoff + l] = rax;

                        ppray[pgoff + l] = ray;

                        ppraz[pgoff + l] = raz;

                        // store B center coordinates

                        pprbx[pgoff + l] = krx[loff + l];

                        pprby[pgoff + l] = kry[loff + l];

                        pprbz[pgoff + l] = krz[loff + l];

                        // store R(AB) distances

                        pprabx[pgoff + l] = rax - krx[loff + l];;

                        ppraby[pgoff + l] = ray - kry[loff + l];

                        pprabz[pgoff + l] = raz - krz[loff + l];
                            
                        // store exponents
                            
                        ppbexp[pgoff + l] = fae;
                            
                        ppkexp[pgoff + l] = kexp[loff + l];
                            
                        // store normalization factors
                            
                        ppss[pgoff + l] = fan * knorm[loff + l];
                    }
                }
                    
                // set up atomic indexes
                    
                ppbaids[icgto] = batmids[i];
                    
                ppkaids[icgto] = katmids[j];
                    
                // update contracted GTO pairs counter
                    
                icgto++;
            }
        }
            
        // compute Obara-Saika factors

        const auto fpi = mathconst::getPiValue();
            
        const auto ntpairs = getNumberOfContrPairs() * _nPrimitivePairs;
            
        #pragma omp simd
        for (int32_t i = 0; i < ntpairs; i++)
        {
            const auto fx = ppbexp[i] + ppkexp[i];
                
            const auto fi = 1.0 / fx;
                
            const auto fz = ppbexp[i] * ppkexp[i] * fi;
                
            // store Obara-Saika factors
                
            ppfx[i] = fx;
                
            ppfi[i] = fi;
                
            ppfz[i] = fz;
                
            ppss[i] *= std::pow(fpi * fi, 1.50)
                
                     * std::exp(-fz * (pprabx[i] * pprabx[i] +
                                          
                                       ppraby[i] * ppraby[i] +
                                          
                                       pprabz[i] * pprabz[i]));
                
            // store P center coordinates
                
            pprpx[i] = fi * (ppbexp[i] * pprax[i] + ppkexp[i] * pprbx[i]);
                
            pprpy[i] = fi * (ppbexp[i] * ppray[i] + ppkexp[i] * pprby[i]);
                
            pprpz[i] = fi * (ppbexp[i] * ppraz[i] + ppkexp[i] * pprbz[i]);
                
            // store R(PA) = P - A distances
                
            pprpax[i] = pprpx[i] - pprax[i];
                
            pprpay[i] = pprpy[i] - ppray[i];
                
            pprpaz[i] = pprpz[i] - ppraz[i];
                
            // store R(PB) = P - B distances
                
            pprpbx[i] = pprpx[i] - pprbx[i];
                
            pprpby[i] = pprpy[i] - pprby[i];
                
            pprpbz[i] = pprpz[i] - pprbz[i];
        }
            
        // set up angular identifiers for bra side

        for (int32_t i = 0; i < bcomps; i++)
        {
            auto gidx = braGtoBlock.getIdentifiers(i);

            auto bidx = _pairIndexes.data(i);
                
            int32_t idx = 0;
                
            for (int32_t j = 0; j < bcgtos; j++)
            {
                for (int32_t k = 0; k < kcgtos; k++)
                {
                    bidx[idx] = gidx[j];
                           
                    idx++;
                }
            }
        }

        // set up angular identifiers for ket side

        for (int32_t i = 0; i < kcomps; i++)
        {
            auto gidx = ketGtoBlock.getIdentifiers(i);

            auto kidx = _pairIndexes.data(bcomps + i);
                
            int32_t idx = 0;
                
            for (int32_t j = 0; j < bcgtos; j++)
            {
                for (int32_t k = 0; k < kcgtos; k++)
                {
                    kidx[idx] = gidx[k];
                           
                    idx++;
                }
            }
        }
    }
    
    /**
     Creates a binned GTOs pairs block object.

     @param gtoBlock the binned GTO block object for bra and ket sides.
     */
    CBinnedGtoPairBlock(const CBinnedGtoBlock<T, B>& gtoBlock)
    {
        // set up dimensions of GTO pairs
           
        _braAngularMomentum = gtoBlock.getAngularMomentum();
           
        _ketAngularMomentum = gtoBlock.getAngularMomentum();
           
        const auto pgtos = gtoBlock.getNumberOfPrimGtos();
           
        _nPrimitivePairs = pgtos * pgtos;
           
        // allocate GTO pairs data
           
        const auto cgtos = gtoBlock.getNumberOfContrGtos();
          
        const auto acomps = angmom::to_SphericalComponents(_braAngularMomentum);
          
        const auto npairs = cgtos * (cgtos + 1) / 2;
        
        _pairIndexes = BufferXY<int32_t, B>(2 * acomps, npairs);
           
        _atomicIndexes = BufferMY<int32_t, B, 2>(npairs);
           
        _pairFactors = BufferMY<T, B, 24>(npairs * _nPrimitivePairs);
        
        // set up pointers to GTOs block data on bra side

        auto gexp = gtoBlock.getExponents();

        auto gnorm = gtoBlock.getNormFactors();

        auto grx = gtoBlock.getCoordinatesX();

        auto gry = gtoBlock.getCoordinatesY();

        auto grz = gtoBlock.getCoordinatesZ();
        
        auto gatmids = gtoBlock.getAtomicIdentifiers();
        
        // set up pointers to Obara-Saika prefactors

        auto ppfx = _pairFactors.data(0);

        auto ppfi = _pairFactors.data(1);

        auto ppfz = _pairFactors.data(2);

        auto ppss = _pairFactors.data(3);

        // set up pointers to P centers

        auto pprpx = _pairFactors.data(4);

        auto pprpy = _pairFactors.data(5);

        auto pprpz = _pairFactors.data(6);

        // set up pointers to R(PA) distances

        auto pprpax = _pairFactors.data(7);

        auto pprpay = _pairFactors.data(8);

        auto pprpaz = _pairFactors.data(9);

        // set up pointers to R(PB) distances

        auto pprpbx = _pairFactors.data(10);

        auto pprpby = _pairFactors.data(11);

        auto pprpbz = _pairFactors.data(12);

        // set up pointers to A centers

        auto pprax = _pairFactors.data(13);

        auto ppray = _pairFactors.data(14);

        auto ppraz = _pairFactors.data(15);

        // set up pointers to B centers

        auto pprbx = _pairFactors.data(16);

        auto pprby = _pairFactors.data(17);

        auto pprbz = _pairFactors.data(18);

        // set up pointers to R(AB) distances

        auto pprabx = _pairFactors.data(19);

        auto ppraby = _pairFactors.data(20);

        auto pprabz = _pairFactors.data(21);
        
        // set up pointers to exponents
        
        auto ppbexp = _pairFactors.data(22);
        
        auto ppkexp = _pairFactors.data(23);
        
        // set up pointers to atomic identifiers
        
        auto ppbaids = _atomicIndexes.data(0);
        
        auto ppkaids = _atomicIndexes.data(1);
        
        // loop over bra and ket GTOs blocks
        
        int32_t icgto = 0;
        
        for (int32_t i = 0; i < cgtos; i++)
        {
            const auto koff = i * pgtos;
            
            for (int32_t j = i; j < cgtos; j++)
            {
                const auto loff = j * pgtos;
                
                // loop over primitive GTO pairs
                
                for (int32_t k = 0; k < pgtos; k++)
                {
                    // A center coordinates

                    const auto rax = grx[koff + k];

                    const auto ray = gry[koff + k];

                    const auto raz = grz[koff + k];
                    
                    // A center primitive GTO

                    const auto fae = gexp[koff + k];

                    const auto fan = gnorm[koff + k];
                    
                    // primitive index
                    
                    const auto pgoff = (icgto * pgtos + k) * pgtos;

                    #pragma omp simd
                    for (int32_t l = 0; l < pgtos; l++)
                    {
                        // store A center coordinates

                        pprax[pgoff + l] = rax;

                        ppray[pgoff + l] = ray;

                        ppraz[pgoff + l] = raz;

                        // store B center coordinates

                        pprbx[pgoff + l] = grx[loff + l];

                        pprby[pgoff + l] = gry[loff + l];

                        pprbz[pgoff + l] = grz[loff + l];

                        // store R(AB) distances

                        pprabx[pgoff + l] = rax - grx[loff + l];

                        ppraby[pgoff + l] = ray - gry[loff + l];

                        pprabz[pgoff + l] = raz - grz[loff + l];
                        
                        // store exponents
                        
                        ppbexp[pgoff + l] = fae;
                        
                        ppkexp[pgoff + l] = gexp[loff + l];
                        
                        // store normalization factors
                        
                        ppss[pgoff + l] = fan * gnorm[loff + l];
                    }
                }
                
                // set up atomic indexes
                
                ppbaids[icgto] = gatmids[i];
                
                ppkaids[icgto] = gatmids[j];
                
                // update contracted GTO pairs counter
                
                icgto++;
            }
        }
        
        // compute Obara-Saika factors

        const auto fpi = mathconst::getPiValue();
        
        const auto ntpairs = getNumberOfContrPairs() * _nPrimitivePairs;
        
        #pragma omp simd
        for (int32_t i = 0; i < ntpairs; i++)
        {
            const auto fx = ppbexp[i] + ppkexp[i];
            
            const auto fi = 1.0 / fx;
            
            const auto fz = ppbexp[i] * ppkexp[i] * fi;
            
            // store Obara-Saika factors
            
            ppfx[i] = fx;
            
            ppfi[i] = fi;
            
            ppfz[i] = fz;
            
            ppss[i] *= std::pow(fpi * fi, 1.50)
            
                     * std::exp(-fz * (pprabx[i] * pprabx[i] +
                                      
                                       ppraby[i] * ppraby[i] +
                                      
                                       pprabz[i] * pprabz[i]));
            
            // store P center coordinates
            
            pprpx[i] = fi * (ppbexp[i] * pprax[i] + ppkexp[i] * pprbx[i]);
            
            pprpy[i] = fi * (ppbexp[i] * ppray[i] + ppkexp[i] * pprby[i]);
            
            pprpz[i] = fi * (ppbexp[i] * ppraz[i] + ppkexp[i] * pprbz[i]);
            
            // store R(PA) = P - A distances
            
            pprpax[i] = pprpx[i] - pprax[i];
            
            pprpay[i] = pprpy[i] - ppray[i];
            
            pprpaz[i] = pprpz[i] - ppraz[i];
            
            // store R(PB) = P - B distances
            
            pprpbx[i] = pprpx[i] - pprbx[i];
            
            pprpby[i] = pprpy[i] - pprby[i];
            
            pprpbz[i] = pprpz[i] - pprbz[i];
        }
        
        // set up angular identifiers for bra and ket sides

        for (int32_t i = 0; i < acomps; i++)
        {
            auto gidx = gtoBlock.getIdentifiers(i);

            auto bidx = _pairIndexes.data(i);
            
            auto kidx = _pairIndexes.data(acomps + i);
            
            int32_t idx = 0;
            
            for (int32_t j = 0; j < cgtos; j++)
            {
                for (int32_t k = j; k < cgtos; k++)
                {
                    bidx[idx] = gidx[j];
                    
                    kidx[idx] = gidx[k];
                    
                    idx++;
                }
            }
        }
    }

    /**
     Creates a binned GTOs pair block object by copying other binned GTOs pair block object.

     @param source the binned GTOs pair block object.
     */
    CBinnedGtoPairBlock(const CBinnedGtoPairBlock<T, B>& source)
    
        : _braAngularMomentum(source._braAngularMomentum)

        , _ketAngularMomentum(source._ketAngularMomentum)

        , _nPrimitivePairs(source._nPrimitivePairs)

        , _pairIndexes(source._pairIndexes)

        , _atomicIndexes(source._atomicIndexes)

        , _pairFactors(source._pairFactors)

    {
        
    }

    /**
     Creates a binned GTOs pair block object by moving other binned GTOs pair block object.

     @param source the binned GTOs pair block object.
     */
    CBinnedGtoPairBlock(CBinnedGtoPairBlock<T, B>&& source) noexcept
    
        : _braAngularMomentum(std::move(source._braAngularMomentum))

        , _ketAngularMomentum(std::move(source._ketAngularMomentum))

        , _nPrimitivePairs(std::move(source._nPrimitivePairs))

        , _pairIndexes(std::move(source._pairIndexes))

        , _atomicIndexes(std::move(source._atomicIndexes))

        , _pairFactors(std::move(source._pairFactors))
    {
        
    }

    /**
     Destroys a binned GTOs pair block object.
     */
    ~CBinnedGtoPairBlock()
    {
        
    }

    /**
     Assigns a binned GTOs pair block object by copying other binned GTOs pair block object.

     @param source the binned GTOs pair block object.
     */
    auto
    operator=(const CBinnedGtoPairBlock<T, B>& source) -> CBinnedGtoPairBlock<T, B>&
    {
        if (this == &source) return *this;

        _braAngularMomentum = source._braAngularMomentum;

        _ketAngularMomentum = source._ketAngularMomentum;
        
        _nPrimitivePairs = source._nPrimitivePairs;

        _pairIndexes = source._pairIndexes;

        _atomicIndexes = source._atomicIndexes;

        _pairFactors = source._pairFactors;

        return *this;
    }

    /**
     Assigns a binned GTOs pair block object by moving other binned GTOs pair block object.

     @param source the binned GTOs pair block object.
     */
    auto
    operator=(CBinnedGtoPairBlock<T, B>&& source) noexcept -> CBinnedGtoPairBlock<T, B>&
    {
        if (this == &source) return *this;

        _braAngularMomentum = std::move(source._braAngularMomentum);

        _ketAngularMomentum = std::move(source._ketAngularMomentum);
        
        _nPrimitivePairs = std::move(source._nPrimitivePairs);

        _pairIndexes = std::move(source._pairIndexes);

        _atomicIndexes = std::move(source._atomicIndexes);

        _pairFactors = std::move(source._pairFactors);

        return *this;
    }

    /**
     Compares binned GTOs pair block object with other binned GTOs pair block object.

     @param other the binned GTOs pair block object.
     @return true if binned GTOs pair block objects are equal, false otherwise.
     */
    auto
    operator==(const CBinnedGtoPairBlock<T, B>& other) const -> bool
    {
        if (_braAngularMomentum != other._braAngularMomentum) return false;

        if (_ketAngularMomentum != other._ketAngularMomentum) return false;
        
        if (_nPrimitivePairs != other._nPrimitivePairs) return false;

        if (_pairIndexes != other._pairIndexes) return false;

        if (_atomicIndexes != other._atomicIndexes) return false;

        if (_pairFactors != other._pairFactors) return false;

        return true;
    }

    /**
     Compares binned GTOs pair block object with other binned GTOs pair block object.

     @param other the binned GTOs pair block object.
     @return true if binned GTOs pair block objects are not equal, false otherwise.
     */
    auto
    operator!=(const CBinnedGtoPairBlock<T, B>& other) const -> bool
    {
        return !(*this == other);
    }
    
    /**
    Compares binned GTOs pair block object with other binned GTOs pair block object.

    @param indexes the vector of indexes.
    @param ivalue the value of index to be selected.
    @return the binned GTOs pair block object.
    */
    auto
    compress(const BufferX<int32_t, B>& indexes,
             const int32_t              ivalue) const -> CBinnedGtoPairBlock<T, B>
    {
        // determine number of contracted GTO pairs
        
        int32_t ncgtos = 0;
        
        const auto ndims = (indexes.getMDSpan()).extent(0);
        
        for (size_t i = 0; i < ndims; i++)
        {
            if (indexes(i) == ivalue) ncgtos++;
        }
        
        // copy selectedd GTO pairs data
        
        if (ncgtos > 0)
        {
            // allocate compresed data buffers
            
            const auto bkcomps = (_pairIndexes.getMDSpan()).extent(0);
            
            BufferXY<int32_t, B> pairids(bkcomps, ncgtos);
            
            BufferMY<int32_t, B, 2> atomids(ncgtos);
            
            BufferMY<T, B, 24> pairdat(ncgtos * _nPrimitivePairs);
            
            // loop over GTO pairs data
            
            int32_t ipgto = 0;
            
            int32_t icgto = 0;
            
            for (size_t i = 0; i < ndims; i++)
            {
                if (indexes(i) == ivalue)
                {
                    // copy primitibe GTO pairs data
                    
                    const auto koff = i * _nPrimitivePairs;
                    
                    for (int32_t j = 0; j < 24; j++)
                    {
                        for (int32_t k = 0; k < _nPrimitivePairs; k++)
                        {
                            pairdat(j, ipgto + k) = _pairFactors(j, koff + k);
                        }
                    }
                    
                    // copy angular indexes
                    
                    for (int32_t j = 0; j < static_cast<int32_t>(bkcomps); j++)
                    {
                        pairids(j, icgto) = _pairIndexes(j, i);
                    }
                    
                    // copy atomic indexes
                    
                    for (int32_t j = 0; j < 2; j++)
                    {
                        atomids(j, icgto) = _atomicIndexes(j, i);
                    }
                    
                    // update GTO pairs counter
                    
                    ipgto += _nPrimitivePairs;
                    
                    icgto++;
                }
            }
            
            return CBinnedGtoPairBlock<T, B>(_braAngularMomentum, _ketAngularMomentum, _nPrimitivePairs,
                                             pairids, atomids, pairdat);
        }
        
        return CBinnedGtoPairBlock<T, B>();
    }
    
    /**
     Gets angular momentum of bra side in binned GTOs pair block object.

     @return the angular momentum of bra side.
     */
    auto
    getBraAngularMomentum() const -> int32_t
    {
        return _braAngularMomentum;
    }

    /**
     Gets angular momentum of ket side in binned GTOs pair block object.

     @return the angular momentum of ket side.
     */
    auto
    getKetAngularMomentum() const -> int32_t
    {
        return _ketAngularMomentum;
    }
    
    /**
     Gets number of contracted GTO pairs of binned GTOs pair block object.

     @return the angular momentum of bra side.
     */
    auto
    getNumberOfContrPairs() const -> int32_t
    {
        if (_atomicIndexes.empty())
        {
            return 0; 
        }
        else
        {
            return  static_cast<int32_t>((_atomicIndexes.getMDSpan()).extent(1));
        }
    }
    
    /**
     Gets number of primitive pairs in contracted GTO pair of binned GTOs pair block object.

     @return the angular momentum of bra side.
     */
    auto
    getNumberOfPrimPairs() const -> int32_t
    {
        return _nPrimitivePairs;
    }
    
    /**
    Gets number of spherical components  in contracted GTO pair of binned GTOs pair block object.

    @return the number of spherical components.
    */
    auto
    getNumberOfComponents() const -> int32_t
    {
        return angmom::to_SphericalComponents(_braAngularMomentum) *
            
               angmom::to_SphericalComponents(_ketAngularMomentum);
    }

    /**
     Checks if binned GTOs pair block is empty.

     @return true if binned GTOs pair block is empty, false otherwise.
     */
    auto
    empty() const -> bool
    {
        return _atomicIndexes.empty();
    }

    /**
     Gets constant vector of Obara-Saika factors: alpha_a + alpha_b.

     @return the vector of Obara-Saika factors.
     */
    auto
    getFactorsXi() const -> const T*
    {
        return _pairFactors.data(0);
    }

    /**
     Gets constant vector of Obara-Saika factors: 1 / (alpha_a + alpha_b).

     @return the vector of Obara-Saika factors.
     */
    auto
    getFactorsOneOverXi() const -> const T*
    {
        return _pairFactors.data(1);
    }

    /**
     Gets constant vector of Obara-Saika factors: alpha_a * alpha_b / (alpha_a + alpha_b).

     @return the vector of Obara-Saika factors.
     */
    auto
    getFactorsZeta() const -> const T*
    {
        return _pairFactors.data(2);
    }

    /**
     Gets constant vector of primitive overlaps c_mu * c_nu  (s_mu | s_nu).

     @return the vector of primitive overlaps.
     */
    auto
    getOverlaps() const -> const T*
    {
        return _pairFactors.data(3);
    }

    /**
     Gets constant vector of Cartesian X coordinates of Gaussian product center P.

     @return the vector of Cartesian X coordinates.
     */
    auto
    getCoordinatesPX() const -> const T*
    {
        return _pairFactors.data(4);
    }

    /**
     Gets constant vector of Cartesian Y coordinates of Gaussian product center P.

     @return the vector of Cartesian Y coordinates.
     */
    auto
    getCoordinatesPY() const -> const T*
    {
        return _pairFactors.data(5);
    }

    /**
     Gets constant vector of Cartesian Z coordinates of Gaussian product center P.

     @return the vector of Cartesian Z coordinates.
     */
    auto
    getCoordinatesPZ() const -> const T*
    {
        return _pairFactors.data(6);
    }

    /**
     Gets constant vector of Cartesian X component of R(PA) = P - A distances.

     @return the vector of Cartesian R(PA)_X distances.
     */
    auto
    getDistancesPAX() const -> const T*
    {
        return _pairFactors.data(7);
    }

    /**
     Gets constant vector of Cartesian Y component of R(PA) = P - A distances.

     @return the vector of Cartesian R(PA)_Y distances.
     */
    auto
    getDistancesPAY() const -> const T*
    {
        return _pairFactors.data(8);
    }

    /**
     Gets constant vector of Cartesian Z component of R(PA) = P - A distances.

     @return the vector of Cartesian R(PA)_Y distances.
     */
    auto
    getDistancesPAZ() const -> const T*
    {
        return _pairFactors.data(9);
    }

    /**
     Gets constant vector of Cartesian X component of R(PB) = P - B distances.

     @return the vector of Cartesian R(PB)_Z distances.
     */
    auto
    getDistancesPBX() const -> const T*
    {
        return _pairFactors.data(10);
    }

    /**
     Gets vector of Cartesian Y component of R(PB) = P - B distances.

     @return the vector of Cartesian R(PB)_Y distances.
     */
    auto
    getDistancesPBY() const -> const T*
    {
        return _pairFactors.data(11);
    }

    /**
     Gets constant vector of Cartesian Z component of R(PB) = P - B distances.

     @return the vector of Cartesian R(PB)_Z distances.
     */
    auto
    getDistancesPBZ() const -> const T*
    {
        return _pairFactors.data(12);
    }

    /**
     Gets constant vector of Cartesian X coordinates of center A.

     @return the vector of Cartesian X coordinates.
     */
    auto
    getCoordinatesAX() const -> const T*
    {
        return _pairFactors.data(13);
    }

    /**
     Gets constant vector of Cartesian Y coordinates of center A.

     @return the vector of Cartesian Y coordinates.
     */
    auto
    getCoordinatesAY() const -> const T*
    {
        return _pairFactors.data(14);
    }

    /**
     Gets constant vector of Cartesian Z coordinates of center A.

     @return the vector of Cartesian Z coordinates.
     */
    auto
    getCoordinatesAZ() const -> const T*
    {
        return _pairFactors.data(15);
    }

    /**
     Gets constant vector of Cartesian X coordinates of center B.

     @return the vector of Cartesian X coordinates.
     */
    auto
    getCoordinatesBX() const -> const T*
    {
        return _pairFactors.data(16);
    }

    /**
     Gets constant vector of Cartesian Y coordinates of center B.

     @return the vector of Cartesian Y coordinates.
     */
    auto
    getCoordinatesBY() const -> const T*
    {
        return _pairFactors.data(17);
    }

    /**
     Gets constant vector of Cartesian Z coordinates of center B.

     @return the vector of Cartesian Z coordinates.
     */
    auto
    getCoordinatesBZ() const -> const T*
    {
        return _pairFactors.data(18);
    }

    /**
     Gets constant vector of Cartesian X component of R(AB) = A - B distances.

     @return the vector of Cartesian R(AB)_Z distances.
     */
    auto
    getDistancesABX() const -> const T*
    {
        return _pairFactors.data(19);
    }

    /**
     Gets vector of Cartesian Y component of R(AB) = A - B distances.

     @return the vector of Cartesian R(AB)_Y distances.
     */
    auto
    getDistancesABY() const -> const T*
    {
        return _pairFactors.data(20);
    }

    /**
     Gets constant vector of Cartesian Z component of R(AB) = A - B distances.

     @return the vector of Cartesian R(AB)_Z distances.
     */
    auto
    getDistancesABZ() const -> const T*
    {
        return _pairFactors.data(21);
    }

    /**
     Gets exponents of primitive GTOs on bra side.

     @return the vector of primitive GTOs exponents.
     */
    auto
    getBraExponents() const -> const T*
    {
        return _pairFactors.data(22);
    }
    
    /**
     Gets exponents of primitive GTOs on ket side.

     @return the vector of primitive GTOs exponents.
     */
    auto
    getKetExponents() const -> const T*
    {
        return _pairFactors.data(23);
    }
    
    /**
     Gets constant pointer to contracted pairs bra indexes in full AO basis for
     specific angular momentum component.

     @param iComponent the component of angular momentum.
     @return the bra indexes in full AO basis.
     */
    auto
    getBraIdentifiers(const int32_t iComponent) const -> const int32_t*
    {
        return _pairIndexes.data(iComponent);
    }

    /**
     Gets constant pointer to contracted pairs ket indexes in full AO basis for
     specific angular momentum component.

     @param iComponent the component of angular momentum.
     @return the ket indexes in full AO basis.
     */
    auto
    getKetIdentifiers(const int32_t iComponent) const -> const int32_t*
    {
        return _pairIndexes.data(angmom::to_SphericalComponents(_braAngularMomentum) + iComponent);
    }
    
    /**
     Gets constant pointer to atomic indexes of bra side.

     @return the atomic indexes of bra side.
     */
    auto
    getBraAtomicIdentifiers() const -> const int32_t*
    {
        return _atomicIndexes.data(0);
    }
    
    /**
     Gets constant pointer to atomic indexes of ket side.

     @return the atomic indexes of ket side.
     */
    auto
    getKetAtomicIdentifiers() const -> const int32_t*
    {
        return _atomicIndexes.data(1);
    }
    
    /**
    Gets vector of R(AB) = A - B distances for given range of contracted
    GTO pairs.
     
    @param iStartPair the starting contracted GTO pair.
    @param iEndPair the last contracted GTO pair.
    @return the vector of R(AB) = A - B distances.
    */
    auto
    getDistancesAB(const int32_t iStartPair,
                   const int32_t iEndPair) const -> BufferMY<T, B, 3>
    {
        BufferMY<T, B, 3> rab(iEndPair - iStartPair);
        
        // set up pointer to distances R(AB) = A - B
        
        auto abx = rab.data(0);
        
        auto aby = rab.data(1);
        
        auto abz = rab.data(2);
        
        // set up pointer to distances data
        
        auto pabx = getDistancesABX();
        
        auto paby = getDistancesABY();
        
        auto pabz = getDistancesABZ();
        
        // copy distances

        for (int32_t i = iStartPair; i < iEndPair; i++)
        {
            abx[i - iStartPair] = pabx[i * _nPrimitivePairs];
            
            aby[i - iStartPair] = paby[i * _nPrimitivePairs];
            
            abz[i - iStartPair] = pabz[i * _nPrimitivePairs];
        }
        
        return rab;
    }
};


#endif /* BinnedGtoPairBlock_hpp */
