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

#ifndef DiagEriDriver_hpp
#define DiagEriDriver_hpp

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "PackedGtoPairContainer.hpp"
#include "BinnedGtoContainer.hpp"
#include "EriBlockDims.hpp"

#include <iostream>

/**
 Class CGtoPairsMapper maps packed GTOs pairs container.
 
 @author Z. Rinkevicius
 */
template <typename T, typename B = mem::Host>
class CDiagEriDriver
{
    
public:
    
    /**
     Creates a diagonal electron repulsion integrals driver object.
     */
    CDiagEriDriver()
    {
        
    }
    
    /**
     Destroys  a diagonal electron repulsion integrals driver object.
     */
    ~CDiagEriDriver()
    {
        
    }
    
    /**
    Computes screening integrals and connstructs packed GTOs pair container for given
    molecule and basis set.
    
    @param intsBuffer The integrals buffer.
    @param gtoPairBlock The pointer to GTOs pairs block.
    @param bPosition The start position of contracted GTOs batch.
    @param ePosition The endposition of contracted GTOs batch.
    */
    auto
    evaluate(      T*                         intsBuffer,
             const CBinnedGtoPairBlock<T, B>* gtoPairBlock,
             const int32_t                    bPosition,
             const int32_t                    ePosition) const -> void
    {
        
    }
    
    /**
    Computes screening integrals and connstructs packed GTOs pair container for given
    molecule and basis set.
    
    @param gtoPairBlock The GTOs pairs block.
    @return The packed GTOs pairs block.
    */
    auto
    compute(const CBinnedGtoPairBlock<T, B>& gtoPairBlock) const -> CPackedGtoPairBlock<T, B>
    {
        // set up dimension and angular momentum
        
        const auto bsize = eridims::getDiagBatchSize();
        
        const auto ncpairs = gtoPairBlock.getNumberOfContrPairs();
        
        const auto nppairs = gtoPairBlock.getNumberOfPrimPairs();
        
        const auto bang = gtoPairBlock.getBraAngularMomentum();
        
        const auto kang = gtoPairBlock.getKetAngularMomentum();
        
        std::cout << "(" << bang << "," << kang << "|" << bang << "," << kang << "): " << ncpairs << " <-> " << nppairs << std::endl;
        
        // allocate integrals buffer
        
        auto buffer = BufferX<T, B>::Zero(ncpairs);
        
        // set up pointers
        
        auto pbuffer = buffer.data();
        
        auto gtopairs = &gtoPairBlock;
        
        #pragma omp parallel shared(pbuffer, gtopairs, bsize, ncpairs)
        {
            #pragma omp single nowait
            {
                const auto nblocks =  ncpairs / bsize;
                
                for (int32_t i = 0; i < nblocks; i++)
                {
                    const auto bstart = i * bsize;
                    
                    const auto bend = bstart + bsize;
                    
                    #pragma omp task firstprivate(bstart, bend)
                    {
                        evaluate(pbuffer, gtopairs, bstart, bend);
                    }
                }
                
                if ((ncpairs % bsize) != 0)
                {
                    const auto bstart = nblocks * bsize;
                    
                    const auto bend = ncpairs;
                    
                    #pragma omp task firstprivate(bstart, bend)
                    {
                        evaluate(pbuffer, gtopairs, bstart, bend);
                    }
                }
            }
        }
        
        return CPackedGtoPairBlock<T, B>();
    }
    
    /**
    Computes screening integrals and connstructs packed GTOs pair container for given
    molecule and basis set.
    
     @param molecule The molecule.
     @param aoBasis The molecular AO basis.
     @return The packed GTOs pairs container.
    */
    auto
    compute(const CMolecule&       molecule,
            const CMolecularBasis& aoBasis) const -> CPackedGtoPairContainer<T, B>
    {
        // create binned GTOs container
        
        const CBinnedGtoContainer<T,B> bkgtos(molecule, aoBasis);
            
        const auto nblocks = bkgtos.getNumberOfBinnedGtoBlocks();
            
        const auto mang = bkgtos.getMaxAngularMomentum();
            
        // ordered loop over angular momentum pairs
        
        CPackedGtoPairContainer<T, B> gtodata;
            
        for (int32_t i = 0; i <= mang; i++)
        {
            for (int32_t j = i; j <= mang; j++)
            {
                for (int32_t k = 0; k < nblocks; k++)
                {
                    if (bkgtos.getAngularMomentum(k) == i)
                    {
                        const auto bgtos = bkgtos.getBinnedGtoBlock(k);
                            
                        for (int32_t l = k; l < nblocks; l++)
                        {
                            if (bkgtos.getAngularMomentum(l) == j)
                            {
                                const auto kgtos = bkgtos.getBinnedGtoBlock(l);
                                    
                                const auto gtopairs = (k == l) ? CBinnedGtoPairBlock<T, B>(bgtos)
                                    
                                                               : CBinnedGtoPairBlock<T, B>(bgtos, kgtos);
                                
                                gtodata.add(compute(gtopairs));
                            }
                        }
                    }
                }
            }
        }
        
        return gtodata;
    }
    
};


#endif /* DiagEriDriver_hpp */
