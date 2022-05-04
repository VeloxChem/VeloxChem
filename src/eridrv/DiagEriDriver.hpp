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
    
    @param molecule The molecule.
    @param aoBasis The molecular AO basis.
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
                                    
                                std::cout << " Bra : " <<  gtopairs.getBraAngularMomentum();
                                
                                std::cout << " Ket : " <<  gtopairs.getKetAngularMomentum();
                                
                                std::cout << " No. contr. pairs: " << gtopairs.getNumberOfContrPairs();
                                
                                std::cout << " No. prim. pairs: " << gtopairs.getNumberOfPrimPairs() << std::endl;
                            }
                        }
                    }
                }
            }
        }
        
        return CPackedGtoPairContainer<T, B>();
    }
    
    
};


#endif /* DiagEriDriver_hpp */
