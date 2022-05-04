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

#ifndef BinnedGtoContainer_hpp
#define BinnedGtoContainer_hpp

#include <cstdint>
#include <vector>
#include <string>
#include <sstream>

#include "BinnedGtoBlock.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "StringFormat.hpp"
#include "AngularMomentum.hpp"

/**
 Class CBinnedGtoContainer stores vector of binned GTOs block objects and provides set of
 methods for manipulating with basis functions of various angular momentum.

 @author Z. Rinkevicius
 */

template <typename T, typename B = mem::Host>
class CBinnedGtoContainer
{
    /**
     The vector of binned GTOs block objects.
     */
    std::vector<CBinnedGtoBlock<T, B>> _gtoBlocks;

   public:
    /**
     Creates an empty binned GTOs container object.
     */
    CBinnedGtoContainer()
    {
        
    }
    

    /**
     Creates a binned GTOs container object.

     @param gtoBlocks the vector of binned GTOs block objects.
     */
    CBinnedGtoContainer(const std::vector<CBinnedGtoBlock<T, B>>& gtoBlocks)
    
        : _gtoBlocks(gtoBlocks)
    {
        
    }

    /**
     Creates a binned GTOs container object.

     @param molecule the molecule.
     @param basis the molecular basis.
     */
    CBinnedGtoContainer(const CMolecule&       molecule,
                        const CMolecularBasis& basis)
    
        : CBinnedGtoContainer<T, B>(molecule, basis, 0, molecule.getNumberOfAtoms())
    {
        
    }

    /**
     Creates a binned GTOs container object from list of atoms in molecule.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param iAtom the index of first atom in list of atoms.
     @param nAtoms the number of atoms in list of atoms.
     */
    CBinnedGtoContainer(const CMolecule&       molecule,
                        const CMolecularBasis& basis,
                        const int32_t          iAtom,
                        const int32_t          nAtoms)
    {
        for (int32_t i = 0; i <= basis.getMaxAngularMomentum(); i++)
        {
            for (const auto& cnum : basis.getContractionDepths(i))
            {
                if (const CBinnedGtoBlock<T, B> gblock(molecule, basis, iAtom, nAtoms, i, cnum); !gblock.empty())
                {
                    _gtoBlocks.push_back(gblock);
                }
            }
        }
    }

    /**
     Creates a binned GTOs container object by copying other binned GTOs container object.

     @param source the binned GTOs container object.
     */
    CBinnedGtoContainer(const CBinnedGtoContainer<T, B>& source)
        
        : _gtoBlocks(source._gtoBlocks)
    {
        
    }

    /**
     Creates a binned GTOs container object by moving other binned GTOs container object.

     @param source the binned GTOs container object.
     */
    CBinnedGtoContainer(CBinnedGtoContainer<T, B>&& source) noexcept
    
        : _gtoBlocks(std::move(source._gtoBlocks))
    {
        
    }

    /**
     Destroys a binned GTOs container object.
     */
    ~CBinnedGtoContainer()
    {
        
    }

    /**
     Assigns a binned GTOs container object by copying other binned GTOs container object.

     @param source the binned GTOs container object.
     */
    auto
    operator=(const CBinnedGtoContainer& source) -> CBinnedGtoContainer<T, B>&
    {
        if (this == &source) return *this;

        _gtoBlocks = source._gtoBlocks;

        return *this;
    }

    /**
     Assigns a binned GTOs container object by moving other binned GTOs container object.

     @param source the binned GTOs container object.
     */
    auto
    operator=(CBinnedGtoContainer<T, B>&& source) noexcept -> CBinnedGtoContainer<T, B>&
    {
        if (this == &source) return *this;

        _gtoBlocks = std::move(source._gtoBlocks);

        return *this;
    }

    /**
     Compares binned GTOs container object with other binned GTOs container object.

     @param other the binned GTOs container object.
     @return true if binned GTOs container objects are equal, false otherwise.
     */
    auto
    operator==(const CBinnedGtoContainer<T, B>& other) const -> bool
    {
        if (_gtoBlocks.size() != other._gtoBlocks.size()) return false;

        for (size_t i = 0; i < _gtoBlocks.size(); i++)
        {
            if (_gtoBlocks[i] != other._gtoBlocks[i]) return false;
        }

        return true;
    }

    /**
     Compares binned GTOs container object with other binned GTOs container object.

     @param other the binned GTOs container object.
     @return true if binned GTOs container objects are not equal, false otherwise.
     */
    auto
    operator!=(const CBinnedGtoContainer<T, B>& other) const -> bool
    {
        return !(*this == other);
    }
    
    /**
     Gets constant pointer to this binned GTOs container object.

     @return the pointer to this binned GTOs container objects.
    */
    auto getPointer() const -> const CBinnedGtoContainer<T, B>*
    {
        return this;
    }

    /**
     Gets maximum angular momentum of GTOs block objects in GTOs container.

     @return the maxumum angular momentum.
     */
    auto
    getMaxAngularMomentum() const -> int32_t
    {
        int32_t maxmom = 0;

        for (const auto& gblock : _gtoBlocks)
        {
            if (const auto cmom = gblock.getAngularMomentum(); maxmom < cmom)
            {
                maxmom = cmom;
            }
        }

        return maxmom;
    }
    
    /**
     Gets angular momentum of specific binned GTOs block from binned GTOs container.

     @param iBlock the index of binned GTOs block in binned GTOs container.
     @return the angular momentum of  binned GTOs block.
     */
    auto
    getAngularMomentum(const int32_t iBlock) const -> int32_t
    {
        return _gtoBlocks[iBlock].getAngularMomentum();
    }
    
    /**
     Gets maximum number of contracted GTOs in GTOs block objects inside GTOs container.

     @return the maximum number of contracted GTOs in GTOs block object.
     */
    auto
    getMaxNumberOfContrGtos() const -> int32_t
    {
        int32_t maxgtos = 0;

        for (const auto& gblock : _gtoBlocks)
        {
            if (const auto cgtos = gblock.getNumberOfContrGtos(); maxgtos < cgtos)
            {
                maxgtos = cgtos;
            }
        }

        return maxgtos;
    }
    
    /**
     Gets number of atomic orbitals in  GTOs container.

     @return the number of atomic orbitals.
     */
    auto
    getNumberOfAtomicOrbitals() const -> int32_t
    {
        int32_t naos = 0;

        for (const auto& gblock : _gtoBlocks)
        {
            naos += angmom::to_SphericalComponents(gblock.getAngularMomentum())
                
                  * gblock.getNumberOfContrGtos();
        }

        return naos;
    }
    
    /**
     Gets maximum number of contracted GTOs in GTOs block with specific angular momentum
     inside GTOs container.

     @param angularMomentum the angular momentum of GTOs block.
     @return the maximum number of contracted GTOs pairs.
     */
    auto
    getMaxNumberOfContrGtos(const int32_t angularMomentum) const -> int32_t
    {
        int32_t maxgtos = 0;

        for (const auto& gblock : _gtoBlocks)
        {
            if (angularMomentum == gblock.getAngularMomentum())
            {
                if (const auto cgtos = gblock.getNumberOfContrGtos(); maxgtos < cgtos)
                {
                    maxgtos = cgtos;
                }
            }
        }

        return maxgtos;
    }

    /**
     Gets number of binned GTOs blocks in binned GTOs container.

     @return the number of binned GTOs blocks.
     */
    auto
    getNumberOfBinnedGtoBlocks() const -> int32_t
    {
        return static_cast<int32_t>(_gtoBlocks.size());
    }

    /**
     Gets specific binned GTOs block object from binned GTOs container.

     @param iBlock the index of binned GTOs block object in binned GTOs container.
     @return the binned GTOs block object.
     */
    auto
    getBinnedGtoBlock(const int32_t iBlock) const -> CBinnedGtoBlock<T, B>
    {
        return _gtoBlocks[iBlock];
    }
    
    /**
     Prints summary of binned GTOs blocks to output string.

     @return the output string.
     */
    auto
    printSummary() const -> std::string
    {
        std::stringstream ss;

        ss << "Binned GTOs Container\n";

        ss << std::string(21, '=') << "\n\n";

        ss << "Size: " << _gtoBlocks.size() << "\n\n";

        for (const auto& gblock : _gtoBlocks)
        {
            ss << "Binned GTOs Block: L = ";
                
            ss << fstr::to_AngularMomentum(gblock.getAngularMomentum());
                
            ss << " Number of GTOs: ";
                
            ss << fstr::to_string(gblock.getNumberOfContrGtos(), 8, fmt::right);
                
            ss << " Contraction Depth: ";
                
            ss << fstr::to_string(gblock.getNumberOfPrimGtos(), 4, fmt::right);
                
            ss << "\n";
        }
            
        return ss.str();
    }
};

#endif /* BinnedGtoContainer_hpp */
