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

#ifndef PackedGtoPairContainer_hpp
#define PackedGtoPairContainer_hpp

#include <cstdint>
#include <vector>
#include <string>
#include <sstream>

#include "StringFormat.hpp"
#include "PackedGtoPairBlock.hpp"

/**
 Class CPackedGtoPairContainer stores vector of packed GTOs pair block objects and provides set of
 methods for manipulating with basis function pairs of various angular momentum.

 @author Z. Rinkevicius
 */
template <typename T, typename B = mem::Host>
class CPackedGtoPairContainer
{
    /**
     The vector of packed GTOs pair block objects.
     */
    VPackedGtoPairBlocks<T, B> _gtoPairBlocks;

   public:
    /**
     Creates an empty packed GTOs  pair container object.
     */
    CPackedGtoPairContainer()
    
        : _gtoPairBlocks(VPackedGtoPairBlocks<T, B>())
    {
        
    }

    /**
     Creates a packed GTOs pair container object.

     @param gtoPairBlocks the vector of packed GTOs pair block objects.
     */
    CPackedGtoPairContainer(const VPackedGtoPairBlocks<T, B>& gtoPairBlocks)
    
        : _gtoPairBlocks(gtoPairBlocks)
    {
        
    }

    /**
     Creates a packed GTOs pair container object by copying other packed GTOs pair container object.

     @param source the packed GTOs pair container object.
     */
    CPackedGtoPairContainer(const CPackedGtoPairContainer<T, B>& source)
    
        : _gtoPairBlocks(source._gtoPairBlocks)
    {
        
    }

    /**
     Creates a packed GTOs pair container object by moving other packed GTOs pair container object.

     @param source the packed GTOs pair container object.
     */
    CPackedGtoPairContainer(CPackedGtoPairContainer<T, B>&& source) noexcept
    
        : _gtoPairBlocks(std::move(source._gtoPairBlocks))
    {
        
    }

    /**
     Destroys a packed GTOs pair container object.
     */
    ~CPackedGtoPairContainer()
    {
        
    }

    /**
     Assigns a packed GTOs pair container object by copying other packed GTOs pair container object.

     @param source the packed GTOs pair container object.
     */
    auto
    operator=(const CPackedGtoPairContainer<T, B>& source) -> CPackedGtoPairContainer<T, B>&
    {
        if (this == &source) return *this;

        _gtoPairBlocks = source._gtoPairBlocks;

        return *this;
    }

    /**
     Assigns a packed GTOs pair container object by moving other packed GTOs pair container object.

     @param source the packed GTOs pair container object.
     */
    auto
    operator=(CPackedGtoPairContainer<T, B>&& source) noexcept -> CPackedGtoPairContainer<T, B>&
    {
        if (this == &source) return *this;

        _gtoPairBlocks = std::move(source._gtoPairBlocks);

        return *this;
    }

    /**
     Compares packed GTOs pair container object with other packed GTOs pair container object.

     @param other the packed GTOs pair container object.
     @return true if packed GTOs pair container objects are equal, false otherwise.
     */
    auto
    operator==(const CPackedGtoPairContainer<T, B>& other) const -> bool
    {
        if (_gtoPairBlocks.size() != other._gtoPairBlocks.size()) return false;

        for (size_t i = 0; i < _gtoPairBlocks.size(); i++)
        {
            if (_gtoPairBlocks[i] != other._gtoPairBlocks[i]) return false;
        }

        return true;
    }

    /**
     Compares packed GTOs pair container object with other packed GTOs pair container object.

     @param other the packed GTOs pair container object.
     @return true if packed GTOs pair container objects are not equal, false otherwise.
     */
    auto
    operator!=(const CPackedGtoPairContainer<T, B>& other) const -> bool
    {
        return !(*this == other);
    }
    
    /**
    Adds packed GTOs pair block object to packed GTOs pair container object.

    @param gtoPairBlock the packed GTOs pair block object.
    */
    auto
    add(const CPackedGtoPairBlock<T, B>& gtoPairBlock) -> void
    {
        _gtoPairBlocks.push_back(gtoPairBlock);
    }
    
    /**
    Gets constant  pointer to packed GTOs pair container  object.

    @return the constant  pointer to packed GTOs pair container  object.
    */
    auto
    to_pointer() const -> const CPackedGtoPairContainer*
    {
        return this;
    }
    
    /**
     Gets number of packed GTOs pair blocks in packed GTOs pair container object.

     @return the number of packed GTOs pair blocks.
     */
    auto
    getNumberOfBlocks() const -> int32_t
    {
        return static_cast<int32_t>(_gtoPairBlocks.size());
    }

    /**
     Gets specific packed GTOs block object from packed GTOs apir container object.

     @param iBlock the index of packed GTOs pair block object in packed GTOs pair container.
     @return the packed GTOs pair block object.
     */
    auto
    getBlock(const int32_t iBlock) const -> CPackedGtoPairBlock<T, B>
    {
        return _gtoPairBlocks[iBlock];
    }

    /**
     Prints summary of packed GTOs pair container to output string.

     @return the output string.
     */
    auto
    printSummary() const -> std::string
    {
        std::stringstream ss;

        ss << "Packed GTOs Pair Container\n";

        ss << std::string(26, '=') << "\n\n";

        ss << "Size: " << _gtoPairBlocks.size() << "\n\n";

        for (const auto& gpblock : _gtoPairBlocks)
        {
            ss << gpblock.printSummary() << std::endl;
        }
            
        return ss.str();
    }
};

// short hand notation for vector of packed GTOs pair container objects

template <typename T, typename B = mem::Host>
using VPackedGtoPairContainers = std::vector<CPackedGtoPairContainer<T, B>>;


#endif /* PackedGtoPairContainer_hpp */
