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

#ifndef PackedGtoPairBlock_hpp
#define PackedGtoPairBlock_hpp

#include <cstdint>
#include <array>
#include <vector>
#include <cmath>
#include <sstream>

#include "BinnedGtoPairBlock.hpp"
#include "Buffer.hpp"
#include "MathFunc.hpp"
#include "StringFormat.hpp"

/**
 Class CPackedGtoPairBlock stores data about packed basis function pairs and provides set of
 methods for manipulating with packed basis function pairs.

 @author Z. Rinkevicius
 */
template <typename T, typename B = mem::Host>
class CPackedGtoPairBlock
{
    /**
     The array of binned GTO pair blocks.
     */
    std::array<CBinnedGtoPairBlock<T, B>, 16> _gtoPairBlocks;

   public:
    /**
     Creates an empty packed GTOs pair block object.
     */
    CPackedGtoPairBlock()
        
        : _gtoPairBlocks(std::array<CBinnedGtoPairBlock<T, B>, 16>())
    {
        
    }
    
    /**
     Creates a packed GTO pair block object.

     @param gtoPairBlocks the vector of binned GTO pair block objects.
     @param indexes the vector of indexes to distribute the vector of binned GTO pair block objects.
     */
    CPackedGtoPairBlock(const std::vector<CBinnedGtoPairBlock<T, B>>& gtoPairBlocks,
                        const std::vector<int32_t>&                   indexes)
    
        : _gtoPairBlocks(std::array<CBinnedGtoPairBlock<T, B>, 16>())
    {
        for (size_t i = 0; i < indexes.size(); i++)
        {
            _gtoPairBlocks[indexes[i]] = gtoPairBlocks[i];
        }
    }
    
    /**
     Creates a packed GTO pair block object.

     @param gtoPairBlock the vector of binned GTO pair block objects.
     @param intsBuffer  the integrals buffer used to partition binned GTI.
     */
    CPackedGtoPairBlock(const CBinnedGtoPairBlock<T, B>& gtoPairBlock,
                        const BufferX<T, B>&             intsBuffer)
    {
        // compute index mapping
        
        const auto nints = (intsBuffer.getMDSpan()).extent(0);
       
        auto idxints = BufferX<int32_t, B>::Zero(nints);
            
        for (int32_t i = 0; i < nints; i++)
        {
            if (const auto tval = std::sqrt(intsBuffer(i)); tval < 1.0)
            {
                idxints(i) = static_cast<int32_t>(-std::log10(tval));
            }
        }
            
        // partition binned GTO pair block
            
        std::set<int32_t> primidx;
            
        for (int32_t i = 0; i < nints; i++)
        {
            primidx.insert(idxints(i));
        }
            
        for (const auto& i : primidx)
        {
            if (i < _gtoPairBlocks.size())
            {
                _gtoPairBlocks[i] = gtoPairBlock.compress(idxints, i);
            }
        }
    }

    /**
     Creates a packed GTOs pair block object by copying other packed GTOs pair block object.

     @param source the packed GTOs pair block object.
     */
    CPackedGtoPairBlock(const CPackedGtoPairBlock<T, B>& source)
    
        : _gtoPairBlocks(source._gtoPairBlocks)
    {
        
    }

    /**
     Creates a packed GTOs pair block object by moving other packed GTOs pair block object.

     @param source the packed GTOs pair block object.
     */
    CPackedGtoPairBlock(CPackedGtoPairBlock<T, B>&& source) noexcept
        
        : _gtoPairBlocks(std::move(source._gtoPairBlocks))
    {
        
    }

    /**
     Destroys a packed GTOs pair block object.
     */
    ~CPackedGtoPairBlock()
    {
        
    }

    /**
     Assigns a packed GTOs pair block object by copying other packed GTOs pair block object.

     @param source the packed GTOs pair block object.
     */
    auto
    operator=(const CPackedGtoPairBlock<T, B>& source) -> CPackedGtoPairBlock<T, B>&
    {
        if (this == &source) return *this;

        _gtoPairBlocks = source._gtoPairBlocks;

        return *this;
    }

    /**
     Assigns a packed GTOs pair block object by moving other packed GTOs pair block object.

     @param source the packed GTOs pair block object.
     */
    auto
    operator=(CPackedGtoPairBlock<T, B>&& source) noexcept -> CPackedGtoPairBlock<T, B>&
    {
        if (this == &source) return *this;

        _gtoPairBlocks = std::move(source._gtoPairBlocks);

        return *this;
    }

    /**
     Compares packed GTOs pair block object with other packed GTOs pair block object.

     @param other the packed GTOs pair block object.
     @return true if packed GTOs pair block objects are equal, false otherwise.
     */
    auto
    operator==(const CPackedGtoPairBlock<T, B>& other) const -> bool
    {
        return _gtoPairBlocks == other._gtoPairBlocks;
    }

    /**
     Compares packed GTOs pair block object with other packed GTOs pair block object.

     @param other the packed GTOs pair block object.
     @return true if packed GTOs pair block objects are not equal, false otherwise.
     */
    auto
    operator!=(const CPackedGtoPairBlock& other) const -> bool
    {
        return !(*this == other);
    }
    
    /**
     Determines number of contracted GTOs pairs subblocks.

     @return the number of  contracted GTOs pairs subblocks.
    */
    auto
    size() const -> int32_t
    {
        return static_cast<int32_t>(_gtoPairBlocks.size());
    }
    
    /**
     Determines size of selected binned GTOs pair block.

     @param iGtoPairBlock the index of GTOs pair block.
     @return the number of contracted GTOs pairs in selected binned GTOs pair block.
    */
    auto
    size(const int32_t iGtoPairBlock) const -> int32_t
    {
        return _gtoPairBlocks[iGtoPairBlock].getNumberOfContrPairs();
    }
    
    /**
     Gets selected binned GTOs pair block.

     @param iGtoPairBlock the index of GTOs pair block.
     @return the binned GTOs pair  block.
    */
    auto
    getGtoPairBlock(const int32_t iGtoPairBlock) const -> CBinnedGtoPairBlock<T, B>
    {
        return _gtoPairBlocks[iGtoPairBlock];
    }
    
    /**
     Gets bra angular momentum of packed GTOs pair block object.

     @return the bra angular momentum of packed GTOs pair block object.
    */
    auto
    getBraAngularMomentum() const -> int32_t
    {
        for (const auto& gpblock : _gtoPairBlocks)
        {
            if (!gpblock.empty())
            {
                return gpblock.getBraAngularMomentum();
            }
        }
            
        return -1;
    }
    
    /**
     Gets ket angular momentum of packed GTOs pair block object.

     @return the ket angular momentum of packed GTOs pair block object.
    */
    auto
    getKetAngularMomentum() const -> int32_t
    {
        for (const auto& gpblock : _gtoPairBlocks)
        {
            if (!gpblock.empty())
            {
                return gpblock.getKetAngularMomentum();
            }
        }
            
        return -1;
    }
    
    /**
     Gets number of primitive Gaussian pairs in contracted pair.

     @return the number of primitive Gaussian pairs.
     */
    auto
    getNumberOfPrimGtoPairs() const -> int32_t
    {
        for (const auto& gpblock : _gtoPairBlocks)
        {
            if (!gpblock.empty())
            {
                return gpblock.getNumberOfPrimPairs();
            }
        }
            
        return -1;
    }
    
    /**
     Gets number of subblocks of GTOs pairs devided by given numner of elements.

     @param iGtoPairBlock the index of GTOs pair block.
     @param nElements the number of GTOs pairs in subblock.
     @return the number of subblocks.
     */
    auto
    getNumberOfSubBlocks(const int32_t iGtoPairBlock,
                         const int32_t nElements) const -> int32_t
    {
        if (const auto nblocks = _gtoPairBlocks[iGtoPairBlock].getNumberOfContrPairs(); (nblocks % nElements) != 0)
        {
            return nblocks / nElements + 1;
        }
        else
        {
            return nblocks / nElements;
        }
    }

    /**
     Prints summary of binned GTOs blocks to output string.

     @return the output string.
     */
    auto
    printSummary() const -> std::string
    {
        std::stringstream ss;

        ss << "Packed GTOs Pair Block:\n";

        ss << std::string(23, '=') << "\n\n";

        ss << "Angular Momentum: |";
        
        ss << fstr::to_AngularMomentum(getBraAngularMomentum());
        
        ss << fstr::to_AngularMomentum(getKetAngularMomentum()) << ">\n";
        
        ss << "Number of Contracted Pairs: ";
        
        ss << getNumberOfPrimGtoPairs() << "\n\n";
        
        for (size_t i = 0; i < _gtoPairBlocks.size(); i++)
        {
            const auto ncpairs = _gtoPairBlocks[i].getNumberOfContrPairs();
            
            if (ncpairs > 0)
            {
                ss << "Index = " << i << " No. Pairs = ";
            
                ss << ncpairs << "\n";
            }
        }
        
        return ss.str();
    }
};

template <typename T, typename B = mem::Host>
using VPackedGtoPairBlocks = std::vector<CPackedGtoPairBlock<T, B>>;

#endif /* PackedGtoPairBlock_hpp */
