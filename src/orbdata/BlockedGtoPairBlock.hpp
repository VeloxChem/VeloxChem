//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef BlockedGtoPairBlock_hpp
#define BlockedGtoPairBlock_hpp

#include <array>
#include <vector>

#include "GtoPairBlock.hpp"

/// @brief Class CBlockedGtoPairBlock stores data about packed basis function pairs and provides set of
/// methods for manipulating with packed basis function pairs.
class CBlockedGtoPairBlock
{
   public:
    /// Creates an empty  blocked basis function pairs block.
    CBlockedGtoPairBlock() = default;
    
    /// @brief Creates a blocked basis function pairs block.
    /// @param gto_pair_blocks  The array of  basis function pairs blocks.
    CBlockedGtoPairBlock(const std::array<CGtoPairBlock, 16>& gto_pair_blocks);

    /// @brief Creates a blocked basis function pairs block.
    /// @param gto_pair_blocks  The vector of  basis function pairs blocks.
    /// @param block_indices  The vector of indexes to distribute the vector of blocked basis function pairs blocks.
    CBlockedGtoPairBlock(const std::vector<CGtoPairBlock> &gto_pair_blocks, const std::vector<int> &block_indices);

    /// Creates a blocked basis function pairs block.
    /// @param gto_pair_block  The basis function pairs block.
    /// @param integrals The vector of integrals to partition  basis function pairs block.
    CBlockedGtoPairBlock(const CGtoPairBlock &gto_pair_block, const std::vector<double> &integrals);

    /// @brief The default copy constructor.
    /// @param other The blocked basis function pairs block to be copied.
    CBlockedGtoPairBlock(const CBlockedGtoPairBlock &other);

    /// @brief The default move constructor.
    /// @param other The blocked basis functions pairs block to be moved.
    CBlockedGtoPairBlock(CBlockedGtoPairBlock &&other) noexcept;

    /// @brief The default destructor.
    ~CBlockedGtoPairBlock() = default;

    /// @brief The default copy assignment operator.
    /// @param other The blocked basis function pairs block to be copy assigned.
    /// @return The assigned basis function pairs block.
    auto operator=(const CBlockedGtoPairBlock &other) -> CBlockedGtoPairBlock &;

    /// @brief The default move assignment operator.
    /// @param other The blocked basis function pairs block to be move assigned.
    /// @return The assigned basis function pairs block.
    auto operator=(CBlockedGtoPairBlock &&other) noexcept -> CBlockedGtoPairBlock &;

    /// @brief The equality operator.
    /// @param other The blocked basis function pairs block to be compared.
    /// @return True if blocked basis function pairs blocks are equal, False otherwise.
    auto operator==(const CBlockedGtoPairBlock &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The blocked basis function pairs  block to be compared.
    /// @return True if blocked basis function pairs blocks are not equal, False
    /// otherwise.
    auto operator!=(const CBlockedGtoPairBlock &other) const -> bool;

    /// @brief Gets selected basis function pairs block.
    /// @param index  The index of basis function pairs block.
    /// @return The selected basis function pairs blocks.
    auto gto_pair_block(const int index) const -> CGtoPairBlock;

    /// @brief Checks if selected basis function pairs block is empty.
    /// @param index  The index of basis function pairs block.
    /// @return True if selected basis function pairs block is empty, False otherwise.
    auto is_empty_gto_pair_block(const int index) const -> bool;
    
    /// @brief Gets array of basis function pairs blocks.
    /// @return The array of basis function pairs blocks.
    auto gto_pair_blocks() const -> std::array<CGtoPairBlock, 16>;
    
    /// @brief Gets number of unique terms (integral bra or ket side of integral) generated
    /// by selected basis function pairs block.
    /// @param index  The index of basis function pairs block.
    /// @return The number of unique terms.
    auto unique_terms(const int index) const -> size_t;

   private:
    /// @brief The array of basis function pairs blocks.
    std::array<CGtoPairBlock, 16> _gto_pair_blocks;
};

#endif /* BlockedGtoPairBlock_hpp */
