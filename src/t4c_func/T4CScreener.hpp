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

#ifndef T4CScreener_hpp
#define T4CScreener_hpp

#include <cstdint>
#include <string>
#include <vector>

#include "BlockedGtoPairBlock.hpp"
#include "GtoPairBlock.hpp"

/// @brief Class CT4CScreener provides methods for storing blocked GTOs pair blocks
/// partitioned according to Cauchy–Schwarz relationship.
class CT4CScreener
{
   public:
    /// Creates an four center integrals screener.
    CT4CScreener() = default;

    /// @brief The creates an four center integral screener.
    /// @param gto_pair_blocks The vector of basis function pair blocks.
    CT4CScreener(const std::vector<CBlockedGtoPairBlock> &gto_pair_blocks);

    /// @brief The default copy constructor.
    /// @param other The integrals screener to be copied.
    CT4CScreener(const CT4CScreener &other);

    /// @brief The default move constructor.
    /// @param other The integrals screener  to be moved.
    CT4CScreener(CT4CScreener &&other) noexcept;

    /// @brief The default destructor.
    ~CT4CScreener() = default;

    /// @brief The default copy assignment operator.
    /// @param other The integrals screener to be copy assigned.
    /// @return The assigned integrals screener.
    auto operator=(const CT4CScreener &other) -> CT4CScreener &;

    /// @brief The default move assignment operator.
    /// @param other The integrals screener to be move assigned.
    /// @return The assigned integrals screener .
    auto operator=(CT4CScreener &&other) noexcept -> CT4CScreener &;

    /// @brief The equality operator.
    /// @param other The integrals screener  to be compared.
    /// @return True if integrals screeners  are equal, False otherwise.
    auto operator==(const CT4CScreener &other) const -> bool;

    /// @brief The non-equality operator.
    /// @param other The integrals screener to be compared.
    /// @return True if integrals screeners  are not equal, False otherwise.
    auto operator!=(const CT4CScreener &other) const -> bool;

    /// @brief Partitions basis function pair blocks for given molecule and molecular basis for given
    /// type of four center integral.
    /// @param basis  The molecular basis.
    /// @param molecule  The molecule.
    /// @param label  The label of four center integral.
    auto partition(const CMolecularBasis &basis, const CMolecule &molecule, const std::string &label) -> void;

    auto partition_atom(const CMolecularBasis& basis, const CMolecule& molecule, const std::string& label, const int iatom) -> void;

    auto partition_atom_pair(const CMolecularBasis& basis, const CMolecule& molecule, const std::string& label, const int iatom, const int jatom) -> void;

    /// @brief Gets vector of blocked basis function pair blocks.
    /// @return The vector of  blocked basis function pair blocks.
    auto gto_pair_blocks() const -> std::vector<CBlockedGtoPairBlock>;

   private:
    /// @brief Vector of partitioned basis function pair blocks.
    std::vector<CBlockedGtoPairBlock> _gto_pair_blocks;
};

#endif /* T4CScreener_hpp */
